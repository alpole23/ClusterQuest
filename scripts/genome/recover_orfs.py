#!/usr/bin/env python3
"""Recover genes that a genome's submitted annotation left out, as GFF3 for antiSMASH.

Partial annotation is common and silently destructive. Measured on the Erwiniaceae
run: 327 of 333 BGC regions contain an intergenic gap of 300+ bp, the median region is
78.4% coding where a well-annotated operon is ~90%, and 914 kb of sequence inside BGC
regions carries no gene call at all. The worst case is *Pantoea ananatis* LMG 5342
region 2 -- a lab-confirmed phosphonolipid cluster -- where 9 of 15 genes were never
called, including the 2-AEP transaminase that decides the product's headgroup.

antiSMASH will not fix this itself. `common/record_processing.py` runs gene finding
only when a record has ZERO CDS features:

    if not sequence.get_cds_features():
        if not options.genefinding_gff3 and options.genefinding_tool != "none":
            genefinding(sequence, options)

One CDS anywhere on a 4.6 Mb chromosome and prodigal never runs. But `--genefinding-gff3`
is checked on the same line and its features are merged in `parse_input_sequence` with
`record.features.extend(...)`, *before* that gate -- so GFF3 is a supported way to add
gene calls to an already-annotated record. No patched antiSMASH required.

TWO RECOVERY METHODS, because neither is sufficient alone:

  prodigal   finds ORFs de novo, including genes no relative has annotated. On
             LMG 5342 it recovered 8 of the 9 missing genes at 99.6-100% identity.
  homology   blastx of the unannotated gaps against a pooled protein set. This is
             what catches what prodigal fragments -- on LMG 5342 prodigal truncated
             aepZ to 162 aa (GTG start, score 3.22, conf 67.7) where the real gene is
             238 aa, and that is precisely the gene the whole exercise is about.

Where the two disagree on the same stretch of DNA, the homology call wins: it is
anchored to a real protein of known length, whereas prodigal's boundary is a
prediction. Neither is allowed to overlap a CDS the submitter already called --
this stage only ever fills gaps, never overrides existing annotation.
"""
import argparse
import bisect
import collections
import shutil
import subprocess
import sys
from pathlib import Path

MIN_GAP = 300            # smallest gap worth searching; a 100-codon protein needs ~300 bp
MIN_AA = 50              # shorter recovered ORFs are mostly noise
MIN_PCT = 40.0           # homology floor; these are within-run relatives, not distant refs
STOPS = {'TAA', 'TAG', 'TGA'}
STARTS = {'ATG', 'GTG', 'TTG'}


def gff_escape(text):
    """Percent-encode the characters GFF3 reserves inside an attribute value.

    Not cosmetic. An unescaped ';' in a note splits the attribute, so
    "recovered by homology; 100.0% identity to X" became a bogus qualifier named
    '100.0% identity to X' with the value "true", and the GenBank antiSMASH then
    wrote could not be parsed at all.
    """
    for ch, code in (('%', '%25'), (';', '%3B'), ('=', '%3D'),
                     (',', '%2C'), ('&', '%26'), ('\t', '%09')):
        text = text.replace(ch, code)
    return text


def existing_cds(rec):
    spans = sorted((int(f.location.start), int(f.location.end))
                   for f in rec.features if f.type == 'CDS')
    return spans, [s for s, _ in spans]


# Adjacent genes in a bacterial operon commonly overlap by a few bases (the ATGA
# motif overlaps stop and start by 4). Rejecting any overlap at all threw out the
# prodigal call for aepZ over an 11 bp overlap with the upstream Ppd.
MAX_OVERLAP = 45


def overlaps(spans, starts, s, e, tol=MAX_OVERLAP):
    i = bisect.bisect_right(starts, e)
    for s2, e2 in spans[max(0, i - 50):i]:
        ov = min(e, e2) - max(s, s2)
        if ov > tol:
            return True
    return False


def gaps_of(rec, spans):
    """Unannotated stretches of at least MIN_GAP, plus the record's two ends."""
    out = []
    if not spans:
        return [(0, len(rec.seq))]
    if spans[0][0] >= MIN_GAP:
        out.append((0, spans[0][0]))
    for (_, e1), (s2, _) in zip(spans, spans[1:]):
        if s2 - e1 >= MIN_GAP:
            out.append((max(0, e1 - MAX_OVERLAP), s2 + MAX_OVERLAP))
    if len(rec.seq) - spans[-1][1] >= MIN_GAP:
        out.append((spans[-1][1], len(rec.seq)))
    return out


def run_prodigal(binary, fasta, out_gff):
    proc = subprocess.run([binary, '-i', str(fasta), '-f', 'gff', '-o', str(out_gff), '-q'],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f'prodigal failed:\n{proc.stderr[:2000]}')
    calls = collections.defaultdict(list)
    for line in Path(out_gff).read_text().splitlines():
        if line.startswith('#'):
            continue
        p = line.split('\t')
        if len(p) < 8 or p[2] != 'CDS':
            continue
        calls[p[0]].append((int(p[3]) - 1, int(p[4]), p[6]))
    return calls


def extend_to_codons(seq, s, e, strand):
    """Grow a blastx HSP out to a start codon and a stop codon, in frame.

    A blastx hit covers the aligned part of a protein, not the gene. Without this the
    recovered feature is a fragment and antiSMASH's translation is wrong at both ends.
    Returns None when no stop codon is reachable, which usually means the gap is
    truncated at a contig edge.
    """
    n = len(seq)
    if strand == '+':
        while e + 3 <= n:
            if str(seq[e:e + 3]).upper() in STOPS:
                e += 3
                break
            e += 3
        else:
            return None
        best = s
        while s - 3 >= 0:
            cod = str(seq[s - 3:s]).upper()
            if cod in STOPS:
                break
            s -= 3
            if cod in STARTS:
                best = s
        return best, e
    while s - 3 >= 0:
        if str(seq[s - 3:s].reverse_complement()).upper() in STOPS:
            s -= 3
            break
        s -= 3
    else:
        return None
    best = e
    while e + 3 <= n:
        cod = str(seq[e:e + 3].reverse_complement()).upper()
        if cod in STOPS:
            break
        e += 3
        if cod in STARTS:
            best = e
    return s, best


def run_blastx(diamond, gaps_fasta, pool_db, threads):
    """One diamond call. Returns raw outfmt-6 lines.

    Separate from the HSP processing below so a BATCH can pay this once. Measured
    on 20 genomes against the 19.8 M-protein pool: 20 separate calls cost 2,536
    CPU-s, one call over the concatenated gaps cost 265 -- about 120 CPU-s of that
    is loading the 7.19 GB database, so 19 of 20 loads were pure overhead.
    """
    proc = subprocess.run(
        [diamond, 'blastx', '-d', str(pool_db), '-q', str(gaps_fasta), '--quiet',
         '--threads', str(threads), '--evalue', '1e-10', '--max-target-seqs', '25',
         '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'qstart', 'qend', 'bitscore'],
        capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f'diamond blastx failed:\n{proc.stderr[:2000]}')
    return proc.stdout.splitlines()


def orfs_from_hits(hit_lines, records):
    """Square blastx HSPs to codon boundaries and drop the ones already annotated."""
    # Every non-overlapping HSP per gap, not just the best one. A gap wide enough to be
    # worth searching is usually wide enough to hold more than one gene: the 1,278 bp
    # gap in LMG 5342 region 2 holds aepZ AND a hypothetical, and keeping one hit per
    # gap silently dropped whichever scored lower.
    per_gap = collections.defaultdict(list)
    for line in hit_lines:
        q, sub, pid, qs, qe, bits = line.split('\t')
        if float(pid) < MIN_PCT:
            continue
        per_gap[q].append((float(bits), sub, float(pid), int(qs), int(qe)))
    hsps = []
    for q, rows in per_gap.items():
        claimed = []
        for _, sub, pid, qs, qe in sorted(rows, key=lambda r: -r[0]):
            lo, hi = sorted((qs, qe))
            if any(lo < h and hi > l for l, h in claimed):
                continue
            claimed.append((lo, hi))
            hsps.append((q, sub, pid, qs, qe))
    out = collections.defaultdict(list)
    for q, sub, pid, qs, qe in hsps:
        rid, off = q.rsplit('|', 1)
        off = int(off.split('-')[0])
        strand = '+' if qs < qe else '-'
        lo, hi = sorted((qs, qe))
        s, e = off + lo - 1, off + hi
        ext = extend_to_codons(records[rid].seq, s, e, strand)
        if not ext:
            continue
        s, e = ext
        if (e - s) // 3 >= MIN_AA:
            out[rid].append((s, e, strand, sub, pid))
    return out


def homology_orfs(diamond, gaps_fasta, pool_db, threads, records):
    """Single-genome path: one blastx, then square its hits. Unchanged behaviour."""
    return orfs_from_hits(run_blastx(diamond, gaps_fasta, pool_db, threads), records)


def merge_and_write(records, spans, prod, hom, out_path, label):
    """Combine homology and prodigal calls into one GFF3.

    Shared by the single-genome and batch paths so the two cannot drift: the
    batch mode exists only to hoist the blastx out, and must produce byte-
    identical GFF3 otherwise.
    """
    # Merge. Homology first so it claims contested DNA: its boundaries come from a real
    # protein, prodigal's are a prediction, and prodigal is the one that truncates.
    added = collections.Counter()
    lines = ['##gff-version 3']
    for rid, rec in records.items():
        sp, st = spans[rid]
        taken = []
        picked = []
        for s, e, strand, sub, pid in sorted(hom.get(rid, [])):
            if overlaps(sp, st, s, e) or any(
                    min(e, e2) - max(s, s2) > MAX_OVERLAP for s2, e2 in taken):
                continue
            taken.append((s, e))
            picked.append((s, e, strand, 'homology', sub, pid))
            added['homology'] += 1
        for s, e, strand in sorted(prod.get(rid, [])):
            if (e - s) // 3 < MIN_AA:
                continue
            if overlaps(sp, st, s, e) or any(
                    min(e, e2) - max(s, s2) > MAX_OVERLAP for s2, e2 in taken):
                continue
            taken.append((s, e))
            picked.append((s, e, strand, 'prodigal', '', 0.0))
            added['prodigal'] += 1
        for i, (s, e, strand, src, sub, pid) in enumerate(sorted(picked), 1):
            tag = f'recovered_{rid.replace("|", "_")}_{i:04d}'
            # A homology call already knows what the gene is, so name it rather than
            # emitting "hypothetical protein" and making a later stage rediscover it.
            # The note keeps the provenance: this is an inference from a homologue.
            product = 'hypothetical protein'
            note = f'recovered by {src}'
            if src == 'homology':
                ref_tag, _, ref_product = sub.partition('|')
                if ref_product.strip():
                    product = ref_product.strip()
                note += f' - {pid:.1f}% identity to {ref_tag}'
            # CDS only, and childless: check_gff_suitability rejects a GFF where CDS
            # appears as a parent in the parent/child map.
            lines.append('\t'.join([rid, 'recover_orfs', 'CDS', str(s + 1), str(e),
                                    '.', strand, '0',
                                    f'ID={tag};locus_tag={tag};'
                                    f'product={gff_escape(product)};'
                                    f'note={gff_escape(note)}']))
    out_path.write_text('\n'.join(lines) + '\n')
    total = added['prodigal'] + added['homology']
    print(f'{label}: recovered {total} genes '
          f'({added["homology"]} by homology, {added["prodigal"]} by prodigal) -> {out_path}')


SEP = '@@'          # genome tag separator in a batched query id


def run_batch(a):
    """Recover ORFs for many genomes with ONE blastx over all their gaps.

    Why: diamond pays ~120 CPU-s loading the 7.19 GB pool database on every
    invocation, before searching anything. Per genome that dominates -- measured
    on 20 genomes, 20 calls cost 2,536 CPU-s against 265 for one, so 90% of the
    work was reloading the same database nineteen extra times.

    Three passes, not two, and the reason is memory: 200 parsed GenBanks is
    ~10 GB of sequence held at once. Pass 1 keeps only the small per-genome state
    (gene calls and CDS spans) and re-parses in pass 3, which costs ~0.5 s per
    genome against the ~110 CPU-s per genome this saves.

    A genome that fails stays a per-genome failure: it gets an empty GFF3 and the
    batch continues, as the shell loop this replaces did. Losing 200 genomes to
    one bad GenBank is the thing batching must not introduce.
    """
    from Bio import SeqIO
    a.workdir.mkdir(parents=True, exist_ok=True)
    gaps_path = a.workdir / 'all_gaps.fna'
    state, skipped = {}, []

    # ── pass 1: per genome, emit gaps and run prodigal ────────────────────
    with gaps_path.open('w') as gf:
        for g in a.genomes:
            stem = g.name[:-5] if g.name.endswith('.gbff') else g.stem
            # Both of these would mis-assign hits silently: a repeated stem makes two
            # genomes share an output path and a state slot, and a stem containing the
            # separator makes partition() split the tag in the wrong place.
            if SEP in stem:
                sys.exit(f'genome name contains {SEP!r}, which batching uses as a tag '
                         f'separator: {g.name}')
            if stem in state or stem in skipped:
                sys.exit(f'two genomes in one batch share the name {stem!r}')
            out_path = a.outdir / f'{stem}.gff3'
            try:
                records = {r.id: r for r in SeqIO.parse(str(g), 'genbank')}
                if not records:
                    raise ValueError('no records')
                n_cds = sum(1 for r in records.values()
                            for f in r.features if f.type == 'CDS')
                if n_cds == 0:
                    out_path.write_text('##gff-version 3\n')
                    print(f'{g.name}: no existing annotation, leaving it to antiSMASH '
                          f'(which runs its own gene finder on zero-CDS records)')
                    continue
                wd = a.workdir / f'_w_{stem}'
                wd.mkdir(parents=True, exist_ok=True)
                fasta = wd / 'genome.fna'
                SeqIO.write(list(records.values()), str(fasta), 'fasta')
                spans = {rid: existing_cds(r) for rid, r in records.items()}
                prod = run_prodigal(a.prodigal, fasta, wd / 'prodigal.gff')
                for rid, rec in records.items():
                    for st, en in gaps_of(rec, spans[rid][0]):
                        gf.write(f'>{stem}{SEP}{rid}|{st}-{en}\n{str(rec.seq[st:en])}\n')
                state[stem] = (g, spans, prod, out_path)
                # prod is parsed, so the FASTA and prodigal GFF are dead weight --
                # and unlike the per-genome loop this replaces, there is no
                # rm -rf between genomes to reclaim them.
                shutil.rmtree(wd, ignore_errors=True)
            # SystemExit too: run_prodigal and run_blastx exit rather than raise, and
            # one genome prodigal cannot call must not take the other 199 with it.
            except (Exception, SystemExit) as exc:        # noqa: BLE001
                print(f'  WARN: {g.name} failed in preparation ({exc}); empty GFF3')
                out_path.write_text('##gff-version 3\n')
                skipped.append(stem)
            finally:
                records = None

    if not state:
        # Same rule as every other batched process here: a batch exits non-zero
        # only when EVERY item failed, because that means the environment is
        # broken rather than the input. Reporting success with 200 empty GFF3s is
        # the silent failure this codebase keeps paying for.
        if skipped:
            sys.exit(f'every genome in the batch failed ({len(skipped)})')
        print('no genomes to recover')
        return

    # ── pass 2: ONE blastx over every gap in the batch ────────────────────
    by_stem = collections.defaultdict(list)
    db = a.pool_db
    if db is None and a.pool and a.pool.exists() and gaps_path.stat().st_size:
        # The pipeline always passes --pool_db. This is the standalone fallback, and
        # unlike the single-genome path it is paid once for the whole batch.
        db = a.workdir / 'pool'
        subprocess.run([a.diamond, 'makedb', '--in', str(a.pool), '-d', str(db),
                        '--quiet'], check=True, capture_output=True)
    if db and gaps_path.stat().st_size:
        for line in run_blastx(a.diamond, gaps_path, db, a.threads):
            stem, sep, rest = line.partition(SEP)
            if not sep:
                continue                                  # not ours; skip rather than mis-assign
            by_stem[stem].append(rest)
    hit_total = sum(len(v) for v in by_stem.values())
    print(f'  one blastx over {len(state)} genomes: {hit_total} hits')

    # A hit whose tag names no genome in this batch would be silently attached to
    # the wrong assembly downstream. Refuse instead.
    unknown = set(by_stem) - set(state)
    if unknown:
        sys.exit(f'blastx returned hits tagged with unknown genomes: {sorted(unknown)[:3]}')

    # ── pass 3: per genome, re-parse and merge ────────────────────────────
    # Wrapped for the same reason pass 1 is: every genome that got this far must
    # end up with a GFF3, or Nextflow's `path "*.gff3"` output is short one file
    # and the whole batch fails over a single bad record.
    for stem, (g, spans, prod, out_path) in state.items():
        try:
            records = {r.id: r for r in SeqIO.parse(str(g), 'genbank')}
            hom = orfs_from_hits(by_stem.get(stem, []), records)
            merge_and_write(records, spans, prod, hom, out_path, g.name)
        except (Exception, SystemExit) as exc:            # noqa: BLE001
            print(f'  WARN: {g.name} failed while merging ({exc}); empty GFF3')
            out_path.write_text('##gff-version 3\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genome', type=Path, help='GenBank to recover genes for')
    ap.add_argument('--out', type=Path, help='GFF3 to write')
    ap.add_argument('--genomes', type=Path, nargs='+',
                    help='batch mode: several GenBanks sharing ONE blastx. Requires '
                         '--outdir; writes <stem>.gff3 into it for each.')
    ap.add_argument('--outdir', type=Path, help='batch mode: directory for the GFF3s')
    ap.add_argument('--pool', type=Path, help='protein FASTA for homology recovery')
    ap.add_argument('--pool_db', type=Path,
                    help='prebuilt DIAMOND database for --pool. Built here if omitted, '
                         'which is correct but slow -- see the note in main().')
    ap.add_argument('--prodigal', default='prodigal')
    ap.add_argument('--diamond', default='diamond')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--workdir', type=Path, default=Path('.'))
    a = ap.parse_args()

    if a.genomes:
        if not a.outdir:
            ap.error('--genomes requires --outdir')
        a.outdir.mkdir(parents=True, exist_ok=True)
        run_batch(a)
        return
    if not (a.genome and a.out):
        ap.error('give --genome with --out, or --genomes with --outdir')

    from Bio import SeqIO
    records = {r.id: r for r in SeqIO.parse(str(a.genome), 'genbank')}
    if not records:
        sys.exit(f'no records in {a.genome}')

    # A genome with NO annotation at all needs nothing from us: that is precisely the
    # case antiSMASH handles itself, since it runs its own gene finder on records with
    # zero CDS features. Recovering them here duplicates that work at 20-130x the cost,
    # because the whole genome counts as a gap and the homology pass then scans ~4.9 Mb
    # against the full protein pool.
    #
    # Measured on the Erwiniaceae screened set: annotated genomes take 8-11 s each,
    # zero-CDS genomes 194-1,320 s. They are 154 of 307 genomes (50%) and would have
    # taken the stage from ~26 min to ~14.5 h.
    n_cds = sum(1 for r in records.values() for f in r.features if f.type == 'CDS')
    if n_cds == 0:
        a.out.write_text('##gff-version 3\n')
        print(f'{a.genome.name}: no existing annotation, leaving it to antiSMASH '
              f'(which runs its own gene finder on zero-CDS records)')
        return

    a.workdir.mkdir(parents=True, exist_ok=True)

    fasta = a.workdir / 'genome.fna'
    SeqIO.write(list(records.values()), str(fasta), 'fasta')

    spans = {rid: existing_cds(r) for rid, r in records.items()}
    prod = run_prodigal(a.prodigal, fasta, a.workdir / 'prodigal.gff')

    # gaps, for the homology pass
    hom = {}
    if a.pool and a.pool.exists():
        gf = a.workdir / 'gaps.fna'
        n = 0
        with gf.open('w') as fh:
            for rid, rec in records.items():
                for s, e in gaps_of(rec, spans[rid][0]):
                    fh.write(f'>{rid}|{s}-{e}\n{str(rec.seq[s:e])}\n')
                    n += 1
        if n:
            # The database is identical for every genome in a run -- the pool is one
            # gather step's single output -- so BUILD_PROTEIN_POOL builds it once and
            # hands it over. Building it here instead was 52% of this script's wall
            # time, the same 7.19 GB database rebuilt once per genome: 600.5 CPU-s per
            # genome against 126, measured on the whole-Bacteria pilot's 19.8 M-protein
            # pool. See docs/comparisons/recover_orfs_recalibration/.
            #
            # Note the gate above cannot save the work: `if n` only asks whether the
            # genome has ONE unannotated stretch of MIN_GAP bp, and every bacterial
            # genome has intergenic DNA at 300 bp. A pristine reference genome paid the
            # full build to search a handful of small gaps.
            #
            # The fallback stays so the script still runs standalone against a bare
            # FASTA, which is how it is used outside the pipeline.
            db = a.pool_db
            if db is None:
                db = a.workdir / 'pool'
                subprocess.run([a.diamond, 'makedb', '--in', str(a.pool), '-d', str(db),
                                '--quiet'], check=True, capture_output=True)
            hom = homology_orfs(a.diamond, gf, db, a.threads, records)

    merge_and_write(records, spans, prod, hom, a.out, a.genome.name)


if __name__ == '__main__':
    main()
