#!/usr/bin/env python3
"""Predict the phosphonate intermediate each GCF makes: 2-AEP, 2-HEP, or neither.

2-AEP and 2-HEP are the two common intermediates in phosphonate biosynthesis, and
which one a cluster makes is decided by a single enzyme. The first two steps are
shared -- pepM gives phosphonopyruvate, Ppd decarboxylates it to phosphonoacetaldehyde
-- and then:

    PnAA + an aepZ-family transaminase  ->  2-AEP   (2-aminoethylphosphonate)
    PnAA + a  PnAA reductase            ->  2-HEP   (2-hydroxyethylphosphonate)

Evidence is TIERED, and the tier is part of the call rather than a footnote --
"2-AEP by homology to aepZ" and "2-AEP route available on domain evidence" are
different claims and must not read alike:

    homologue        sequence hit to a characterised reference at >= 40% identity
    weak homologue   the same at >= 25%
    class V transaminase / Fe-ADH with Ppd
                     domain-family evidence only, no usable sequence hit
    unknown          Ppd present but the third enzyme is absent, or is a class I/II
                     transaminase, which is NOT the aepZ family
    not via PnAA     no Ppd, so neither route is open

Sequence alone was tried first and called 18 of 19 families "none". The reference
set is simply too thin: aepZ is the only characterised 2-AEP transaminase available,
so a genuine orthologue in a distant lineage sits under any usable floor. Four
families carry Ppd with a class-V transaminase and none hits aepZ at 25%/60%
coverage. Reporting those as "none" was precision masquerading as knowledge.

Domain evidence alone is not enough either, which is why both are kept: PF00266
covers the whole class-V PLP family and a cluster can carry the domain without the
function. The references are:

    2-AEP   aepZ from B. fragilis NCTC 9343 PS B locus (AF285774), whose product
            carries 2-AEP on O-4 of the GlcNAc residue of the PS B repeat unit
    2-HEP   phosphonoacetaldehyde reductase from the B. fragilis 2-HEP locus, and
            the Fe-ADH from Glycomyces sp. NRRL B-16210 (KJ125437), whose
            phosphonoglycan carries 2-HEP confirmed by 31P NMR and MS

A worked negative: neither lab-characterised Erwiniaceae cluster hits aepZ at
E<1e-3. GCF-18's transaminase is class I/II (PF00155), not the class V aepZ family.
The confirmed phosphonolipid therefore does NOT use the canonical 2-AEP machinery.
It is reported as "unknown (Ppd + class I/II transaminase, not aepZ family)" rather
than a bare "unknown", because a bare unknown invites the reader to assume 2-AEP --
which is the assumption the evidence actually rules out.

Carrier class is reported separately and only when a branch point is called, because
the intermediate and its carrier are independent: 2-AEP and 2-HEP both appear on glycans, and
the same intermediate can go onto a lipid instead by host machinery the cluster does
not encode.
"""
import argparse
import collections
import csv
import json
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.domain_functions import category  # noqa: E402
from utils.antismash_parser import genome_dir_map  # noqa: E402

# Identity floors. Deliberately low: these are the only characterised references for
# either route and they come from two phyla, so a genuine orthologue in a third
# lineage can sit well under 40%. A hit between the floors is reported as `possible`
# rather than promoted, and the identity is always published beside the call.
# Profile bitscores, not percent identities. Calls were identical at 30, 50, 80
# and 150 across 367 regions -- the lowest true hit scored 77.4 against medians
# of 194-720 -- so these are floors, not tuned values. STRONG is where a hit is
# unambiguous; POSSIBLE admits a distant orthologue and is reported as such.
STRONG_BITS, POSSIBLE_BITS = 120.0, 50.0


def build_profiles(refs, workdir, hmmbuild, hmmalign):
    """One profile HMM per reference class, built at run time from the FASTA.

    Profiles rather than pairwise identity, and the difference is not cosmetic.
    The DIAMOND version this replaced called 18 of 19 families "none", and its
    own note blamed a reference set too thin for a distant orthologue to clear
    any usable floor. That diagnosis was half right: the set was thin, but the
    method was the larger problem. Measured on the same references over 367
    regions, pairwise identity found 2 of 26 Enterobacterial AEP clusters and a
    profile found all 26. A profile scores what the family conserves; identity
    scores only how close you are to one sequence.

    Built by seed-and-refine, which is why this repository needs no
    multiple-alignment tool: hmmbuild from the longest single member (a lone
    sequence is a trivial alignment), hmmalign the rest of the class to that
    model, hmmbuild again from the resulting alignment.
    """
    by_class = collections.defaultdict(list)
    for header, seq in read_fasta(refs):
        by_class[header.split('|')[0]].append((header.split('|')[1], seq))

    hmms = []
    for cls, members in sorted(by_class.items()):
        faa = workdir / f'{cls}.faa'
        faa.write_text(''.join(f'>{n}\n{s}\n' for n, s in members))
        seed_name, seed_seq = max(members, key=lambda m: len(m[1]))
        seed_faa = workdir / f'{cls}_seed.faa'
        seed_faa.write_text(f'>{seed_name}\n{seed_seq}\n')
        seed_hmm = workdir / f'{cls}_seed.hmm'
        _run([hmmbuild, '--amino', '-n', f'{cls}_seed', str(seed_hmm), str(seed_faa)])
        if len(members) == 1:
            hmms.append((cls, seed_hmm))
            continue
        sto = workdir / f'{cls}.sto'
        with sto.open('w') as fh:
            _run([hmmalign, '--amino', '--trim', str(seed_hmm), str(faa)], stdout=fh)
        hmm = workdir / f'{cls}.hmm'
        _run([hmmbuild, '--amino', '-n', cls, str(hmm), str(sto)])
        hmms.append((cls, hmm))
    return hmms


def read_fasta(path):
    """Line-based on purpose: splitting on '>' breaks a header containing '->'."""
    recs, cur = [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            cur = [line[1:], []]
            recs.append(cur)
        elif cur is not None:
            cur[1].append(line.strip())
    return [(h, ''.join(b)) for h, b in recs]


def _run(cmd, stdout=None):
    proc = subprocess.run(cmd, stdout=stdout or subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.exit(f'{Path(cmd[0]).name} failed:\n{proc.stderr[:2000]}')


def run_hmmsearch(refs, query, threads, workdir, hmmbuild, hmmalign, hmmsearch):
    """{query name: {class: (bitscore, class, evalue)}} -- one search per class.

    The return shape matches what the DIAMOND version produced so the call logic
    below is untouched, except that the score is now a profile bitscore rather
    than a percent identity.
    """
    best = collections.defaultdict(dict)
    for cls, hmm in build_profiles(refs, workdir, hmmbuild, hmmalign):
        tbl = workdir / f'{cls}.tbl'
        _run([hmmsearch, '--cpu', str(threads), '--tblout', str(tbl),
              str(hmm), str(query)])
        for line in tbl.read_text().splitlines():
            if line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            name, evalue, score = cols[0], float(cols[4]), float(cols[5])
            cur = best[name].get(cls)
            if cur is None or score > cur[0]:
                best[name][cls] = (score, cls, evalue)
    return best


def load_regions(antismash_paths, db_path, cutoff):
    dirs = genome_dir_map(antismash_paths)
    from Bio import SeqIO
    db = sqlite3.connect(db_path)
    members = collections.defaultdict(list)
    for path, fid in db.execute(
            """SELECT g.path, bf.family_id FROM bgc_record_family bf
               JOIN bgc_record r ON r.id=bf.record_id JOIN gbk g ON g.id=r.gbk_id
               JOIN family f ON f.id=bf.family_id WHERE f.cutoff=?""", (cutoff,)):
        p = Path(path)
        members[fid].append((p.parent.name, p.name))
    out = {}
    for fid, ms in members.items():
        regions = []
        for genome, fname in ms:
            base = dirs.get(genome)
            fp = (base / fname) if base else None
            if not fp or not fp.exists():
                continue
            label, doms, genes = fname[:-4], collections.defaultdict(list), []
            for rec in SeqIO.parse(str(fp), 'genbank'):
                for f in rec.features:
                    if f.type == 'PFAM_domain':
                        t = (f.qualifiers.get('locus_tag') or ['?'])[0]
                        for x in f.qualifiers.get('db_xref', []):
                            if x.startswith('PF'):
                                doms[t].append(x.split('.')[0])
            for rec in SeqIO.parse(str(fp), 'genbank'):
                for f in rec.features:
                    if f.type != 'CDS':
                        continue
                    q = f.qualifiers
                    t = (q.get('locus_tag') or q.get('gene') or ['?'])[0]
                    seq = (q.get('translation') or [''])[0]
                    if seq:
                        genes.append(dict(tag=t, seq=seq, accessions=doms.get(t, []),
                                          product=(q.get('product') or [''])[0]))
            regions.append(dict(label=label, genome=genome, genes=genes))
        if regions:
            out[fid] = regions
    return out


# Route enzymes by domain family. Weaker evidence than a sequence hit, but not
# negligible, and the sequence references are thin: aepZ is the ONLY characterised
# 2-AEP transaminase available, so a genuine orthologue in a distant lineage can sit
# below any usable identity floor. Measured here: four families carry Ppd plus a
# class-V transaminase and none of them hits aepZ at 25%/60% coverage. Reporting
# those as "none" would have been precision masquerading as knowledge.
PPD_DOMS = {'PF02775', 'PF02776'}
AEP_TRANSAMINASE = {'PF00266'}   # Aminotran_5, class V PLP — the family aepZ is in
OTHER_TRANSAMINASE = {'PF00155'}  # Aminotran_1_2, class I/II — NOT the aepZ family
HEP_REDUCTASE = {'PF00465', 'PF25137'}

CARRIER_DOMS = {'PF00534', 'PF00535', 'PF13439', 'PF13579', 'PF00953'}
CARRIER_TXT = re.compile(r'glycosyltransferase|glycosyl transferase|flippase|'
                         r'polysaccharide|nucleotide sugar dehydrogenase', re.I)


def carrier_class(region):
    """How many carrier-building genes sit in this region."""
    n = 0
    for g in region['genes']:
        if CARRIER_DOMS & set(g['accessions']) or CARRIER_TXT.search(g['product'] or ''):
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--antismash', type=Path, required=True, nargs='+',
                    help='antiSMASH results: either the taxon directory or the '
                         'genome directories themselves, as Nextflow stages them')
    ap.add_argument('--db', type=Path, required=True)
    ap.add_argument('--references', type=Path,
                    default=Path(__file__).resolve().parents[2] /
                    'assets' / 'reference_sequences' / 'reference_branch_point_enzymes.faa')
    ap.add_argument('--cutoff', type=float, default=0.3)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--hmmbuild', default='hmmbuild')
    ap.add_argument('--hmmalign', default='hmmalign')
    ap.add_argument('--hmmsearch', default='hmmsearch')
    ap.add_argument('--threads', type=int, default=4)
    a = ap.parse_args()

    fams = load_regions(a.antismash, a.db, a.cutoff)
    work = a.out.parent / '_branch_point'
    work.mkdir(parents=True, exist_ok=True)
    q = work / 'query.faa'
    index = {}
    with q.open('w') as fh:
        for fid, regions in fams.items():
            for r in regions:
                for g in r['genes']:
                    key = f'{fid}::{r["label"]}::{g["tag"]}'
                    index[key] = (fid, r['label'], g)
                    fh.write(f'>{key}\n{g["seq"]}\n')
    best = run_hmmsearch(a.references, q, a.threads, work,
                         a.hmmbuild, a.hmmalign, a.hmmsearch)

    rows = []
    for fid, regions in sorted(fams.items()):
        # per-region call, then the family's majority
        calls, evid = [], []
        for r in regions:
            hits = collections.defaultdict(list)
            for g in r['genes']:
                b = best.get(f'{fid}::{r["label"]}::{g["tag"]}', {})
                for cls, (pid, ref, ev) in b.items():
                    hits[cls].append((pid, g['tag'], ref))
            has_ppd = 'PPD' in hits
            accs = set()
            for g in r['genes']:
                accs |= set(g['accessions'])
            has_ppd = has_ppd or bool(accs & PPD_DOMS)

            # Tiered, most specific first. The tier is part of the call, not a
            # footnote: "2-AEP by homology to aepZ" and "2-AEP route available on
            # domain evidence" are different claims and must not read alike.
            call, pct, ref = 'not via PnAA', None, ''
            if has_ppd:
                cand = ([(p, 'AEP', rf) for p, _, rf in hits.get('AEP', [])] +
                        [(p, 'HEP', rf) for p, _, rf in hits.get('HEP', [])])
                cand = [c for c in cand if c[0] >= POSSIBLE_BITS]
                if cand:
                    p, cls, rf = max(cand)
                    tier = 'homologue' if p >= STRONG_BITS else 'weak homologue'
                    call = f'2-{cls} ({tier})'
                    # A profile has no single best-matching sequence -- the model
                    # is the whole class -- so the reference column names the
                    # profile, not a member. Reading it as one reference would
                    # overstate what a profile hit tells you.
                    pct, ref = p, f'{cls} profile'
                elif accs & AEP_TRANSAMINASE:
                    call = '2-AEP (class V transaminase, no profile hit)'
                elif accs & HEP_REDUCTASE:
                    call = '2-HEP (Fe-ADH with Ppd, no profile hit)'
                elif accs & OTHER_TRANSAMINASE:
                    # This branch used to fire on LMG 5342 region 2 and call it
                    # "not aepZ family". That was wrong, and instructively so: the
                    # region carries TWO transaminases, a deposited class I/II
                    # aspartate transaminase AND a recovered aepZ-family
                    # 2-aminoethylphosphonate--pyruvate transaminase. Reading the
                    # deposited annotation alone saw only the first, so the
                    # conclusion was an artefact of the very annotation gap
                    # RECOVER_ORFS exists to close. Laboratory work confirms that
                    # cluster runs the 2-AEP pathway. The branch is kept for a
                    # region that genuinely has only a class I/II transaminase,
                    # but it no longer asserts the aepZ family is absent -- it
                    # cannot know that from what it can see.
                    call = 'unknown (Ppd + class I/II transaminase only)'
                else:
                    call = 'unknown (Ppd, no third enzyme found)'
            calls.append(call)
            evid.append((pct, ref))
        top, n = collections.Counter(calls).most_common(1)[0]
        pcts = [p for p, _ in evid if p is not None]
        refs = collections.Counter(r for _, r in evid if r)
        rows.append(dict(
            gcf=fid, members=len(regions), branch_point=top,
            support=f'{n}/{len(regions)}',
            best_bitscore=round(max(pcts), 1) if pcts else '',
            reference=refs.most_common(1)[0][0] if refs else '',
            carrier_genes=round(sum(carrier_class(r) for r in regions) / len(regions), 1),
            carrier_class=('glycan' if sum(carrier_class(r) for r in regions) / len(regions) >= 4
                           else 'none in cluster')))

    cols = ['gcf', 'members', 'branch_point', 'support', 'best_bitscore', 'reference',
            'carrier_genes', 'carrier_class']
    with a.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    for f in work.glob('*'):
        f.unlink()
    try:
        work.rmdir()
    except OSError:
        pass
    print(f'wrote {a.out}')
    c = collections.Counter(r['branch_point'] for r in rows)
    print('branch point calls: ' + ', '.join(f'{k} {v}' for k, v in c.most_common()))


if __name__ == '__main__':
    main()
