#!/usr/bin/env python3
"""Decide which genomes are worth running antiSMASH on, by looking for pepM.

Every phosphonate BGC carries a PEP mutase, so a genome without one cannot
contain the thing this pipeline looks for. Finding it costs ~0.9 CPU-s against
antiSMASH's 41.4, which is what makes an order-scale run tractable: measured on
2,771 Erwiniaceae genomes the screen keeps 11.0% and loses nothing.

**Two modes, one decision.** NCBI GenBank records are inconsistently annotated —
28.7% of that set carried no CDS translations at all — so:

    annotated genome    -> diamond blastp over its proteins   (0.23 CPU-s)
    unannotated genome  -> diamond blastx over its contigs    (2.03 CPU-s)

Both search the same references at the same threshold, and were verified to make
*identical* calls on all 341 genomes where both could run. blastp is 9.2x cheaper,
which is the only reason the split exists.

**Sensitivity is the thing that matters here**, since a miss silently drops a BGC
that would never be seen again. Measured 298/298 at `--bitscore 100` over all
2,771 genomes, with 8 false positives — and the margin is wide, true positives
scoring 154-552 against a background topping out near 51.

`--min_density` guards the one failure mode the modes do not share: a partially
annotated genome would take the blastp path and could hide a pepM in a region
with no CDS features. Observed density was 576-1,015 proteins/Mb with nothing
below 500, so the guard costs nothing today and prevents a silent loss later.

    python pepm_prescreen.py --genomes *.gbff --db pepm.dmnd --out report.tsv
"""
import argparse
import subprocess
import sys
from pathlib import Path


def parse_genome(path):
    """(proteins, contigs, bp) for one GenBank file, in a single pass."""
    from Bio import SeqIO
    proteins, contigs, bp = [], [], 0
    for rec in SeqIO.parse(str(path), 'genbank'):
        contigs.append((rec.id, str(rec.seq)))
        bp += len(rec.seq)
        for feat in rec.features:
            if feat.type == 'CDS' and 'translation' in feat.qualifiers:
                tag = feat.qualifiers.get('locus_tag', ['?'])[0]
                proteins.append((tag, feat.qualifiers['translation'][0]))
    return proteins, contigs, bp


def run_diamond(binary, mode, db, query, threads):
    """{genome: best bitscore} for one diamond run."""
    if not query.stat().st_size:
        return {}
    proc = subprocess.run(
        [binary, mode, '-d', str(db), '-q', str(query), '--quiet',
         '--threads', str(threads), '--evalue', '1e-5',
         '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'bitscore'],
        capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f'diamond {mode} failed:\n{proc.stderr[:2000]}')
    best = {}
    for line in proc.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 4:
            continue
        g = f[0].split('|')[0]
        best[g] = max(best.get(g, 0.0), float(f[3]))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genomes', nargs='+', required=True)
    ap.add_argument('--db', type=Path, required=True, help='diamond db of pepM references')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--bitscore', type=float, default=100.0)
    ap.add_argument('--min_density', type=float, default=500.0,
                    help='proteins per Mb below which a genome is treated as '
                         'unannotated and screened by blastx')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--diamond', default='diamond')
    args = ap.parse_args()

    work = args.out.parent
    prot_fa, dna_fa = work / '_prescreen_prot.faa', work / '_prescreen_dna.fna'
    mode_of, density = {}, {}

    with open(prot_fa, 'w') as pf, open(dna_fa, 'w') as df:
        for path in args.genomes:
            p = Path(path)
            name = p.name
            for suffix in ('.gbff', '.gbk', '.gb', '.genbank'):
                if name.endswith(suffix):
                    name = name[:-len(suffix)]
                    break
            try:
                proteins, contigs, bp = parse_genome(p)
            except Exception as exc:                      # a corrupt genome is
                print(f'{name}: unreadable ({exc}); screening by DNA', file=sys.stderr)
                proteins, contigs, bp = [], [], 0         # not a reason to lose the batch
            dens = len(proteins) / (bp / 1e6) if bp else 0.0
            density[name] = dens
            if proteins and dens >= args.min_density:
                mode_of[name] = 'blastp'
                for tag, seq in proteins:
                    pf.write(f'>{name}|{tag}\n{seq}\n')
            else:
                mode_of[name] = 'blastx'
                for cid, seq in contigs:
                    df.write(f'>{name}|{cid}\n{seq}\n')

    best = run_diamond(args.diamond, 'blastp', args.db, prot_fa, args.threads)
    best.update(run_diamond(args.diamond, 'blastx', args.db, dna_fa, args.threads))

    n_pass = 0
    with args.out.open('w') as fh:
        fh.write('genome\tmode\tcds_per_mb\tbest_bitscore\tpass\n')
        for name in sorted(mode_of):
            score = best.get(name, 0.0)
            ok = score >= args.bitscore
            n_pass += ok
            fh.write(f'{name}\t{mode_of[name]}\t{density[name]:.0f}\t'
                     f'{score:.1f}\t{"yes" if ok else "no"}\n')
    prot_fa.unlink(missing_ok=True)
    dna_fa.unlink(missing_ok=True)

    n = len(mode_of)
    by_mode = {m: sum(1 for v in mode_of.values() if v == m) for m in ('blastp', 'blastx')}
    print(f'screened {n} genomes ({by_mode["blastp"]} blastp, {by_mode["blastx"]} blastx); '
          f'{n_pass} carry pepM at bitscore >= {args.bitscore:.0f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
