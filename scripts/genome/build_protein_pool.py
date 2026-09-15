#!/usr/bin/env python3
"""Pool every annotated protein across a set of genomes, for ORF recovery.

This is the reference set the homology half of `recover_orfs.py` searches. The
premise is that annotation quality is uneven between submissions: P. ananatis
LMG 5342's 2012 deposit never called the 2-AEP transaminase in its phosphonolipid
cluster, while P. ananatis VY148's 2021 assembly calls it and 8 other genes in the
same cluster. A gene absent from one genome is usually present in a sibling.

Headers are `>locus_tag|product`, which is what recover_orfs.py splits to name a
recovered gene rather than leaving it "hypothetical protein".

Proteins are deduplicated by sequence. Within a screened set most genomes are close
relatives, so the same protein appears many times over; collapsing them cuts the
DIAMOND database substantially without losing any reference.
"""
import argparse
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genomes', nargs='+', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--min_aa', type=int, default=30,
                    help='ignore very short proteins; they generate spurious blastx hits')
    a = ap.parse_args()

    from Bio import SeqIO
    seen = {}
    genomes = skipped = total = 0
    for path in a.genomes:
        try:
            recs = list(SeqIO.parse(str(path), 'genbank'))
        except Exception as err:                      # a malformed deposit is not fatal
            print(f'  WARN: {path.name}: {err}', file=sys.stderr)
            skipped += 1
            continue
        genomes += 1
        for rec in recs:
            for f in rec.features:
                if f.type != 'CDS':
                    continue
                q = f.qualifiers
                seq = (q.get('translation') or [''])[0]
                if len(seq) < a.min_aa:
                    continue
                total += 1
                if seq in seen:
                    continue
                tag = (q.get('locus_tag') or q.get('gene') or ['?'])[0]
                product = (q.get('product') or [''])[0].strip() or 'hypothetical protein'
                seen[seq] = f'{tag}|{product}'

    with a.out.open('w') as fh:
        for seq, header in seen.items():
            fh.write(f'>{header}\n{seq}\n')
    print(f'pooled {len(seen):,} unique proteins from {total:,} CDS across {genomes} genomes'
          + (f' ({skipped} unreadable)' if skipped else ''))


if __name__ == '__main__':
    main()
