#!/usr/bin/env python3
"""Pool every annotated protein across a set of genomes, for ORF recovery.

This is the reference set the homology half of `recover_orfs.py` searches. The
premise is that annotation quality is uneven between submissions: P. ananatis
LMG 5342's 2012 deposit never called the 2-AEP transaminase in its phosphonolipid
cluster, while P. ananatis VY148's 2021 assembly calls it and 8 other genes in the
same cluster. A gene absent from one genome is usually present in a sibling.

Headers are `>locus_tag|product`, which is what recover_orfs.py splits to name a
recovered gene rather than leaving it "hypothetical protein".

Proteins are deduplicated by sequence -- tracked as a sha256 digest and streamed
straight to the output, so the pool is never resident in memory. Within a screened
set most genomes are close relatives, so the same protein appears many times over;
collapsing them cuts the DIAMOND database substantially without losing any reference.
"""
import argparse
import hashlib
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
    # Remember proteins by DIGEST, and write each one out the moment it is first
    # seen. Both halves are needed and the second is the one that does the work:
    # keying on a hash saves nothing if the sequence still has to be held until a
    # final write pass, so the sequence is streamed to the file and then dropped.
    #
    # This used to key the dict on the sequence itself, which kept the whole pool
    # resident: measured 10.2 GB over 19,818,123 proteins on the whole-Bacteria
    # sweep, ~550 bytes each, linear in the pool. That is the single thing that
    # stopped this machine going wider -- ~30,000 retained genomes implies ~142 M
    # proteins and ~78 GB, against 56 GB of RAM. Measured here at 516 -> 108
    # bytes per protein, 4.8x, which takes that case to ~14 GB.
    #
    # Full sha256, not a truncation. A collision drops a protein SILENTLY, so the
    # failure mode argues for the wide digest; at 1e8 proteins the probability
    # against 2^256 is far below a cosmic-ray bit flip, and hashing costs ~1 us
    # per protein against a stage that runs for hours.
    #
    # Output order is unchanged: the old code relied on dict insertion order,
    # which IS first-seen order, and so is writing as we go. First header wins
    # either way. The pool is byte-identical, which is what makes this testable.
    seen = set()
    genomes = skipped = total = 0
    with a.out.open('w') as fh:
        for path in a.genomes:
            try:
                recs = list(SeqIO.parse(str(path), 'genbank'))
            except Exception as err:                  # a malformed deposit is not fatal
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
                    digest = hashlib.sha256(seq.encode()).digest()
                    if digest in seen:
                        continue
                    seen.add(digest)
                    tag = (q.get('locus_tag') or q.get('gene') or ['?'])[0]
                    product = (q.get('product') or [''])[0].strip() or 'hypothetical protein'
                    fh.write(f'>{tag}|{product}\n{seq}\n')

    print(f'pooled {len(seen):,} unique proteins from {total:,} CDS across {genomes} genomes'
          + (f' ({skipped} unreadable)' if skipped else ''))


if __name__ == '__main__':
    main()
