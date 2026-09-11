#!/usr/bin/env python3
"""Collect one representative GenBank per GCF, for an exact centre-to-centre tree.

Under partitioning the merged `distance` table holds only within-partition
comparisons, so a global tree over all BGCs would have to substitute a constant
for every cross-partition pair — 43% of the matrix on Erwiniaceae, which leaves
the deep topology arbitrary.

A tree over *family centres* avoids that entirely: run BiG-SCAPE again on just
the centres and every centre pair is measured. The set is small and grows with
diversity rather than with BGC count — 19 centres for Erwiniaceae, 81 for
Streptomyces, 100 combined — so this stays cheap where the all-BGCs comparison
does not.

    python scripts/clustering/extract_family_centers.py \\
        --db merged.db --outdir centre_input --cutoff 0.30
"""
import argparse
import shutil
import sqlite3
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--db', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--cutoff', type=float, default=0.30)
    ap.add_argument('--manifest', type=Path, default=None,
                    help='TSV mapping family id -> centre record and GBK')
    args = ap.parse_args()

    con = sqlite3.connect(f'file:{args.db}?mode=ro', uri=True)
    rows = con.execute(
        """SELECT f.id, f.center_id, g.path, COUNT(rf.record_id) AS n
           FROM family f
           JOIN bgc_record_family rf ON rf.family_id = f.id
           JOIN bgc_record br ON br.id = rf.record_id
           JOIN bgc_record c ON c.id = f.center_id
           JOIN gbk g ON g.id = c.gbk_id
           WHERE f.cutoff = ? AND br.record_type = 'region'
           GROUP BY f.id ORDER BY n DESC""", (args.cutoff,)).fetchall()
    con.close()

    if not rows:
        raise SystemExit(f'no families at cutoff {args.cutoff}')

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for fam_id, centre_id, path, n in rows:
        src = Path(path)
        if not src.exists():
            print(f'  centre GBK missing for family {fam_id}: {src}', file=sys.stderr)
            continue
        # One directory per centre so BiG-SCAPE's recursive glob sees the same
        # layout it gets from antiSMASH, and so names cannot collide.
        d = args.outdir / f'family_{fam_id:05d}'
        d.mkdir(exist_ok=True)
        shutil.copy2(src, d / src.name)
        manifest.append((fam_id, centre_id, n, src.name))

    print(f'{len(manifest)} family centres from {len(rows)} families '
          f'(cutoff {args.cutoff})')
    if args.manifest:
        with args.manifest.open('w') as fh:
            fh.write('family_id\tcenter_record_id\tn_members\tgbk\n')
            for fam_id, centre_id, n, name in manifest:
                fh.write(f'{fam_id}\t{centre_id}\t{n}\t{name}\n')
    # Two centres are the minimum for a distance; one family means no tree.
    if len(manifest) < 3:
        print('fewer than 3 centres — a centre tree is not meaningful',
              file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
