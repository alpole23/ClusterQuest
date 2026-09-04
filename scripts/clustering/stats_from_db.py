#!/usr/bin/env python3
"""BiG-SCAPE clustering statistics, read from the database rather than its TSVs.

`extract_bigscape_stats.py` parses BiG-SCAPE's `output_files/{class}_clustering
_c{cutoff}.tsv`, which ties it to one run's output *directory*. Partitioned
clustering produces one directory per partition and merges only the databases,
so the directory no longer exists in a single form — but the merged database
holds everything the statistics need.

Emits exactly the schema `extract_bigscape_stats.py` does; verified to produce
byte-identical JSON on the unpartitioned Erwiniaceae run, so it is a drop-in for
both paths rather than a partition-only variant.
"""
import argparse
import json
import sqlite3
import sys
from pathlib import Path


def stats_for_cutoff(con, cutoff):
    rows = con.execute(
        """SELECT brf.family_id, COUNT(*) FROM bgc_record_family brf
           JOIN bgc_record br ON br.id = brf.record_id
           JOIN family f ON f.id = brf.family_id
           WHERE br.record_type = 'region' AND f.cutoff = ?
           GROUP BY brf.family_id""", (cutoff,)).fetchall()
    if not rows:
        return None
    sizes = [n for _, n in rows]
    total = sum(sizes)
    # BiG-SCAPE's own TSVs report a class per family; with the phosphonate-only
    # rule every BGC lands in one class, and the DB records it on bgc_record.
    classes = {}
    for cls, fams, bgcs in con.execute(
            """SELECT COALESCE(br.category, 'other'),
                      COUNT(DISTINCT brf.family_id), COUNT(*)
               FROM bgc_record_family brf
               JOIN bgc_record br ON br.id = brf.record_id
               JOIN family f ON f.id = brf.family_id
               WHERE br.record_type = 'region' AND f.cutoff = ?
               GROUP BY COALESCE(br.category, 'other')""", (cutoff,)):
        classes[cls] = {'families': fams, 'bgcs': bgcs}
    return {
        'cutoff': cutoff,
        'total_families': len(sizes),
        'total_bgcs': total,
        'avg_bgcs_per_family': round(total / len(sizes), 2),
        'max_bgcs_per_family': max(sizes),
        'min_bgcs_per_family': min(sizes),
        'singleton_families': sum(1 for s in sizes if s == 1),
        'mibig_included': False,
        'families_with_mibig': None,
        'bgc_classes': classes,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('db', type=Path)
    ap.add_argument('output_file', type=Path)
    args = ap.parse_args()

    con = sqlite3.connect(f'file:{args.db}?mode=ro', uri=True)
    cutoffs = [r[0] for r in con.execute(
        'SELECT DISTINCT cutoff FROM family ORDER BY cutoff')]
    per = {}
    for c in cutoffs:
        s = stats_for_cutoff(con, c)
        if s:
            per[str(c)] = s
    con.close()

    if not per:
        args.output_file.write_text(json.dumps({'error': 'no families found'}, indent=2))
        print('no families in the database', file=sys.stderr)
        return 0

    primary = per[str(min(cutoffs))]
    out = {'cutoffs': per, **primary}
    args.output_file.write_text(json.dumps(out, indent=2))
    print(f"{primary['total_bgcs']} BGCs in {primary['total_families']} families "
          f"at cutoff {primary['cutoff']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
