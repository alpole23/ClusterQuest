#!/usr/bin/env python3
"""Compare two BiG-SCAPE databases family-for-family, at every shared cutoff.

Written for the 8 kb vs 10 kb neighbourhood question, but it takes any two databases
over the same regions. Records are matched by region filename, not by record id --
ids are autoincrements and mean nothing across runs, the same reason family ids do
(see "GCF numbering stays BiG-SCAPE's").

Reports split and merged separately rather than only an index, because they are
different problems: a merge joins chemistry the other run kept apart, a split breaks
a family you may already have written about.
"""
import argparse
import collections
import os
import sqlite3
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'scripts'))
from viz.cutoff_comparison import agreement, shape          # noqa: E402


def memberships(db_path):
    """{cutoff: {region filename: family_id}} -- keyed by name so ids need not match."""
    db = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        out = collections.defaultdict(dict)
        for cut, path, fam in db.execute(
                """SELECT f.cutoff, g.path, bf.family_id
                     FROM bgc_record_family bf
                     JOIN family f ON f.id = bf.family_id
                     JOIN bgc_record r ON r.id = bf.record_id
                     JOIN gbk g ON g.id = r.gbk_id
                    WHERE r.record_type = 'region'"""):
            out[cut][os.path.basename(path)] = fam
    finally:
        db.close()
    return {c: out[c] for c in sorted(out)}


def movers(a, b):
    """Members whose family differs, after matching families by best overlap.

    The metric bigscape_reproducibility/ settled on: counting changed co-membership
    PAIRS inflates one moving member into ~200 pairs.
    """
    ct = collections.Counter((a[r], b[r]) for r in a if r in b)
    best = {}
    for (x, y), _ in sorted(ct.items(), key=lambda kv: -kv[1]):
        best.setdefault(x, y)
    return sum(1 for r in a if r in b and best.get(a[r]) != b[r])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('new', help='the database to assess (e.g. the 8 kb run)')
    ap.add_argument('ref', help='the database to compare against (e.g. the 10 kb run)')
    args = ap.parse_args()

    new, ref = memberships(args.new), memberships(args.ref)
    shared_regions = set.intersection(
        *[set(m) for m in list(new.values()) + list(ref.values())])
    print(f'{len(shared_regions):,} regions present in both runs\n')

    hdr = (f'{"cutoff":>7} {"families":>18} {"largest":>14} {"singletons":>14}'
           f' {"ARI":>8} {"split":>9} {"merged":>9} {"moved":>7}')
    print(hdr)
    print('-' * len(hdr))
    for cut in sorted(set(new) & set(ref)):
        a = {r: f for r, f in ref[cut].items() if r in shared_regions}
        b = {r: f for r, f in new[cut].items() if r in shared_regions}
        sa, sb = shape(a), shape(b)
        ag = agreement(a, b) or {'ari': float('nan'), 'split': 0, 'merged': 0}
        print(f'{cut:>7.2f} {sa["families"]:>8} -> {sb["families"]:<7}'
              f' {sa["largest"]:>6} -> {sb["largest"]:<5}'
              f' {sa["singletons"]:>6} -> {sb["singletons"]:<5}'
              f' {ag["ari"]:>8.4f} {ag["split"]:>9,} {ag["merged"]:>9,}'
              f' {movers(a, b):>7}')
    print('\nleft of each arrow is the reference run, right is the new one')


if __name__ == '__main__':
    main()
