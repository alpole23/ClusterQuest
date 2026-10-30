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

# Deliberately self-contained. An earlier version imported agreement() and shape()
# from scripts/viz/cutoff_comparison.py, which lives on an unmerged branch -- so this
# script raised ModuleNotFoundError on main and would have failed at exactly the
# moment it was needed. A comparison record has to run from a clean checkout.


def _comb2(n):
    return n * (n - 1) // 2


def shape(member):
    """Family-count statistics for one {record: family} map."""
    sizes = collections.Counter(member.values())
    vals = sorted(sizes.values())
    return {
        'regions': len(member),
        'families': len(sizes),
        'largest': vals[-1] if vals else 0,
        'singletons': sum(1 for v in vals if v == 1),
        'median': vals[len(vals) // 2] if vals else 0,
    }


def agreement(a, b):
    """Co-membership change from clustering `a` to clustering `b`.

    `split` is pairs together in A and apart in B; `merged` is the reverse. ARI is
    the standard contingency-table form -- 1,476 regions is 1.09 million pairs and
    the closed form is O(n + k1*k2).
    """
    shared = sorted(set(a) & set(b))
    if len(shared) < 2:
        return None
    ct = collections.Counter((a[r], b[r]) for r in shared)
    ra = collections.Counter(a[r] for r in shared)
    rb = collections.Counter(b[r] for r in shared)
    n = len(shared)
    both = sum(_comb2(v) for v in ct.values())
    sa = sum(_comb2(v) for v in ra.values())
    sb = sum(_comb2(v) for v in rb.values())
    total = _comb2(n)
    expected = sa * sb / total if total else 0.0
    denom = 0.5 * (sa + sb) - expected
    return {
        'n': n,
        'split': sa - both,
        'merged': sb - both,
        'ari': (both - expected) / denom if denom else 1.0,
    }


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
