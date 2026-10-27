#!/usr/bin/env python3
"""Decide how to split a taxon's metadata fetch, and print the shards.

Why this exists: `datasets` cannot hold one stream open long enough to deliver a
whole-domain manifest. Measured on Bacteria GenBank (3,350,527 records), four
consecutive attempts died with `stream error ... INTERNAL_ERROR; received from
peer` at 0.29 / 1.30 / 0.53 / 1.31 GB of an expected ~1.72 GB -- up to 76% in,
after one to two hours each, with a retry restarting from zero.

Smaller fetches are fine: 149,539 records (Bacteroidota) completes in ~4.5 min.
So the fix is to ask for less at a time.

SPLIT BY RELEASE DATE, not by sub-taxon. Date ranges need no taxonomy database
and no new tool in this process's environment, they apply uniformly to any taxon
rather than depending on which ranks exist beneath it, and they subdivide as far
as needed. Verified on Pantoea that `--released-before X` and `--released-after
X` are an EXACT partition: 310 + 1,497 = 1,807, with zero accessions shared and
zero missing.

BELOW THE THRESHOLD NOTHING CHANGES. One shard is emitted with no date flags, so
a small taxon runs the identical single call it runs today. That matters because
an extra `datasets` invocation costs ~0.9 s of fixed overhead (measured: 8
records 0.9 s, 345 records 2.3 s, 1,807 records 7.9 s -- 0.9 s + 3.9 ms/record),
so sharding a genus unconditionally would make it ~4x slower for no benefit.
Same shape as `bigscape_partition_threshold`, which is off below 10,000 BGCs
because partitioning measured 2-3.7x SLOWER at 185-518.
"""
import argparse
import json
import sys
import urllib.parse
import urllib.request
from datetime import date

API = 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon'


def count(taxon, source, atypical, released_after=None, timeout=120):
    """Records matching these filters, from the REST report's total_count.

    `released_after` is the number released ON OR AFTER that date, so the count
    in [a, b) is count(a) - count(b). The filter name is `first_release_date`,
    established by probing: it returns 1,497 for Pantoea at 2020-01-01, which is
    exactly what the CLI's `--released-after 2020-01-01` returns. Four other
    plausible spellings were tried and every one FAILED OPEN -- the filter was
    ignored and the full 1,807 came back. That failure mode is why this is
    pinned by test rather than by guess: a silently ignored filter would size
    every shard as the whole taxon and defeat the sharding entirely.

    Deliberately NOT `datasets summary --limit 1`: that caps total_count at the
    limit and returns 1, which would make every taxon look tiny. Verified
    against the CLI -- Pantoea 1,807, Bacteroidota 149,539, Bacteria 3,350,527.
    """
    q = {'page_size': 1, 'filters.assembly_source': source.lower()}
    if atypical:
        q['filters.exclude_atypical'] = 'true'
    if released_after:
        q['filters.first_release_date'] = f'{released_after}T00:00:00Z'
    url = f'{API}/{urllib.parse.quote(taxon)}/dataset_report?' + urllib.parse.urlencode(q)
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return int(json.load(fh).get('total_count', 0))


def split(cnt, total, threshold, start_year=1995):
    """Date boundaries chosen from MEASURED counts, not from assumed uniformity.

    Assemblies skew hard toward recent years, so equal-width year bins would put
    most of the taxon in the last few shards -- the exact thing being avoided.
    Each boundary count costs ~1 s, so ~30 of them is a rounding error against a
    fetch measured in hours, and it is worth paying to know rather than assume.

    A year that exceeds the threshold on its own is subdivided by month. If a
    single MONTH still exceeds it, that is accepted and reported: the threshold
    sits far below where fetches actually start failing (149,539 observed good
    against 3,350,527 bad), so overshooting it is not by itself a failure.
    """
    years = list(range(start_year, date.today().year + 2))
    marks = [(f'{y}-01-01', cnt(f'{y}-01-01')) for y in years]

    def emit(cuts):
        """N cut points make N+1 ranges, open at both ends.

        An earlier version returned N ranges by pairing each cut with the next,
        which collapsed a single cut into one unbounded shard -- i.e. no split at
        all, silently. The sum check in main() is what caught it.
        """
        out, prev = [], None
        for c in cuts:
            out.append((prev, c))
            prev = c
        out.append((prev, None))
        return out

    bounds, carried = [], 0
    for i, (d, remaining) in enumerate(marks):
        in_bin = remaining - (marks[i + 1][1] if i + 1 < len(marks) else 0)
        if in_bin > threshold:                             # this year alone is too big
            if carried:
                bounds.append(d); carried = 0
            for m in range(1, 13):
                md = f'{d[:4]}-{m:02d}-01'
                if md not in bounds:
                    bounds.append(md)
            continue
        carried += in_bin
        if carried > threshold:
            bounds.append(d); carried = in_bin
    bounds = sorted(set(bounds))
    return emit(bounds) if bounds else [(None, None)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--taxon', required=True)
    ap.add_argument('--assembly_source', required=True)
    ap.add_argument('--threshold', type=int, required=True)
    ap.add_argument('--exclude_atypical', action='store_true')
    a = ap.parse_args()

    if a.threshold <= 0:                        # sharding disabled outright
        print('-\t-')
        return

    try:
        total = count(a.taxon, a.assembly_source, a.exclude_atypical)
    except Exception as exc:                    # noqa: BLE001
        # Fail OPEN to today's behaviour: one unsharded call. The count is an
        # optimisation, and a network blip deciding it must not turn into a
        # changed fetch strategy.
        print(f'WARNING: could not count {a.taxon} ({exc}); fetching unsharded',
              file=sys.stderr)
        print('-\t-')
        return

    if total <= a.threshold:
        print(f'{a.taxon}: {total:,} records, at or under the {a.threshold:,} '
              f'threshold -- one fetch, unchanged', file=sys.stderr)
        print('-\t-')
        return

    cnt = lambda d: count(a.taxon, a.assembly_source, a.exclude_atypical, released_after=d)
    shards = split(cnt, total, a.threshold)
    print(f'{a.taxon}: {total:,} records over the {a.threshold:,} threshold -- '
          f'{len(shards)} date shards', file=sys.stderr)
    # Report the implied sizes. A shard set that does not sum to the total means
    # the boundaries are wrong, and losing genomes silently is the failure this
    # whole change exists to avoid.
    sizes = [(cnt(af) if af else total) - (cnt(bf) if bf else 0) for af, bf in shards]
    print(f'  shard sizes: min {min(sizes):,} max {max(sizes):,} sum {sum(sizes):,} '
          f'({"EXACT" if sum(sizes) == total else "MISMATCH vs %d" % total})',
          file=sys.stderr)
    # "-" for an open end, never an empty field. Tab is an IFS *whitespace*
    # character, so `IFS=$'\t' read -r A B` STRIPS a leading tab instead of
    # yielding an empty first field -- which silently turned the first shard's
    # --released-before into a --released-after and lost every record before the
    # first cut. Caught by the sharded-vs-unsharded accession comparison.
    for after, before in shards:
        print(f'{after or "-"}\t{before or "-"}')


if __name__ == '__main__':
    main()
