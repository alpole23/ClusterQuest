#!/usr/bin/env python3
"""Does partitioning BiG-SCAPE by pepM identity rebuild the same GCFs?

`pepm_all_by_all.py` reports how many same-GCF pairs a pepM cut would separate.
That is a necessary check but not a sufficient one: it asks whether the cut
*could* split a family, using the unpartitioned run's own family labels. It
cannot tell you what BiG-SCAPE actually does when re-run on smaller inputs.

The difference matters because family assignment is a global step over the
distance matrix, not a per-pair threshold. Removing BGCs changes the clustering
context, so a partitioned run can in principle produce different families even
though every distance *within* a partition is identical — BiG-SCAPE recomputes
them from the same domain content either way.

So this actually does it: split the region GBKs by pepM component, run BiG-SCAPE
on each, and compare the resulting families to the unpartitioned run.

Families are compared by **co-membership, not by label**. GCF ids are per-run
SQLite autoincrements (see the GCF numbering note in CLAUDE.md), so "family 3"
means nothing across runs. The question is only ever whether BGCs x and y landed
together in both.

    python scripts/bench_bigscape_partitioned.py \\
        --db results/bigscape_results/Streptomyces/Streptomyces.db \\
        --pairs <pepm_vs_neighbourhood.tsv> --threshold 0.60 \\
        --bigscape <path> --pfam <Pfam-A.hmm> --workdir results/bench_partition
"""
import argparse
import collections
import csv
import itertools
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
from pathlib import Path


def load_reference(db):
    """{record_id: family_id} and {record_id: gbk path} from the whole-set run."""
    con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    fam, paths = {}, {}
    for r in con.execute(
            """SELECT br.id AS rid, g.path AS path, brf.family_id AS fid
               FROM bgc_record br
               JOIN gbk g ON g.id = br.gbk_id
               LEFT JOIN bgc_record_family brf ON brf.record_id = br.id
               WHERE br.record_type = 'region'"""):
        paths[r['rid']] = r['path']
        if r['fid'] is not None:
            fam[r['rid']] = r['fid']
    con.close()
    return fam, paths


def components(pairs_tsv, threshold, nodes):
    """Single-linkage components of pepM identity at the threshold."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    with open(pairs_tsv) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if float(row['pepm_identity']) >= threshold:
                a, b = int(row['record_a']), int(row['record_b'])
                if a in parent and b in parent:
                    ra, rb = find(a), find(b)
                    if ra != rb:
                        parent[ra] = rb
    groups = collections.defaultdict(list)
    for n in nodes:
        groups[find(n)].append(n)
    return sorted(groups.values(), key=len, reverse=True)


def run_partition(members, paths, workdir, bigscape, pfam, cores, cutoff, tag, gbk_root=None):
    """Run BiG-SCAPE on one partition; return {record_id: local family id}.

    A partition of one BGC has nothing to compare against, so BiG-SCAPE is not
    invoked — a lone BGC is a singleton family by definition, and running it
    would only pay the fixed Pfam-loading cost.
    """
    if len(members) < 2:
        return {members[0]: f'{tag}_solo'}, 0.0

    indir = workdir / f'in_{tag}'
    outdir = workdir / f'out_{tag}'
    for d in (indir, outdir):
        shutil.rmtree(d, ignore_errors=True)
    # Mirror the antiSMASH layout (one directory per genome) so BiG-SCAPE's
    # recursive glob and its "cluster"/"region" filename filter behave as they
    # do in the pipeline.
    for rid in members:
        src = Path(paths[rid])
        # A BiG-SCAPE database records gbk.path as it was at clustering time,
        # which for a pipeline run points into work/. That directory is routinely
        # cleaned (`nextflow clean`, or just reclaiming disk), so a database can
        # outlive the files it names. --gbk_root re-resolves by the layout every
        # arrangement of these files shares: <genome>/<region>.gbk.
        if not src.exists() and gbk_root:
            alt = Path(gbk_root) / src.parent.name / src.name
            if not alt.exists():
                found = list(Path(gbk_root).glob(f'*/{src.name}'))
                alt = found[0] if found else alt
            if not alt.exists():
                raise SystemExit(
                    f'{src} is gone and not found under {gbk_root}. '
                    f'Point --gbk_root at the antismash_results directory for this run.')
            src = alt
        gdir = indir / src.parent.name
        gdir.mkdir(parents=True, exist_ok=True)
        dest = gdir / src.name
        if not dest.exists():
            shutil.copy2(src, dest)

    env = dict(os.environ)
    env['PATH'] = f"{Path(bigscape).resolve().parent}{os.pathsep}{env.get('PATH', '')}"
    t0 = time.monotonic()
    proc = subprocess.run(
        [bigscape, 'cluster', '-i', str(indir), '-o', str(outdir),
         '--pfam-path', str(pfam), '--alignment-mode', 'auto',
         '--gcf-cutoffs', cutoff, '--include-singletons',
         '--cores', str(cores), '--classify', 'category'],
        capture_output=True, text=True, env=env)
    elapsed = time.monotonic() - t0
    if proc.returncode != 0:
        tail = (proc.stderr or '').strip().splitlines()
        raise SystemExit(f'bigscape failed on partition {tag}: '
                         f'{tail[-1] if tail else proc.returncode}')

    db = next(outdir.rglob('*.db'), None)
    if db is None:
        raise SystemExit(f'no database produced for partition {tag}')
    con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
    # Map back by gbk path basename: record ids are per-run, paths are not.
    byname = {Path(paths[r]).name: r for r in members}
    out = {}
    for name, fid in con.execute(
            """SELECT g.path, brf.family_id FROM bgc_record br
               JOIN gbk g ON g.id = br.gbk_id
               JOIN bgc_record_family brf ON brf.record_id = br.id
               WHERE br.record_type = 'region'"""):
        rid = byname.get(Path(name).name)
        if rid is not None:
            out[rid] = f'{tag}_{fid}'
    con.close()
    shutil.rmtree(indir, ignore_errors=True)
    return out, elapsed


def compare(ref, test):
    """Adjusted Rand Index plus the pairs the two clusterings disagree on.

    ARI is the right summary because it scores agreement on *co-membership* and
    corrects for the agreement two random clusterings would show by chance. The
    raw pair counts are reported alongside it because they say which direction
    an error goes: splitting a real family is a different failure from merging
    two distinct ones.
    """
    shared = sorted(set(ref) & set(test))
    n = len(shared)
    together_ref = together_test = both = 0
    split, merged = [], []
    for a, b in itertools.combinations(shared, 2):
        r = ref[a] == ref[b]
        t = test[a] == test[b]
        together_ref += r
        together_test += t
        both += r and t
        if r and not t:
            split.append((a, b))
        elif t and not r:
            merged.append((a, b))

    total = n * (n - 1) // 2
    if total == 0:
        return {}
    expected = together_ref * together_test / total
    max_idx = (together_ref + together_test) / 2
    ari = ((both - expected) / (max_idx - expected)) if max_idx != expected else 1.0
    return {
        'bgcs_compared': n,
        'pairs': total,
        'together_reference': together_ref,
        'together_partitioned': together_test,
        'agreed_together': both,
        'split_by_partitioning': len(split),
        'merged_by_partitioning': len(merged),
        'adjusted_rand_index': round(ari, 6),
        'identical': not split and not merged,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--db', type=Path, required=True, help='unpartitioned BiG-SCAPE db')
    ap.add_argument('--pairs', type=Path, required=True,
                    help='pepm_vs_neighbourhood.tsv from pepm_all_by_all.py')
    ap.add_argument('--pfam', type=Path, required=True)
    ap.add_argument('--bigscape', default='bigscape')
    ap.add_argument('--threshold', type=float, default=0.60)
    ap.add_argument('--cores', type=int, default=8)
    ap.add_argument('--cutoff', default='0.30')
    ap.add_argument('--workdir', type=Path, required=True)
    ap.add_argument('--gbk_root', type=Path,
                    help='antismash_results directory to re-resolve region GBKs '
                         'against, when the paths recorded in the database no '
                         'longer exist (work/ cleaned since the run)')
    ap.add_argument('--keep', action='store_true', help='keep per-partition outputs')
    args = ap.parse_args()

    args.workdir.mkdir(parents=True, exist_ok=True)
    ref_fam, paths = load_reference(args.db)
    parts = components(args.pairs, args.threshold, sorted(paths))
    print(f'reference: {len(paths)} BGCs, {len(set(ref_fam.values()))} families')
    print(f'partitions at pepM identity {args.threshold}: {len(parts)}, '
          f'sizes {[len(p) for p in parts][:12]}'
          f'{" ..." if len(parts) > 12 else ""}')

    test_fam, total_s = {}, 0.0
    for i, members in enumerate(parts):
        fam, secs = run_partition(members, paths, args.workdir, args.bigscape,
                                  args.pfam, args.cores, args.cutoff, f'p{i:03d}',
                                  gbk_root=args.gbk_root)
        test_fam.update(fam)
        total_s += secs
        print(f'  partition {i:>3}: {len(members):>4} BGCs -> '
              f'{len(set(fam.values())):>3} families  ({secs:>6.1f}s)')
        if not args.keep:
            shutil.rmtree(args.workdir / f'out_p{i:03d}', ignore_errors=True)

    result = compare(ref_fam, test_fam)
    result.update({'threshold': args.threshold, 'partitions': len(parts),
                   'largest_partition': len(parts[0]),
                   'families_reference': len(set(ref_fam.values())),
                   'families_partitioned': len(set(test_fam.values())),
                   'bigscape_seconds_partitioned': round(total_s, 1)})
    (args.workdir / f'partition_check_{args.threshold}.json').write_text(
        json.dumps(result, indent=2))

    print()
    print(f"  families      reference {result['families_reference']:>4}   "
          f"partitioned {result['families_partitioned']:>4}")
    print(f"  co-membership pairs     reference {result['together_reference']:>6,}   "
          f"partitioned {result['together_partitioned']:>6,}")
    print(f"  split by partitioning   {result['split_by_partitioning']:>6,}"
          "   (real families broken apart — the failure that matters)")
    print(f"  merged by partitioning  {result['merged_by_partitioning']:>6,}"
          "   (distinct families joined)")
    print(f"  adjusted Rand index     {result['adjusted_rand_index']:>6.4f}"
          f"   {'IDENTICAL' if result['identical'] else ''}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
