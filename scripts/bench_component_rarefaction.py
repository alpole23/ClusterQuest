#!/usr/bin/env python3
"""Does the largest pepM component grow with sampling depth, or saturate?

Partitioning BiG-SCAPE by pepM identity is verified lossless
(`bench_bigscape_partitioned.py`, ARI 1.0000 on three sets). What it does *not*
settle is how the partitioned run scales, because partitioning does not break
dense clusters apart — it dilutes them. Across every set measured, the Pantoea
component stayed exactly 236 BGCs while only its share moved: 71% alone, 46%
mixed with Streptomyces.

**Peak memory is set by the largest component in absolute terms, not its share.**
BiG-SCAPE memory fits `GB = 1.14 + 1.29e-7*n^2` on the largest job, so the
projection at a million genomes hinges entirely on whether a cluster like that
one keeps growing as more genomes are sampled, or stops.

This subsamples the observed set at increasing depths and reports how the
largest component and the component count accumulate — the same logic as
species rarefaction, applied to partition size.

**Random subsampling cannot answer the scaling question, and this script says
so.** Drawing a random subset from a pool where 45.6% of BGCs sit in one
component gives ~45.6% in that component at every depth — measured, slope ratio
0.95, share pinned within 44-47% across a 12x range. That is arithmetic, not
biology: the curve describes the pool it was drawn from and cannot see what
adding *new clades* would do, which is what growing to a million genomes
actually means.

Use `--add-from` for that. It holds one taxon fixed and adds increasing amounts
of another, which is the real experiment. Measured on Erwiniaceae + Streptomyces
the largest component stays at **exactly 236 BGCs** while 0 to 185 Streptomyces
BGCs are added — its share falls 70.9% to 45.6% purely because the denominator
grows, and the component count rises 6 to 24. New diversity adds components
beside the dense cluster; it never enlarges it.

The consequence for sizing: **peak memory is set by the most deeply sampled
single clade, not by total dataset size.** A million genomes spread across many
clades does not enlarge any one component; a million genomes concentrated on one
over-sequenced clade would.

    python scripts/bench_component_rarefaction.py \\
        --pairs <pepm_vs_neighbourhood.tsv> --threshold 0.60 \\
        --out results/bench_combined/component_rarefaction
"""
import argparse
import collections
import csv
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.plotting import SVG_METADATA, canonicalise_svg  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402


def load_edges(pairs_tsv, threshold):
    """Pairs joined at or above the pepM identity threshold, plus every node."""
    edges, nodes = [], set()
    with open(pairs_tsv) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            a, b = int(row['record_a']), int(row['record_b'])
            nodes.add(a)
            nodes.add(b)
            if float(row['pepm_identity']) >= threshold:
                edges.append((a, b))
    return sorted(nodes), edges


def components(sample, edges):
    """(largest, count) for the subgraph induced on `sample`."""
    parent = {n: n for n in sample}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    sizes = collections.Counter(find(n) for n in sample)
    return max(sizes.values()), len(sizes)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pairs', type=Path, required=True)
    ap.add_argument('--threshold', type=float, default=0.60)
    ap.add_argument('--replicates', type=int, default=30)
    ap.add_argument('--steps', type=int, default=12)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--add-from', type=Path, default=None,
                    help='BiG-SCAPE db of the combined run; switches to holding one '
                         'source taxon fixed and adding the other, which is the '
                         'experiment that bears on scaling')
    ap.add_argument('--base-tag', default='Erwiniaceae',
                    help='substring identifying the held-fixed taxon in gbk.path')
    args = ap.parse_args()

    nodes, edges = load_edges(args.pairs, args.threshold)
    n = len(nodes)
    print(f'{n} BGCs, {len(edges):,} pepM links at identity >= {args.threshold}')

    rng = random.Random(args.seed)

    if args.add_from:
        import sqlite3
        con = sqlite3.connect(f'file:{args.add_from}?mode=ro', uri=True)
        base, extra = [], []
        for rid, path in con.execute(
                "SELECT br.id, g.path FROM bgc_record br JOIN gbk g ON g.id = br.gbk_id "
                "WHERE br.record_type = 'region'"):
            (base if args.base_tag in path else extra).append(rid)
        con.close()
        print(f'holding {len(base)} {args.base_tag} BGCs fixed, '
              f'adding up to {len(extra)} others')
        print(f"  {'added':>8}{'total':>8}{'largest':>9}{'share':>8}{'components':>12}")
        curve = []
        for k in [round(len(extra) * i / 8) for i in range(9)]:
            big, cnt = [], []
            for _ in range(args.replicates if 0 < k < len(extra) else 1):
                s_ = set(base) | set(rng.sample(extra, k))
                L, C = components(s_, edges)
                big.append(L)
                cnt.append(C)
            L, C = sum(big) / len(big), sum(cnt) / len(cnt)
            curve.append({'added': k, 'total': len(base) + k, 'largest_mean': L,
                          'largest_share': L / (len(base) + k), 'components_mean': C})
            print(f'  {k:>8}{len(base)+k:>8}{L:>9.0f}{L/(len(base)+k):>8.1%}{C:>12.1f}')
        grew = curve[-1]['largest_mean'] - curve[0]['largest_mean']
        print(f'\n  largest component grew by {grew:+.0f} BGCs across the addition')
        print('  -> new diversity adds components beside the dense cluster, '
              'it does not enlarge it' if abs(grew) < 1 else
              '  -> the dense cluster absorbed some of the added BGCs')
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / 'component_addition.json').write_text(json.dumps(
            {'base_tag': args.base_tag, 'threshold': args.threshold,
             'curve': curve, 'largest_growth': grew}, indent=2))
        return 0

    depths = [max(2, round(n * (i + 1) / args.steps)) for i in range(args.steps)]
    rows = []
    for d in depths:
        big, cnt = [], []
        for _ in range(args.replicates if d < n else 1):
            s = set(rng.sample(nodes, d))
            L, C = components(s, edges)
            big.append(L)
            cnt.append(C)
        rows.append({
            'depth': d,
            'largest_mean': sum(big) / len(big),
            'largest_min': min(big), 'largest_max': max(big),
            'components_mean': sum(cnt) / len(cnt),
            'largest_share': (sum(big) / len(big)) / d,
        })
        print(f"  depth {d:>5}   largest {rows[-1]['largest_mean']:>7.1f}"
              f"   share {rows[-1]['largest_share']:>5.1%}"
              f"   components {rows[-1]['components_mean']:>6.1f}")

    # Is the largest component still growing linearly at full depth? Compare the
    # slope over the last third against the slope over the first third: a ratio
    # near 1 means unsaturated growth, well below 1 means it is flattening.
    third = max(2, len(rows) // 3)
    def slope(seg):
        dx = seg[-1]['depth'] - seg[0]['depth']
        return (seg[-1]['largest_mean'] - seg[0]['largest_mean']) / dx if dx else 0.0
    early, late = slope(rows[:third]), slope(rows[-third:])
    ratio = late / early if early else None

    summary = {'threshold': args.threshold, 'bgcs': n, 'replicates': args.replicates,
               'curve': rows, 'slope_early': early, 'slope_late': late,
               'slope_ratio': ratio}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'component_rarefaction.json').write_text(json.dumps(summary, indent=2))

    print(f'\n  growth of the largest component:')
    print(f'    early slope {early:.3f} BGCs per BGC sampled')
    print(f'    late  slope {late:.3f}')
    if ratio is not None:
        verdict = ('still growing linearly — assume it keeps growing'
                   if ratio > 0.8 else
                   'flattening — the dense cluster is well covered'
                   if ratio < 0.5 else 'partially saturating')
        print(f'    ratio       {ratio:.2f}   {verdict}')

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    xs = [r['depth'] for r in rows]
    ax.fill_between(xs, [r['largest_min'] for r in rows],
                    [r['largest_max'] for r in rows],
                    color='#a65628', alpha=0.18, lw=0, label='min-max over replicates')
    ax.plot(xs, [r['largest_mean'] for r in rows], color='#a65628', lw=1.8,
            label='largest component')
    ax.plot(xs, [r['components_mean'] for r in rows], color='#22645f', lw=1.5,
            ls='--', label='number of components')
    ax.plot(xs, xs, color='#98938a', lw=0.8, ls=':', label='y = x (all one component)')
    ax.set_xlabel('BGCs sampled')
    ax.set_ylabel('BGCs')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc='upper left')
    fig.tight_layout()
    for ext in ('png', 'svg'):
        p = str(args.out / f'component_rarefaction.{ext}')
        fig.savefig(p, dpi=200, metadata=SVG_METADATA if ext == 'svg' else None)
        if ext == 'svg':
            canonicalise_svg(p)
    plt.close(fig)
    print(f'  wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
