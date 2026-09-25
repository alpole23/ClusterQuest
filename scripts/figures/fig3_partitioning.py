#!/usr/bin/env python3
"""Figure 3 — BiG-SCAPE partitioning: what it costs, what it buys, what it leaves alone.

The honest claim, and the one the data supports: partitioning does NOT make
clustering faster at any scale we have measured -- it is 1.7-3.7x SLOWER,
because every partition re-pays BiG-SCAPE's fixed Pfam-load cost. What it buys
is memory, and memory is the wall: at the ~121,000 BGCs a million genomes
yields, one job needs ~1.9 TB, which is not a machine most groups have. It buys
that without touching the answer -- ARI 1.0000, identical family counts.

Sources, all measured and committed:
  docs/benchmark_data/bigscape_scaling.tsv       peak RSS and CPU vs BGC count
  docs/benchmark_data/bigscape_fit.json          the fits
  docs/benchmark_data/partition_check_*.json     ARI and family counts
  CLAUDE.md "BiG-SCAPE Partitioning"             the paired runtimes

Usage:
    python scripts/figures/fig3_partitioning.py --outdir docs/figures
"""
import argparse
import csv
import json
from pathlib import Path

from figure_style import (AFTER, BEFORE, ACCENT, FAINT, GOOD, INK,
                          bar_values, panel_label, plt, save)

ROOT = Path(__file__).resolve().parent.parent.parent
BENCH = ROOT / 'docs' / 'benchmark_data'

# Measured peak-RSS fit from the scaling benchmark (CLAUDE.md, BiG-SCAPE scaling).
MEM_INTERCEPT, MEM_COEFF = 1.14, 1.29e-7

# Paired runtimes: (label, BGCs, one-job seconds, partitioned seconds, partitions).
# One-job figures for 185 and 518 are from the partitioning benchmark; 333 is the
# scaling benchmark's own 333-BGC wall time on the same box.
RUNTIMES = [
    ('Streptomyces\n185 BGCs', 185, 30.0, 112.0, 19),
    ('Erwiniaceae\n333 BGCs', 333, 45.1, 77.7, 6),
    ('Combined\n518 BGCs', 518, 93.0, 191.0, 24),
]

# Family counts, unpartitioned vs partitioned at the safe 0.60 cut. All ARI 1.0000.
#
# All three rows were measured at the 5 kb neighbourhood, which is what antiSMASH's
# own phosphonate rule declares. The pipeline default is now 10 kb, and the counts
# there are different -- Erwiniaceae reads 22 families rather than 19. They are left
# as measured: re-running Streptomyces and Combined at 10 kb means re-running
# antiSMASH over both sets, and mixing windows between rows of one figure would be
# worse than a consistent older one.
#
# What matters is that the CLAIM is window-independent, and that was checked rather
# than assumed. Re-run at the 10 kb default on Erwiniaceae (2026-09-25, 334 BGCs,
# 6 partitions sized 236/89/4/2/2/1): 22 families unpartitioned, 22 partitioned,
# 18,568 co-membership pairs either way, 0 split, 0 merged, ARI 1.0000. Identical
# conclusion, different absolute numbers.
FAMILIES = [
    ('Streptomyces\n185 BGCs', 81, 81),
    ('Erwiniaceae\n333 BGCs', 19, 19),
    ('Combined\n518 BGCs', 100, 100),
]

TARGET_BGCS = 121_000          # BGCs a million genomes yields at observed prevalence
PARTITIONED_SHARE = 0.21       # largest component, Streptomyces-like diversity @0.60


def mem_gb(n):
    return MEM_INTERCEPT + MEM_COEFF * n ** 2


def load_measured_memory():
    rows = []
    with open(BENCH / 'bigscape_scaling.tsv') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r.get('exit') == '0':
                rows.append((int(r['bgcs']), float(r['max_rss_gb'])))
    return sorted(rows)


def check_fit(measured):
    """The fit is extrapolated 12x past its data, so verify it on the data it has."""
    worst = max(abs(mem_gb(n) - gb) / gb for n, gb in measured if n >= 4000)
    return worst


def panel_runtime(ax):
    labels = [r[0] for r in RUNTIMES]
    one = [r[2] for r in RUNTIMES]
    part = [r[3] for r in RUNTIMES]
    x = range(len(labels))
    w = 0.36
    b1 = ax.bar([i - w / 2 for i in x], one, w, label='one job', color=BEFORE)
    b2 = ax.bar([i + w / 2 for i in x], part, w, label='partitioned', color=ACCENT)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel('BiG-SCAPE wall time (s)')
    ax.set_title('The cost: partitioning is slower\nat every scale measured', pad=26)
    ax.set_ylim(0, max(part) * 1.42)
    bar_values(ax, list(b1) + list(b2), fmt='{:.0f}')
    for i, (_, _, o, p, nparts) in enumerate(RUNTIMES):
        ax.text(i, max(part) * 1.30, f'{p / o:.1f}× slower',
                ha='center', fontsize=8.5, color=ACCENT, fontweight='bold')
        ax.text(i, max(part) * 1.21, f'{nparts} partitions',
                ha='center', fontsize=7.3, color=FAINT)
    ax.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, -0.13))


def panel_memory(ax, measured):
    ns = [n for n, _ in measured]
    gbs = [gb for _, gb in measured]
    ax.plot(ns, gbs, 'o', color=INK, ms=4.5, label='measured peak RSS', zorder=5)

    # The quadratic was fitted on 1,500-4,000 BGCs and only describes the regime
    # where growth has started -- below ~4,000 peak RSS looks flat and the fit
    # does not apply. Solid over the measured range, dashed where it is being
    # extrapolated, so the figure cannot be read as 12x more measurement.
    fitted = [n for n in range(4000, 10001, 100)]
    extrap = [n for n in range(10000, TARGET_BGCS + 1, 250)]
    ax.plot(fitted, [mem_gb(n) for n in fitted], '-', color=BEFORE, lw=1.8,
            label='one job, fit over measured range')
    ax.plot(extrap, [mem_gb(n) for n in extrap], '--', color=BEFORE, lw=1.6,
            label='one job, extrapolated')

    part_n = int(TARGET_BGCS * PARTITIONED_SHARE)
    ax.plot(extrap, [mem_gb(int(n * PARTITIONED_SHARE)) for n in extrap],
            '--', color=GOOD, lw=1.8, label='largest partition, extrapolated')

    one_gb, part_gb = mem_gb(TARGET_BGCS), mem_gb(part_n)
    ax.plot([TARGET_BGCS], [one_gb], 'o', color=ACCENT, ms=7, zorder=6)
    ax.plot([TARGET_BGCS], [part_gb], 'o', color=GOOD, ms=7, zorder=6)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(280, 420_000)
    ax.set_ylim(0.8, 9000)

    for gb, txt in ((128, '128 GB node'), (1024, '1 TB node')):
        ax.axhline(gb, ls=':', lw=0.9, color=FAINT)
        ax.text(300, gb * 1.15, txt, fontsize=7.2, color=FAINT)

    ax.annotate(f'{one_gb:,.0f} GB', xy=(TARGET_BGCS, one_gb), xytext=(10, 2),
                textcoords='offset points', ha='left', fontsize=9,
                color=ACCENT, fontweight='bold')
    ax.annotate(f'{part_gb:,.0f} GB', xy=(TARGET_BGCS, part_gb), xytext=(10, -3),
                textcoords='offset points', ha='left', fontsize=9,
                color=GOOD, fontweight='bold')

    ax.set_xlabel('BGCs clustered')
    ax.set_ylabel('peak memory (GB)')
    ax.set_title(f'The gain: {one_gb - part_gb:,.0f} GB at the scale\n'
                 'a million genomes reaches', pad=26)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20), fontsize=7.2, ncol=2)
    return one_gb, part_gb


def panel_families(ax):
    labels = [f[0] for f in FAMILIES]
    one = [f[1] for f in FAMILIES]
    part = [f[2] for f in FAMILIES]
    x = range(len(labels))
    w = 0.36
    b1 = ax.bar([i - w / 2 for i in x], one, w, label='one job', color=BEFORE)
    b2 = ax.bar([i + w / 2 for i in x], part, w, label='partitioned', color=AFTER)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel('gene cluster families')
    ax.set_title('The answer is unchanged:\nidentical families, ARI 1.0000', pad=26)
    ax.set_ylim(0, max(one) * 1.42)
    bar_values(ax, list(b1) + list(b2), fmt='{:.0f}')
    for i in x:
        ax.text(i, max(one) * 1.30, 'ARI 1.0000', ha='center', fontsize=8.5,
                color=GOOD, fontweight='bold')
        ax.text(i, max(one) * 1.21, '0 split, 0 merged', ha='center',
                fontsize=7.3, color=FAINT)
    ax.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, -0.13))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='docs/figures')
    args = ap.parse_args()

    measured = load_measured_memory()
    resid = check_fit(measured)
    print(f'memory fit worst residual on measured points >=4000 BGCs: {resid:.1%}')

    # Cross-check the committed ARI claims rather than trusting the constants above.
    for name in ('partition_check_0.6', 'partition_check_combined_0.6'):
        d = json.load(open(BENCH / f'{name}.json'))
        assert d['adjusted_rand_index'] == 1.0, name
        assert d['families_reference'] == d['families_partitioned'], name
        print(f'{name}: {d["bgcs_compared"]} BGCs, '
              f'{d["families_reference"]} families both ways, ARI '
              f'{d["adjusted_rand_index"]:.4f}')

    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.5))
    panel_runtime(axes[0])
    one_gb, part_gb = panel_memory(axes[1], measured)
    panel_families(axes[2])
    for ax, letter in zip(axes, 'ABC'):
        panel_label(ax, letter, dx=-0.13, dy=1.22)

    fig.suptitle(
        'Partitioning BiG-SCAPE: a 1.7–3.7× runtime cost buys '
        f'{one_gb - part_gb:,.0f} GB and changes no cluster assignment',
        fontsize=11.5, fontweight='bold', y=1.02)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, args.outdir, 'fig3_partitioning')


if __name__ == '__main__':
    main()
