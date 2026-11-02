#!/usr/bin/env python3
"""Figure 4 — the two mechanisms that make a large taxon runnable on one workstation.

This pipeline was developed and run on a single 16-core / 56 GB box, with no
cluster. Two features exist for that reason and this figure prices both of them:
the pepM pre-screen, which decides how much data reaches antiSMASH, and
BiG-SCAPE partitioning, which decides whether clustering fits in RAM at all.

WHAT EACH PANEL IS ALLOWED TO CLAIM, and why the three differ in strength:

  A  pre-screen CPU      MEASURED on both arms for Erwiniaceae. The
                         Enterobacterales without-screen arm was never run and is
                         drawn hatched, as a projection from the measured
                         per-genome rate.
  B  partitioning memory PROJECTED from two measured ingredients -- the monolithic
                         memory curve and the largest-partition share. No
                         partitioned run has ever had its peak memory measured.
  C  partitioning cost   MEASURED, but it is CUMULATIVE COMPUTE, not wall time.

Panels B and C are split rather than drawn on one pair of axes. They are two
consequences of one intervention measured in different units, and a dual-axis
chart would invite a reader to compare a gigabyte against a second.

─── THE BASIS, WHICH IS THE WHOLE TRAP ───────────────────────────────────────

Panel A is `realtime x cpus` (cores RESERVED) throughout, including the
Erwiniaceae pair, which is reported elsewhere on `realtime x %cpu` (cores USED).
The two differ by about 2x here and may never be mixed.

The reason for choosing the less flattering basis is the Enterobacterales bar.
That run has no committed trace, and its per-process CPU figures are recoverable
as `realtime x cpus` by three independent tells -- BIGSCAPE reports 1.0 CPU-h
over a 7 m 21 s realtime, which is 8.16 effective cores against its 8 allocated,
a utilisation no real process achieves. So Enterobacterales exists ONLY on
`realtime x cpus`, and putting Erwiniaceae beside it on `x %cpu` would have been
exactly the cross-basis error this project has already published twice. On
`x %cpu` the Erwiniaceae antiSMASH ratio is 5.14x rather than 5.49x.

─── THE CLAIM PANEL A ACTUALLY MAKES ─────────────────────────────────────────

Not "the screen saves 5x". The saving is set by RETENTION, which is a property of
the taxon and not of the pipeline:

    Erwiniaceae      10.76% retention  ->  5.5x less antiSMASH CPU   (both arms measured)
    Enterobacterales  0.113% retention  ->  ~443x                    (one arm projected)

A 95x difference in retention gives an 81x difference in saving. Quoting a single
number for "what the screen saves" would be quoting a property of Erwiniaceae.

─── WHY PANEL B IS A PROJECTION AND MUST SAY SO ──────────────────────────────

BiG-SCAPE's peak memory is flat below ~4,000 BGCs and quadratic above it. Every
run that has ever been partitioned is at or below 1,302 BGCs -- deep inside the
flat region -- so the measured memory saving from partitioning, on every run on
record, is between 0.00 and 0.43 GB. The benefit is real and it is arithmetic
from the measured curve, but it has never been demonstrated, and the panel says
so rather than drawing a measured-looking point.

THE SHARE IS THE WHOLE ANSWER, AND IT IS NOT A CONSTANT. Partitioning cuts the
largest job from N BGCs to share x N, and memory is quadratic, so the saving goes
as 1/share^2. The five measured shares span 20.5% to 72.4%, which is 5.6x against
1.8x -- so a single "partitioning saves Nx" is not a thing this data supports. The
panel draws the band and labels both ends.

A RETRACTION THIS FIGURE EXISTS PARTLY TO RECORD. fig3_partitioning.py uses
PARTITIONED_SHARE = 0.21 and reports a 22.4x memory saving (2,792 -> 125 GB at
117,000 BGCs). 0.21 is the Streptomyces share -- the most favourable of the five,
measured on a 185-BGC test set that has no dominant pepM component -- while the
117,000 BGCs it is applied to are projected from the ERWINIACEAE rate, whose own
measured share is 0.707. Share and population come from different taxa. Carried
consistently at 0.707 the same arithmetic gives 1,396 GB and 2.0x, not 125 GB and
22.4x. The real-taxon figure is the one this panel plots.

─── WHY PANEL C SAYS "COMPUTE" AND NOT "WALL" ────────────────────────────────

`bigscape_seconds_partitioned` in the partition_check files is not a wall time.
scripts/bench_bigscape_partitioned.py accumulates it in a serial loop
(`total_s += secs`), one partition after another, so it is the SUM of N
subprocess elapsed times. The one-job arm is a single process's elapsed time. The
ratio is therefore cumulative compute against wall, and no partitioned WALL time
is committed anywhere in this repository -- partitions run concurrently under
Nextflow, so the real elapsed penalty is smaller than this and unmeasured.

Quoted as cumulative compute the number is sound and the mechanism is clear: each
partition re-pays BiG-SCAPE's fixed Pfam load, which is why 19 partitions of 185
BGCs cost more than one job of 185.

Sources, all committed:
  docs/benchmark_data/bigscape_scaling.tsv            monolithic memory + wall vs BGCs
  docs/benchmark_data/partition_check_*.json          partition counts, largest, seconds
  docs/benchmark_data/pepm_strep.json                 Streptomyces component decomposition
  docs/comparisons/bigscape_partition_equivalence/    Enterobacterales 1,302-BGC split
  docs/comparisons/bigscape_partitioning_10kb/        the 10 kb Erwiniaceae re-run
  docs/comparisons/pepm_prescreen/controlled_erwiniaceae/  the controlled A/B pair
  docs/comparisons/enterobacterales_genbank/summary.json   the 1.35 M genome run

Usage:
    python scripts/figures/fig4_scalability.py --outdir docs/figures
"""
import argparse
import csv
import json
from pathlib import Path

from figure_style import (ACCENT, AFTER, BEFORE, FAINT, GOOD, INK,
                          bar_values, panel_label, plt, save)

ROOT = Path(__file__).resolve().parent.parent.parent
BENCH = ROOT / 'docs' / 'benchmark_data'
COMP = ROOT / 'docs' / 'comparisons'

# ru_maxrss reports the max of any ONE process and BiG-SCAPE runs eleven. Total PSS
# summed across the tree -- which divides each shared page by its sharers, so the
# shared Pfam data is counted once -- is 1.58x higher, stable at 1.51-1.65 over five
# sizes (bigscape_partition_equivalence/memory_scaling_pss.tsv). Without this the
# curve understates what the machine must actually hold by a third.
MEM_RSS_TO_PSS = 1.58

# Where the quadratic starts. Below this peak memory is flat (1.62 -> 3.40 GB across
# the first seven sizes) and partitioning buys nothing at all -- which is the regime
# every partitioned run on record sits in.
QUADRATIC_ONSET = 4000

# This box's BiG-SCAPE allocation. BiG-SCAPE carries `process_high`, which the local
# profile pins to 16 GB (nextflow.config:468-469); the slurm profile overrides
# withName: BIGSCAPE to 64 GB (nextflow.config:519-520).
RAM_LINES = [(16, '16 GB — this box'), (64, '64 GB — slurm profile')]


def measured_memory():
    """The monolithic peak-memory curve, scaled to tree PSS. Measured, 333-10,000."""
    rows = []
    with open(BENCH / 'bigscape_scaling.tsv') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['exit'] == '0':
                rows.append((int(r['bgcs']),
                             float(r['max_rss_gb']) * MEM_RSS_TO_PSS,
                             float(r['wall_s'])))
    return sorted(rows)


def interp(pts, n):
    """Linear interpolation BETWEEN measured points. Returns None outside them.

    Deliberately not the quadratic fit. The fit was estimated on 1,500-4,000 BGCs
    and fig3_partitioning extrapolates it 12x past its data; this panel refuses to
    leave the measured range at all, so every gigabyte it prints sits between two
    points that were actually observed.
    """
    xs = [p[0] for p in pts]
    if n < xs[0] or n > xs[-1]:
        return None
    for (a, ga), (b, gb) in zip([(p[0], p[1]) for p in pts],
                                [(p[0], p[1]) for p in pts][1:]):
        if a <= n <= b:
            return ga + (gb - ga) * (n - a) / (b - a)
    return None


def partition_shares():
    """Every measured largest-partition share, read from its own record.

    The share is what sets the memory saving, and it is the quantity that varies
    most between taxa -- so it is read from the files rather than written here.
    """
    out = []

    def add(label, total, largest, src):
        out.append({'label': label, 'total': total, 'largest': largest,
                    'share': largest / total, 'source': src})

    # `partitioning` is a LIST of threshold sweeps, not a dict keyed by threshold.
    # 0.60 is the shipped default (bigscape_partition_identity); the same file also
    # records 0.5/0.7/0.8/0.9, and the share moves a long way across them --
    # 68.1% at 0.5 against 20.5% at 0.6 -- so the threshold has to be pinned here.
    s = json.loads((BENCH / 'pepm_strep.json').read_text())
    p = next(v for v in s['partitioning'] if abs(v['threshold'] - 0.6) < 1e-9)
    add('Streptomyces\n185 BGCs', s['sequences'], p['largest'],
        'pepm_strep.json @ threshold 0.6')

    d = json.loads((BENCH / 'partition_check_combined_0.6.json').read_text())
    add('Combined\n518 BGCs', d['bgcs_compared'], d['largest_partition'],
        'partition_check_combined_0.6.json')

    d = json.loads((BENCH / 'partition_check_0.6.json').read_text())
    add('Erwiniaceae\n333 BGCs', d['bgcs_compared'], d['largest_partition'],
        'partition_check_0.6.json')

    d = json.loads((COMP / 'bigscape_partitioning_10kb'
                    / 'partition_check_0.6.json').read_text())
    add('Erwiniaceae 10 kb\n334 BGCs', d['bgcs_compared'], d['largest_partition'],
        'bigscape_partitioning_10kb/partition_check_0.6.json')

    # The largest real taxon ever partitioned, and the share the figure trusts.
    add('Enterobacterales\n1,302 BGCs', 1302, 942,
        'bigscape_partition_equivalence/summary.json')
    return out


# (label, BGCs, one-job wall s, cumulative partitioned compute s, partitions)
#
# NOT a wall-time pair -- see the module docstring. The partitioned arm is a serial
# sum over partitions from bench_bigscape_partitioned.py; the one-job arm is a single
# process's elapsed time. The 10 kb Erwiniaceae re-run is excluded: its one-job arm is
# a Nextflow task realtime (1m 30s) that carries task overhead the others do not, so
# its ratio is not comparable with these three.
PARTITION_COST = [
    ('Streptomyces\n185 BGCs', 185, 30.0, 112.0, 19),
    ('Erwiniaceae\n333 BGCs', 333, 45.1, 77.7, 6),
    ('Combined\n518 BGCs', 518, 93.0, 191.0, 24),
]


def prescreen_cpu():
    """antiSMASH CPU-hours with and without the screen, on `realtime x cpus`.

    Erwiniaceae is the controlled pair: one taxon, one commit, one differing
    parameter, both arms measured. Enterobacterales has no without-screen arm --
    1,349,504 genomes were never put through antiSMASH -- so that bar is the
    controlled pair's measured per-genome rate applied to the GenBank genome count,
    and it is drawn hatched.
    """
    erw = json.loads((COMP / 'pepm_prescreen' / 'controlled_erwiniaceae'
                      / 'summary.json').read_text())
    gb = json.loads((COMP / 'enterobacterales_genbank' / 'summary.json').read_text())

    # The committed cpu_minutes block is `realtime x %cpu`. The x cpus figures live in
    # the file's prose (cpu_accounting_basis: "On realtime x cpus this run is ANTISMASH
    # 3,450.8 -> 629.1 CPU-min"), re-derived from both traces during the 2026-11-02
    # audit. Held here rather than parsed out of a sentence.
    erw_off_h, erw_on_h = 3450.8333 / 60, 629.1333 / 60
    erw_genomes = erw['genomes']
    erw_screened = erw['genomes_reaching_antismash']['screen_on']

    # Per-genome antiSMASH cost WITHOUT the screen, on x cpus: 3,450.83 CPU-min over
    # all 2,863 genomes. The with-screen arm's per-genome cost is 1.69x higher, because
    # the genomes it keeps are the BGC-bearing ones -- so using the off-arm rate for
    # the projection is the conservative direction.
    rate_s = erw_off_h * 3600 / erw_genomes

    ent_genomes = gb['run']['genomes']
    ent_passed = int(gb['found']['passed_screen'].split()[0].replace(',', ''))
    ent_on_h = float(gb['the_cpu_shape_has_inverted']['measured']['ANTISMASH']
                     .split()[0].replace(',', ''))
    ent_off_h = rate_s * ent_genomes / 3600

    return [
        {'label': 'Erwiniaceae\n2,863 genomes', 'off': erw_off_h, 'on': erw_on_h,
         'retention': erw_screened / erw_genomes, 'projected': False},
        {'label': f'Enterobacterales\n{ent_genomes:,} genomes',
         'off': ent_off_h, 'on': ent_on_h,
         'retention': ent_passed / ent_genomes, 'projected': True},
    ]


def panel_prescreen(ax, data):
    x = list(range(len(data)))
    w = 0.36
    off = [d['off'] for d in data]
    on = [d['on'] for d in data]
    b1 = ax.bar([i - w / 2 for i in x], off, w, color=BEFORE,
                label='antiSMASH on every genome')
    b2 = ax.bar([i + w / 2 for i in x], on, w, color=AFTER,
                label='pepM pre-screen first')
    # The projected bar is hatched and unfilled rather than coloured, so it cannot be
    # read as a measurement however the caption is skimmed. It also gets its own
    # legend entry: a hatch pattern is not self-explaining, and this is the one mark
    # in the figure that was never observed.
    for bar, d in zip(b1, data):
        if d['projected']:
            bar.set_facecolor('none')
            bar.set_edgecolor(ACCENT)
            bar.set_hatch('////')
            bar.set_linewidth(1.1)

    ax.set_yscale('log')
    ax.set_ylabel('antiSMASH CPU-hours (log)\nrealtime × cpus')
    ax.set_title('The pre-screen: the saving is set by\nretention, not by the pipeline',
                 pad=30)
    bar_values(ax, list(b1) + list(b2), fmt='{:,.0f}', dy=0.0)

    ax.set_ylim(1, max(off) * 30)
    # Annotations are placed in AXES fraction on y and DATA on x. Multiplying a data
    # value by a constant to "put it near the top" is a linear-axis habit: on a log
    # axis it lands wherever the decade happens to fall, which is what collided the
    # ratio label into the retention label in the first draft.
    yt = ax.get_xaxis_transform()
    for i, d in enumerate(data):
        # One decimal below 10x, none above: 5.48x printed as "5x" understated the
        # measured arm by half a multiple, and 443.0x adds nothing over 443x.
        r = d['off'] / d['on']
        ax.text(i, 0.95, f"{r:,.1f}× less" if r < 10 else f"{r:,.0f}× less",
                transform=yt, ha='center', va='top', fontsize=10, color=AFTER,
                fontweight='bold')
        ax.text(i, 0.875, f"{100 * d['retention']:.3g}% retained", transform=yt,
                ha='center', va='top', fontsize=7.2, color=FAINT)
    ax.set_xticks(x)
    ax.set_xticklabels([d['label'] for d in data], fontsize=8)
    ax.set_xlim(-0.8, len(data) - 1 + 0.8)

    # Hand the CONTAINERS to legend(), not their first Rectangle: `ax.bar(label=...)`
    # sets the label on the BarContainer, and a bare child renders as `_nolegend_`.
    from matplotlib.patches import Patch
    handles = [b1, Patch(facecolor='none', edgecolor=ACCENT, hatch='////',
                         linewidth=1.1, label='…projected; this arm was never run'),
               b2]
    ax.legend(handles=handles, loc='upper center', ncol=1,
              bbox_to_anchor=(0.5, -0.17), fontsize=8)


def panel_memory(ax, pts, shares):
    ns = [p[0] for p in pts]
    gbs = [p[1] for p in pts]
    ax.plot(ns, gbs, 'o-', color=BEFORE, ms=4.5, lw=1.8, zorder=5,
            label='one job (measured)')

    lo = min(s['share'] for s in shares)
    hi = max(s['share'] for s in shares)
    real = next(s for s in shares if s['total'] == 1302)   # the real-taxon share

    # Both edges of the band and the real-taxon line stay inside the measured range:
    # a point is drawn only where BOTH the total and share x total are interpolable.
    grid = [n for n in range(ns[0], ns[-1] + 1, 50)
            if interp(pts, n) is not None and interp(pts, int(n * lo)) is not None]
    band_lo = [interp(pts, int(n * lo)) for n in grid]
    band_hi = [interp(pts, int(n * hi)) for n in grid]
    ax.fill_between(grid, band_lo, band_hi, color=GOOD, alpha=0.18, lw=0,
                    label='largest partition, across the five measured pepM shares '
                          '(20.5–72.4%)')
    ax.plot(grid, [interp(pts, int(n * real['share'])) for n in grid], '-',
            color=GOOD, lw=2.0, zorder=5,
            label=f"largest partition at {100 * real['share']:.1f}% "
                  f"— the share real taxa measure")

    # Read the saving off the largest size benchmarked. The arrow is pulled in from
    # the last point so its label has room that is not the axis edge.
    N = ns[-1]
    one, part = interp(pts, N) or gbs[-1], interp(pts, int(N * real['share']))
    ax.annotate('', xy=(N, one), xytext=(N, part),
                arrowprops=dict(arrowstyle='<->', color=ACCENT, lw=1.4))
    ax.text(N * 1.1, (one * part) ** 0.5, f'{one - part:.0f} GB\n{one / part:.1f}×',
            ha='left', va='center', fontsize=9, color=ACCENT,
            fontweight='bold', linespacing=1.4)

    # The flat region, which is where every partitioned run actually sits.
    ax.axvspan(ns[0] * 0.85, QUADRATIC_ONSET, color=FAINT, alpha=0.08, lw=0, zorder=0)
    # Two lines, placed ABOVE the 64 GB reference line rather than across it. Three
    # lines starting lower ran straight through both RAM lines and their labels.
    ax.text(ns[0] * 0.95, 102,
            'flat below ~4,000 BGCs — every run ever\n'
            'partitioned sits in here, saving ≤0.43 GB',
            fontsize=7.0, color=FAINT, va='top', ha='left', linespacing=1.5)

    # RAM lines labelled at the far left and BELOW the line. Above it, the 64 GB label
    # met the underside of the flat-region note; right of the shaded region, it met
    # the note's second line, and the 16 GB label met the saving arrow.
    for gb, txt in RAM_LINES:
        ax.axhline(gb, ls=':', lw=0.9, color=FAINT, zorder=1)
        ax.text(ns[0] * 0.92, gb * 0.94, txt, fontsize=6.9, color=FAINT, va='top')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(ns[0] * 0.80, ns[-1] * 2.3)
    ax.set_ylim(1.9, 110)
    ax.set_xlabel('BGCs clustered')
    ax.set_ylabel('BiG-SCAPE peak memory, GB (log)\ntree PSS = ru_maxrss × 1.58')
    ax.set_title('What partitioning buys: memory, and only\nwhere the curve has turned',
                 pad=30)
    # Below the axes, as panels A and C do. Inside the plot there is no quadrant big
    # enough: the curve climbs through the middle, the flat-region note owns the upper
    # left, and the saving annotation owns the right edge.
    ax.legend(loc='upper center', ncol=1, bbox_to_anchor=(0.5, -0.17), fontsize=8)


def panel_cost(ax):
    labels = [c[0] for c in PARTITION_COST]
    one = [c[2] for c in PARTITION_COST]
    part = [c[3] for c in PARTITION_COST]
    x = list(range(len(labels)))
    w = 0.36
    b1 = ax.bar([i - w / 2 for i in x], one, w, color=BEFORE, label='one job')
    b2 = ax.bar([i + w / 2 for i in x], part, w, color=ACCENT,
                label='summed over partitions')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('BiG-SCAPE compute, seconds\n(partitioned arm is a serial sum)')
    ax.set_title('What it costs: every partition re-pays\nthe fixed Pfam load', pad=30)
    ax.set_ylim(0, max(part) * 1.52)
    bar_values(ax, list(b1) + list(b2), fmt='{:.0f}')
    # "3.7x more compute" over each pair ran into its neighbours at this panel width.
    # The unit belongs to the axis, which already says "BiG-SCAPE compute, seconds",
    # so the label over a bar only has to carry the multiplier.
    yt = ax.get_xaxis_transform()
    for i, (_, _, o, p, nparts) in enumerate(PARTITION_COST):
        ax.text(i, 0.96, f'{p / o:.1f}×', transform=yt, ha='center', va='top',
                fontsize=10, color=ACCENT, fontweight='bold')
        ax.text(i, 0.885, f'{nparts} partitions', transform=yt, ha='center',
                va='top', fontsize=7.3, color=FAINT)
    ax.legend(loc='upper center', ncol=1, bbox_to_anchor=(0.5, -0.17), fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='docs/figures')
    args = ap.parse_args()

    pts = measured_memory()
    shares = partition_shares()
    cpu = prescreen_cpu()

    # Print what the figure is asserting, so a regeneration that changes a number
    # says so on stdout rather than only in the rendered pixels.
    print('pre-screen, antiSMASH CPU-h on realtime x cpus:')
    for d in cpu:
        tag = 'PROJECTED off-arm' if d['projected'] else 'both arms measured'
        print(f"  {d['label'].replace(chr(10), ' '):42s} "
              f"{d['off']:>10,.1f} -> {d['on']:>7,.1f}  "
              f"{d['off'] / d['on']:>7,.1f}x  "
              f"retention {100 * d['retention']:.3g}%  [{tag}]")

    print('\nmeasured largest-partition shares:')
    for s in shares:
        print(f"  {s['label'].replace(chr(10), ' '):34s} "
              f"{s['largest']:>5,}/{s['total']:<6,} = {100 * s['share']:5.1f}%"
              f"   ({s['source']})")

    N = pts[-1][0]
    real = next(s for s in shares if s['total'] == 1302)
    one = pts[-1][1]
    print(f'\nmemory at {N:,} BGCs (largest size benchmarked): one job {one:.1f} GB')
    for s in sorted(shares, key=lambda s: s['share']):
        m = interp(pts, int(N * s['share']))
        print(f"  share {100 * s['share']:5.1f}% -> largest partition "
              f"{int(N * s['share']):>6,} BGCs = {m:5.1f} GB, "
              f"saves {one - m:5.1f} GB ({one / m:.1f}x)")
    print(f"  figure plots the real-taxon share {100 * real['share']:.1f}%")

    # Guard the claim the panel's shaded band makes. If a future run measures a share
    # outside this range the band moves, and the caption's numbers must move with it.
    assert 0.20 <= min(s['share'] for s in shares) <= 0.21, 'share range moved'
    assert 0.72 <= max(s['share'] for s in shares) <= 0.73, 'share range moved'

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 5.2))
    panel_prescreen(axes[0], cpu)
    panel_memory(axes[1], pts, shares)
    panel_cost(axes[2])
    for ax, letter in zip(axes, 'ABC'):
        panel_label(ax, letter, dx=-0.17, dy=1.25)

    fig.suptitle('Running a 1.35-million-genome taxon on one 16-core / 56 GB '
                 'workstation: what each mechanism buys, and what it costs',
                 fontsize=11.5, fontweight='bold', y=1.015)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.text(0.5, -0.055,
             'A: the Enterobacterales without-screen arm was never run — it is the '
             'controlled pair’s measured per-genome rate applied to 1,349,504 genomes. '
             'B: projected from the measured\nmonolithic curve and the measured '
             'partition shares; no partitioned run has had its memory measured. '
             'C: cumulative compute summed over partitions, not elapsed time.',
             ha='center', va='top', fontsize=7.4, color=FAINT, linespacing=1.6)
    save(fig, args.outdir, 'fig4_scalability')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
