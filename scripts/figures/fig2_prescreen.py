#!/usr/bin/env python3
"""Figure 2 — the pepM pre-screen: large runtime saving, no loss of BGCs.

Note what this figure does NOT claim. The screen does not *find* more BGCs --
it cannot, since it only decides which genomes antiSMASH sees. Its design goal
is to be detection-neutral while removing most of the compute, so the evidence
that matters is that the BGC count is UNCHANGED. Panel B is therefore a
negative result presented as one.

Both taxa are controlled A/B pairs -- one taxon, one code version, one differing
parameter, identical antiSMASH task granularity. The last condition is the one
that was missing before 2026-11-01 and that cost the old Erwiniaceae row most of
its headline ratio; see the comment on COST.

Sources, all committed under docs/comparisons/pepm_prescreen/.

Usage:
    python scripts/figures/fig2_prescreen.py --outdir docs/figures
"""
import argparse
from pathlib import Path

from figure_style import (ACCENT, AFTER, BEFORE, FAINT, GOOD, INK,
                          bar_values, panel_label, plt, save)

ROOT = Path(__file__).resolve().parent.parent.parent

# (clade, genomes, % BGC-positive, total CPU-min screen off, screen on, screen cost)
# Ordered by BGC prevalence, lowest first: the saving is monotonic in it, and
# that ordering is the whole point of the panel.
#
# BOTH rows are now controlled A/B pairs: one taxon, one code version, one differing
# parameter, both arms aggregated on the same basis. That was not true before 2026-11-01.
#
# The old Erwiniaceae row read 3,680 -> 569 CPU-min (6.5x). It was internally consistent --
# both its arms are realtime x cpus -- but UNCONTROLLED: its screen-on arm batched antiSMASH
# and ran RECOVER_ORFS, its screen-off arm did neither, and GTDB-Tk differed. Re-run as a
# controlled pair the figure is 3,450.8 -> 629.1 = 5.49x on that same basis, or
# 1,984.4 -> 386.1 = 5.14x on realtime x %cpu, which is what this figure plots.
#
# A RETRACTED EXPLANATION, left here because it was published: this comment briefly blamed
# antiSMASH task granularity (2,760 one-genome tasks against 7 batched ones). That is wrong.
# modules/analysis/antismash.nf loops "for GENOME in ${genomes}" and invokes antiSMASH once
# per genome even inside a batched task, so batching amortises Nextflow task overhead and
# nothing of antiSMASH's own start-up. Like-for-like, antiSMASH per-genome cost barely moved
# between the two eras: 43.3 -> 41.6 CPU-s on x %cpu, 80.1 -> 72.3 on x cpus. The apparent
# "80.4 -> 41.6" collapse was entirely a unit mismatch between two summary.json files.
#
# THE RULE THIS LEAVES: a CPU figure is meaningless until you name its basis, and two ratios
# are comparable only on the same basis and the same code. Panel A is realtime x %cpu.
#
# AND A SECOND RETRACTION, 2026-11-01, of the line that used to sit here: it said "P. ananatis
# needed no correction -- 1.7x is a within-pair ratio on one basis." The ratio is indeed
# within-pair valid. But this panel plots ABSOLUTE CPU-minutes, and controlled_pantoea_ananatis
# measured them as realtime x cpus (its summary.json method field says so), against Erwiniaceae's
# realtime x %cpu. So the BAR HEIGHTS crossed bases even though each ratio did not -- the exact
# error this file's own rule forbids, surviving because the rule was checked against the ratio
# and not against the geometry. Its traces are gone (only tasks_*.tsv, realtime and cpus, no
# %cpu), so it cannot be re-expressed, and P. ananatis is therefore OUT of panel A. Its 1.7x
# stays in manuscript section 2.3 with its basis named. Panel B keeps it: region COUNTS have no
# CPU basis, so the two-taxon detection-neutrality claim is unaffected.
#
# Screening cost is reported as what the screening WORK costs, not as a stage: the
# screen runs inside FETCH_RENAME_SCREEN rather than as its own timed task, so on
# Erwiniaceae it is an arm-to-arm delta and on P. ananatis (measured 2026-09-17, when
# PEPM_PRESCREEN still existed) it is a stage total.
#
# This panel reports ANTISMASH CPU-minutes, not total pipeline CPU-minutes, and the
# change matters. Whole-pipeline CPU for the controlled Erwiniaceae pair is 2,287 -> 498
# (4.60x), which mixes in RECOVER_ORFS -- itself 18x cheaper with the screen on, because
# it too only sees screened genomes. GTDB-Tk is off in both arms. Comparing totals
# attributes to the screen work the screen did not do.
#
# Earlier versions of this panel carried five datasets ordered by prevalence; two of them
# (Actinomycetes 323, B. fragilis 136) were curated probes rather than taxa and supplied the
# endpoints the trend rested on, so the trend is gone and the claim is now cost-at-equal-yield
# on the one taxon whose cost is measured on this basis.
#
# Enterobacterales is NOT here. Its 46x needed a projected arm -- 1.35 M antiSMASH runs
# that were never performed -- and a projected bar beside measured ones invites exactly
# the reading it should not get. The projection belongs in the text. P. ananatis is absent
# from panel A for the same reason in a different variable: a bar on another accounting
# basis beside measured ones invites the same bad reading, however the caption hedges.
#
# (clade, antiSMASH CPU-min screen off, screen on, the screen's own CPU-min)
COST = [
    ('Erwiniaceae\n2,863 genomes\n11% BGC-positive', 1984.4, 386.1, 68.8),
]

# (clade, regions found with the screen OFF, with it ON).
# Panel B carries BOTH taxa, unlike panel A: these are region counts, which carry no CPU
# accounting basis, so P. ananatis is directly comparable here and not there.
#
# B. fragilis was here as the held-out-clade sensitivity probe and is where the
# pseudogene blind spot was found. It is a curated set, not a taxon, so it is written
# up in the manuscript text instead -- a bar beside a taxon measurement cannot carry
# the caveat that makes it meaningful.
DETECT = [
    ('Erwiniaceae\n2,863 genomes', 335, 335),
    ('P. ananatis\n344 genomes', 226, 226),
]

def panel_cost(ax):
    labels = [c[0] for c in COST]
    off = [c[1] for c in COST]
    on = [c[2] for c in COST]
    x = list(range(len(labels)))
    w = 0.36
    b1 = ax.bar([i - w / 2 for i in x], off, w, label='antiSMASH only', color=BEFORE)
    b2 = ax.bar([i + w / 2 for i in x], on, w, label='pre-screen + antiSMASH', color=AFTER)
    ax.set_yscale('log')
    ax.set_ylabel('antiSMASH CPU-minutes (log)')
    ax.set_title('A dilute taxon, run both ways', pad=24)
    bar_values(ax, list(b1) + list(b2), fmt='{:,.0f}', dy=0.0)
    top = max(off)
    for i, (_, o, n, screen) in zip(x, COST):
        ax.text(i, top * 3.4, f'{o / n:.1f}x less', ha='center', fontsize=9.5,
                color=AFTER, fontweight='bold')
        ax.text(i, top * 2.1, f'the screen itself cost {screen:.0f} CPU-min',
                ha='center', fontsize=7.2, color=FAINT)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    # One pair stretched across the full axis reads as a chart with a series missing;
    # hold it to the middle, as panel B already does.
    ax.set_xlim(-1.4, len(labels) - 1 + 1.4)
    ax.set_ylim(40, top * 7.5)
    ax.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, -0.19))

def panel_detection(ax):
    labels = [d[0] for d in DETECT]
    off = [d[1] for d in DETECT]
    on = [d[2] for d in DETECT]
    x = range(len(labels))
    w = 0.36
    b1 = ax.bar([i - w / 2 for i in x], off, w, label='antiSMASH only', color=BEFORE)
    b2 = ax.bar([i + w / 2 for i in x], on, w, label='pre-screen + antiSMASH', color=AFTER)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    # A single pair stretched across the full axis reads as a bar chart with a
    # missing series. Hold the bars to the middle instead.
    ax.set_xlim(-1.4, len(labels) - 1 + 1.4)
    ax.set_ylabel('phosphonate BGC regions detected')
    ax.set_title('and finds every BGC the unscreened run finds', pad=24)
    ax.set_ylim(0, max(off) * 1.46)
    bar_values(ax, list(b1) + list(b2), fmt='{:.0f}')

    for i, (_, o, n) in enumerate(DETECT):
        assert o == n, f'unexpected detection difference in {DETECT[i][0]}'
        ax.text(i, max(off) * 1.31, 'identical', ha='center', fontsize=8.5,
                color=GOOD, fontweight='bold')
        ax.text(i, max(off) * 1.22, '100% sensitivity', ha='center',
                fontsize=7.2, color=FAINT)
    ax.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, -0.19))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='docs/figures')
    args = ap.parse_args()

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8))
    panel_cost(axes[0])
    panel_detection(axes[1])
    for ax, letter in zip(axes, 'AB'):
        panel_label(ax, letter, dx=-0.13, dy=1.19)

    # The headline rests on the MEASURED pair. The 46x is projected and is
    # labelled as such inside panel A; putting it in the title would make the
    # figure's loudest claim the one number that was never run.
    # Reports the measured range on antiSMASH CPU. The Enterobacterales 46x is
    # deliberately absent: its without-screen arm was never run.
    fig.suptitle('The pepM pre-screen cuts antiSMASH compute 5.1\u00d7 on a dilute taxon \u2014 and finds every BGC',
                 fontsize=11.5, fontweight='bold', y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, args.outdir, 'fig2_prescreen')


if __name__ == '__main__':
    main()
