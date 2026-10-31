#!/usr/bin/env python3
"""Figure 2 — the pepM pre-screen: large runtime saving, no loss of BGCs.

Note what this figure does NOT claim. The screen does not *find* more BGCs --
it cannot, since it only decides which genomes antiSMASH sees. Its design goal
is to be detection-neutral while removing most of the compute, so the evidence
that matters is that the BGC count is UNCHANGED. Panel B is therefore a
negative result presented as one, and the Bacteroides pair is included
precisely because it is the case where neutrality initially failed: a
pseudogene-flagged pepM carries no /translation, so the screen dropped two
genomes until parse_genome was taught to translate such a CDS itself.

The screening cost annotations read "screening work", not "screen stage": the
screen was its own stage when these were measured and now runs inside
FETCH_RENAME_SCREEN, so the work is the same and the accounting boundary is not.
See the comment on COST.

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
# The fourth column was measured when the screen was its own stage, PEPM_PRESCREEN,
# and could be read straight off the trace. It no longer is: the screen now runs
# inside FETCH_RENAME_SCREEN, so that a genome the screen rejects is deleted before
# its output is declared and never reaches publishDir. The screen does the same work
# on the same genomes for the same verdicts -- verified identical on P. ananatis
# (193 of 344) and Winslowiella -- but it is no longer a separately timed task, so
# these figures are labelled as what the screening WORK costs rather than as a stage.
# Re-measuring would mean eight runs, the screen-off arms of which are the expensive
# ones (Erwiniaceae alone is 2,771 genomes through antiSMASH, ~61 CPU-h), to move
# numbers by a few percent and change no conclusion.
# TWO measured full taxa, compared on antiSMASH CPU -- the thing the screen acts on.
#
# This panel reports ANTISMASH CPU-minutes, not total pipeline CPU-minutes, and the
# change matters. The Erwiniaceae totals (4,508 -> 1,497) carry 741 and 655 CPU-min of
# GTDB-Tk, which the screen does not touch, and the screen-on arm additionally has
# RECOVER_ORFS (90 min) that the screen-off arm predates. Comparing totals mixes the
# screen with two unrelated changes; docs/comparisons/pepm_prescreen/summary.json says
# so directly -- "antiSMASH CPU-minutes are comparable; task counts are not".
#
# Both rows are full taxa run BOTH ways. Earlier versions of this panel carried five
# datasets ordered by prevalence; two of them (Actinomycetes 323, B. fragilis 136) were
# curated probes rather than taxa and supplied the endpoints the trend rested on, so
# the trend is gone and the claim is now just cost-at-equal-yield.
#
# Enterobacterales is NOT here. Its 46x needed a projected arm -- 1.35 M antiSMASH runs
# that were never performed -- and a projected bar beside measured ones invites exactly
# the reading it should not get. The projection belongs in the text.
#
# (clade, antiSMASH CPU-min screen off, screen on, the screen's own CPU-min)
COST = [
    ('Erwiniaceae\n2,771 genomes\n11% BGC-positive', 3680.1, 568.9, 91.3),
    ('P. ananatis\n344 genomes\n56% BGC-positive', 528.7, 311.3, 10.5),
]

# (clade, regions found with the screen OFF, with it ON).
# Panel B carries the same two taxa as panel A, and only those.
#
# B. fragilis was here as the held-out-clade sensitivity probe and is where the
# pseudogene blind spot was found. It is a curated set, not a taxon, so it is written
# up in the manuscript text instead -- a bar beside a taxon measurement cannot carry
# the caveat that makes it meaningful.
DETECT = [
    ('Erwiniaceae\n2,771 genomes', 333, 333),
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
    ax.set_title('Two full taxa, each run both ways', pad=24)
    bar_values(ax, list(b1) + list(b2), fmt='{:,.0f}', dy=0.0)
    top = max(off)
    for i, (_, o, n, screen) in zip(x, COST):
        ax.text(i, top * 3.4, f'{o / n:.1f}x less', ha='center', fontsize=9.5,
                color=AFTER, fontweight='bold')
        ax.text(i, top * 2.1, f'the screen itself cost {screen:.0f} CPU-min',
                ha='center', fontsize=7.2, color=FAINT)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlim(-0.8, len(labels) - 1 + 0.8)
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
    fig.suptitle('The pepM pre-screen cuts antiSMASH compute 1.7\u20136.5\u00d7 on full taxa \u2014 and finds every BGC',
                 fontsize=11.5, fontweight='bold', y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, args.outdir, 'fig2_prescreen')


if __name__ == '__main__':
    main()
