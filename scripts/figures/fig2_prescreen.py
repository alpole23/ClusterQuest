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
COST = [
    ('Actinomycetes\n323 genomes\n9.6% positive', 1618.8, 179.4, 12.7),
    ('Erwiniaceae\n2,771 genomes\n11% positive', 4508.0, 1497.0, 91.0),
    ('P. ananatis\n344 genomes\n56% positive', 603.0, 346.0, 10.5),
    ('B. fragilis\n136 genomes\n72% positive', 277.4, 247.3, 3.6),
]

# (clade, regions found with the screen OFF, with it ON). The Bacteroides pair is
# the pre-fix number; see PSEUDO_NOTE.
DETECT = [
    ('Actinomycetes\n323 genomes', 32, 32),
    ('P. ananatis\n344 genomes', 226, 226),
    ('B. fragilis\n136 genomes', 143, 141),
]
BACT_FIXED = 143


def panel_cost(ax):
    labels = [c[0] for c in COST]
    off = [c[1] for c in COST]
    on = [c[2] for c in COST]
    x = range(len(labels))
    w = 0.36
    b1 = ax.bar([i - w / 2 for i in x], off, w, label='antiSMASH only', color=BEFORE)
    b2 = ax.bar([i + w / 2 for i in x], on, w, label='pre-screen + antiSMASH', color=AFTER)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yscale('log')
    ax.set_ylabel('total pipeline CPU-minutes (log)')
    ax.set_title('The saving tracks how dilute the taxon is', pad=24)
    ax.set_ylim(40, max(off) * 9)
    bar_values(ax, list(b1) + list(b2), fmt='{:,.0f}', dy=0.0)
    for i, (_, o, n, screen) in enumerate(COST):
        ax.text(i, max(off) * 4.4, f'{o / n:.1f}× less',
                ha='center', fontsize=9, color=AFTER, fontweight='bold')
        ax.text(i, max(off) * 2.1, f'screening work: {screen:.0f}',
                ha='center', fontsize=7.2, color=FAINT)
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
    ax.set_ylabel('phosphonate BGC regions detected')
    ax.set_title('and costs no BGCs — once a blind spot was closed', pad=24)
    ax.set_ylim(0, max(off) * 1.46)
    bar_values(ax, list(b1) + list(b2), fmt='{:.0f}')

    for i, (_, o, n) in enumerate(DETECT):
        if o == n:
            ax.text(i, max(off) * 1.31, 'identical', ha='center', fontsize=8.5,
                    color=GOOD, fontweight='bold')
            ax.text(i, max(off) * 1.22, '100% sensitivity', ha='center',
                    fontsize=7.2, color=FAINT)
        else:
            # Do NOT draw the corrected bar. 143 against 141 is 1.4% of this axis,
            # an invisible mark pretending to be evidence. State it in words.
            ax.text(i, max(off) * 1.31, f'−{o - n} before fix', ha='center',
                    fontsize=8.5, color=ACCENT, fontweight='bold')
            ax.text(i, max(off) * 1.22, f'{BACT_FIXED}/{o} after fix', ha='center',
                    fontsize=7.6, color=GOOD, fontweight='bold')
            ax.annotate('pseudogene pepM carries\nno /translation',
                        xy=(i + w / 2, n), xytext=(i - 0.05, max(off) * 0.74),
                        ha='center', fontsize=7.2, color=ACCENT,
                        arrowprops=dict(arrowstyle='->', color=ACCENT, lw=0.9))
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

    fig.suptitle('The pepM pre-screen removes up to 9× the compute '
                 'without removing a single BGC',
                 fontsize=11.5, fontweight='bold', y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, args.outdir, 'fig2_prescreen')


if __name__ == '__main__':
    main()
