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
# ONE measured example, and it is a FULL TAXON.
#
# This panel used to carry five datasets ordered by prevalence, under the claim that
# the saving tracks how dilute the taxon is. That claim was confounded: two of the
# five -- Actinomycetes (323 genomes) and B. fragilis (136) -- are CURATED TEST SETS,
# assembled to probe screen sensitivity on a held-out clade, not full taxa. They
# supplied the two endpoints the trend rested on, 9.0x and 1.1x, so a monotonic
# reading of five points was really a reading of two hand-picked ones. Removing them
# removes the trend, and the trend was not the claim worth making.
#
# What is left is the claim the screen actually supports: on a full taxon, run both
# ways, it costs a fraction of the compute and finds the same BGCs. Erwiniaceae is
# the only full taxon at meaningful scale where BOTH arms were genuinely run.
COST = [
    ('Erwiniaceae\n2,771 genomes\n11% BGC-positive', 4508.0, 1497.0, 91.0),
]

# Enterobacterales is kept OUT of COST and drawn separately, because one of its two
# bars would be projected. The screened arm is measured: 1,349,507 GenBank genomes, the
# whole run 649 CPU-h. The without-screen arm is 1,349,507 antiSMASH runs at the 80.4
# CPU-s/genome this project measures, ~30,100 CPU-h, which was not spent and will not be.
# Putting a projected bar beside four measured pairs would read as a fifth measurement.
# It is annotated instead, on the same axes, marked as projected.
#
# Updated 2026-10-17 from the RefSeq run (150,690 genomes, 0.84% positive, 32x) to the
# GenBank run at 8.9x the input. The claim gets stronger but NOT in proportion to the
# retention difference: retention fell 0.84% -> 0.11%, yet the saving only goes 32x ->
# 46x, because the screened arm's cost is dominated by FETCH_RENAME_SCREEN (589 of
# 649 CPU-h). The screen saves antiSMASH time; it does not save download time.
#
# CPU-minutes. 649 CPU-h is actual consumption (realtime x %cpu), matching how the
# RefSeq baseline was derived -- the allocation-based figure reads 755 CPU-h and is not
# comparable, because cpus/task changed from 4 to 1 mid-run. See
# docs/comparisons/enterobacterales_genbank/PREDICTIONS.md.
ENTERO = ('Enterobacterales\n1,349,507 genomes\n0.11% positive', 1808119.0, 38940.0)

# (clade, regions found with the screen OFF, with it ON). The Bacteroides pair is
# the pre-fix number; see PSEUDO_NOTE.
# Panel B answers a different question from panel A, so it keeps a second dataset.
#
# Erwiniaceae is the full taxon, run both ways: 333 regions either way, 298/298
# sensitivity, and the screened arm recovered all 333 with identical nucleotide
# sequences (5,845,237 bases). B. fragilis STAYS even though it is a curated
# held-out clade, because that is exactly its job -- it supplies no prevalence
# claim, it is the probe that found the pseudogene blind spot, and dropping it
# would delete the figure's only admission that the screen ever missed anything.
# It is labelled as held-out so it cannot be read as a taxon measurement.
DETECT = [
    ('Erwiniaceae\n2,771 genomes\nfull taxon', 333, 333),
    ('B. fragilis\n136 genomes\nheld-out clade', 143, 141),
]
BACT_FIXED = 143


def panel_cost(ax):
    # Erwiniaceae is the measurement; Enterobacterales is drawn beside it as the
    # scale this is for, hatched and labelled, because its without-screen arm was
    # never run -- 1.35 M antiSMASH invocations at the 80.4 CPU-s/genome measured
    # here, ~30,100 CPU-h that was not spent and will not be. A projected bar beside
    # a measured pair must never read as a second measurement.
    lab_e, eoff, eon = ENTERO
    w = 0.36
    ax.bar(0 - w / 2, eoff, w, color='none', edgecolor=BEFORE, hatch='///', lw=1.2)
    ax.bar(0 + w / 2, eon, w, color=AFTER)
    ax.text(0, eoff * 1.9, f'{eoff / eon:.0f}x less', ha='center', fontsize=9.5,
            color=AFTER, fontweight='bold')
    ax.text(0, eoff * 0.55, 'without-screen arm\nPROJECTED, not run',
            ha='center', va='top', fontsize=7.4, color=BEFORE, style='italic')

    labels = [c[0] for c in COST]
    off = [c[1] for c in COST]
    on = [c[2] for c in COST]
    x = [i + 1 for i in range(len(labels))]
    b1 = ax.bar([i - w / 2 for i in x], off, w, label='antiSMASH only', color=BEFORE)
    b2 = ax.bar([i + w / 2 for i in x], on, w, label='pre-screen + antiSMASH', color=AFTER)
    ax.set_yscale('log')
    ax.set_ylabel('total pipeline CPU-minutes (log)')
    ax.set_title('One full taxon, run both ways', pad=24)
    bar_values(ax, list(b1) + list(b2), fmt='{:,.0f}', dy=0.0)
    for i, (_, o, n, screen) in zip(x, COST):
        ax.text(i, eoff * 1.9, f'{o / n:.1f}x less',
                ha='center', fontsize=9, color=AFTER, fontweight='bold')
        ax.text(i, eoff * 0.75, f'screening work: {screen:.0f}',
                ha='center', fontsize=7.2, color=FAINT)
    ax.set_xticks([0] + x)
    ax.set_xticklabels([lab_e] + labels, fontsize=8)
    ax.set_xlim(-0.6, len(labels) + 0.6)
    ax.set_ylim(40, eoff * 5.5)
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

    # The headline rests on the MEASURED pair. The 46x is projected and is
    # labelled as such inside panel A; putting it in the title would make the
    # figure's loudest claim the one number that was never run.
    fig.suptitle('The pepM pre-screen cuts antiSMASH compute threefold on a full taxon \u2014 and costs no BGCs',
                 fontsize=11.5, fontweight='bold', y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, args.outdir, 'fig2_prescreen')


if __name__ == '__main__':
    main()
