#!/usr/bin/env python3
"""Figures 5 and S1 — what the pipeline does, and the three mechanisms that make it work.

Two figures from one generator, because they share the same facts and must not drift:

  fig5_mechanisms   main body. Three panels, one per mechanism, no stage overview.
  figS1_architecture  supplemental. Stages, data artefacts and the funnel in one view.

Every number is measured on the Enterobacterales RefSeq run (2026-09-27, 150,690
genomes, 43.1 h, 0 failures) and sourced from
docs/comparisons/enterobacterales_refseq/summary.json. Nothing here is illustrative.

The architecture has one defining fact and one surprise, and both figures carry them:

  the funnel    150,690 genomes -> 1,309 pass the pepM screen (0.87%) -> 1,303 regions
                -> 72 families. 99.1% is discarded before anything expensive runs.
  the inversion cost runs opposite to intuition. Fetch-and-screen is 40.4% of CPU over
                6,028 tasks; BiG-SCAPE, the one that is quadratic in both time and
                memory, is 0.7%. Moving and filtering data costs more than analysing it
                at this prevalence, which is why the screen is the architecture rather
                than an optimisation.

Usage:
    python scripts/figures/fig5_architecture.py --outdir docs/figures
"""
import argparse
import json
from pathlib import Path

from figure_style import (ACCENT, AFTER, BEFORE, FAINT, GOOD, INK,
                          panel_label, plt, save)
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parent.parent.parent
MEASURED = ROOT / 'docs' / 'comparisons' / 'enterobacterales_refseq' / 'summary.json'

# Role colours. Deliberately few: a reader should not have to learn a key to follow
# the flow, so hue marks the KIND of thing (compute / artefact / discard) and nothing
# is encoded by colour alone -- every box is also labelled.
C_COMPUTE = '#dce6ee'
C_ARTEFACT = '#f2ead9'
C_DISCARD = '#ebedef'
C_GATE = '#f7ddcb'


def measured():
    """The run's own numbers, so the figure cannot drift from the record."""
    d = json.loads(MEASURED.read_text())
    f, r = d['found'], d['run']
    passed = int(str(f['passed_screen']).split()[0].replace(',', ''))
    return {
        'genomes': r['genomes'], 'screened': f['screened'], 'passed': passed,
        'regions': f['regions'], 'positive': f['bgc_positive_genomes'],
        'gcfs': f['gcfs'], 'wall': r['wall'], 'tasks': r['tasks'],
        'cpu': d['the_cpu_shape_has_inverted']['measured'],
    }


def fit_text(ax, x, y, w, h, text, fs=8.2, weight='normal', colour=None,
             pad=0.90, minimum=5.4):
    """Place text centred in a box, wrapped and shrunk until it measures inside it.

    Hand-tuning font sizes per label does not survive the first edit to a string, and
    the first draft of this figure had five labels overflowing into their neighbours.
    This wraps on the box width and then shrinks until matplotlib's own measurement of
    the rendered extent fits, so a label can only overflow if it cannot fit at all.
    """
    import textwrap
    fig = ax.get_figure()
    fig.canvas.draw()                       # a renderer must exist to measure
    bb = ax.get_window_extent()
    box_px_w, box_px_h = w * bb.width * pad, h * bb.height * pad
    size = fs
    while size >= minimum:
        # Wrap on an estimated character width, then verify by measuring.
        char_px = size * fig.dpi / 72 * 0.58
        ncols = max(6, int(box_px_w / char_px))
        wrapped = (text if '\n' in text else
                   '\n'.join(textwrap.wrap(text, ncols, break_long_words=False,
                                           break_on_hyphens=False)))
        tmp = ax.text(x + w / 2, y + h / 2, wrapped, ha='center', va='center',
                      fontsize=size, color=colour or INK, linespacing=1.3,
                      fontweight=weight, zorder=4)
        ext = tmp.get_window_extent(fig.canvas.get_renderer())
        if ext.width <= box_px_w and ext.height <= box_px_h:
            return tmp
        tmp.remove()
        size -= 0.3
    return ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
                   fontsize=minimum, color=colour or INK, linespacing=1.3,
                   fontweight=weight, zorder=4)


def box(ax, x, y, w, h, text, fill=C_COMPUTE, edge=None, fs=8.2, weight='normal',
        style='round,pad=0.012,rounding_size=0.018'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, linewidth=0.9,
                                facecolor=fill, edgecolor=edge or '#aab4bf',
                                mutation_aspect=1, zorder=3))
    fit_text(ax, x, y, w, h, text, fs=fs, weight=weight)


# box() draws with boxstyle pad=0.012, so the visible edge sits that far outside the
# rect the caller gave. An arrow ending on the rect therefore has its head drawn UNDER
# the box, which is zorder 3 against the arrow's 2. Pulling both ends back by slightly
# more than the pad puts the head in clear air.
#
# Done in data units rather than FancyArrowPatch's shrinkA/shrinkB, which are in points:
# these axes differ in physical width by 3x between the two figures, so a point value
# that cleared the pad in one panel would gouge a hole in the other.
BOX_PAD = 0.012
ARROW_GAP = BOX_PAD + 0.002


def arrow(ax, p, q, colour=None, lw=1.1, style='-|>', ls='-', rad=0.0,
          gap=ARROW_GAP, gap_a=None):
    (x0, y0), (x1, y1) = p, q
    dx, dy = x1 - x0, y1 - y0
    dist = (dx * dx + dy * dy) ** 0.5
    if dist > 1e-9:
        ux, uy = dx / dist, dy / dist
        ga = ARROW_GAP if gap_a is None else gap_a
        # Never eat more than 40% of the segment, so short connectors survive.
        ga, gb = min(ga, dist * 0.4), min(gap, dist * 0.4)
        x0, y0 = x0 + ux * ga, y0 + uy * ga
        x1, y1 = x1 - ux * gb, y1 - uy * gb
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=9, linewidth=lw,
                                 color=colour or FAINT, linestyle=ls, zorder=2,
                                 connectionstyle=f'arc3,rad={rad}'))


# arc3's control point is offset in DATA coordinates, then stretched by the axes
# transform -- so one `rad` renders as two different shapes in two panels of
# different aspect. These panels run 1.12 (architecture, near square) and 2.59
# (mechanisms, wide and short), so an arc tuned in one swept far wider in the other.
# Scaling rad by the aspect ratio keeps the rendered curvature the same in both.
REFERENCE_ASPECT = 1.12


def aspect_rad(ax, rad):
    """`rad` corrected so the drawn curve matches REFERENCE_ASPECT's appearance."""
    bb = ax.get_window_extent()
    if bb.height <= 0:
        return rad
    return rad * REFERENCE_ASPECT / (bb.width / bb.height)


def blank(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')


# ─── Mechanism panels, shared by both figures ───────────────────────────────

def mech_screen(ax, m):
    """Fetch, screen and discard inside ONE task, so transient disk is bounded.

    The architectural point is the task boundary, not the screen. A batch is fetched,
    screened and reduced to survivors before the next batch starts, so a rejected genome
    is never written rather than written and deleted. Measured on the RefSeq run: a
    150 GB peak of which the output directory is 30 GB.
    """
    blank(ax)
    ax.add_patch(Rectangle((0.05, 0.18), 0.90, 0.60, facecolor='none',
                           edgecolor=FAINT, linewidth=0.9, linestyle=(0, (4, 3)),
                           zorder=1))
    ax.text(0.50, 0.825, 'one FETCH_RENAME_SCREEN task  ·  × 6,028',
            ha='center', fontsize=8.0, color=FAINT, style='italic')

    box(ax, 0.075, 0.50, 0.145, 0.20, 'batch of 25\naccessions', C_ARTEFACT)
    box(ax, 0.325, 0.50, 0.145, 0.20, 'fetch\n(NCBI datasets)')
    box(ax, 0.575, 0.50, 0.145, 0.20, 'pepM screen\n(DIAMOND)', C_GATE)
    box(ax, 0.800, 0.555, 0.135, 0.09, 'keep', GOOD, edge=GOOD, fs=8.0)
    box(ax, 0.800, 0.285, 0.135, 0.09, 'delete', C_DISCARD, fs=8.0)
    for a, b in ((0.220, 0.325), (0.470, 0.575)):
        arrow(ax, (a, 0.60), (b, 0.60))
    arrow(ax, (0.720, 0.620), (0.800, 0.612), colour=GOOD)
    arrow(ax, (0.648, 0.500), (0.868, 0.378), colour=BEFORE, rad=-0.40)
    ax.text(0.868, 0.678, '0.87%', ha='center', fontsize=8.6, color=GOOD,
            fontweight='bold')
    ax.text(0.868, 0.232, '99.1%', ha='center', fontsize=8.6, color=FAINT)
    ax.text(0.50, 0.075,
            'Measured: the output directory held 30 GB of a 150 GB peak, because the\n'
            '99.1% that fail the screen are never written rather than written and deleted.',
            ha='center', va='center', fontsize=7.6, color=INK, linespacing=1.5)
    ax.set_title('Screening is inside the fetch', pad=9, fontsize=9.5)


def mech_recovery(ax, m):
    """ORF recovery runs UPSTREAM of detection, which is why it changes detection.

    The ordering is the mechanism. BUILD_PROTEIN_POOL -> RECOVER_ORFS -> a GFF handed
    to antiSMASH, so a recovered CDS is not a post-hoc annotation: it can match the
    phosphonate rule itself and extend the cluster core. Treating recovery as a
    cosmetic annotation step is what makes the blind spot invisible.
    """
    blank(ax)
    box(ax, 0.045, 0.620, 0.175, 0.165, 'screened\ngenomes', C_ARTEFACT, fs=8.0)
    box(ax, 0.045, 0.300, 0.175, 0.175, 'BUILD_\nPROTEIN_POOL')
    box(ax, 0.345, 0.300, 0.175, 0.175, 'RECOVER_ORFS\n(miniprot)')
    box(ax, 0.345, 0.020, 0.175, 0.130, 'recovered\nGFF3', C_ARTEFACT, fs=8.0)
    box(ax, 0.590, 0.620, 0.175, 0.165, 'antiSMASH\nphosphonate rule', C_GATE, fs=8.0)
    box(ax, 0.590, 0.300, 0.175, 0.165, 'region\nGenBanks', C_ARTEFACT, fs=8.0)

    arrow(ax, (0.132, 0.620), (0.132, 0.475))
    arrow(ax, (0.220, 0.388), (0.345, 0.388))
    arrow(ax, (0.432, 0.300), (0.432, 0.150))
    # The feedback arc loops OUTWARD, round the right-hand column, rather than
    # cutting back across RECOVER_ORFS. It enters antiSMASH from the right, which
    # also keeps it clear of the region-GenBanks card it would otherwise cross.
    # One continuous patch. Splitting it into a headless curve plus a straight entry
    # segment put the head exactly where it belonged but left a visible break at the
    # join, which is worse than the problem it solved. A single arc3 at this radius
    # clears the region-GenBanks card and lands its head on antiSMASH's right edge;
    # the small end gap is what keeps that head outside the rounded boxstyle.
    # NOT aspect-corrected, deliberately. Clearing the region-GenBanks card is a
    # DATA-space requirement and does not shrink when the curvature is tightened for
    # a wide panel -- correcting it made the arc cut straight through the card and
    # fall short of antiSMASH. One rad that clears in both panels is the honest
    # trade; the two still differ slightly in appearance because the panels do.
    #
    # A wide outward sweep: out past the right-hand column, then back in to
    # antiSMASH's right edge. rad 0.78 carries the apex to x = 0.886, clear of the
    # region-GenBanks card at 0.777.
    #
    # Do not raise rad much above this. Past ~0.87 matplotlib renders the arc in TWO
    # pieces with a visible gap in the middle -- measured by rasterising the patch
    # and counting runs of coloured columns: 1 run at 0.86, 2 runs with a 29 px hole
    # at 0.88. The path itself stays inside the axes, so it is a curvature artefact
    # in the arrowstyle, not clipping.
    #
    # gap=0 at the head is what makes it visible. Earlier versions trimmed 0.004
    # here, which is less than box()'s 0.012 boxstyle pad -- so the tip landed
    # INSIDE the card's visible edge and the fill drew over it. Ending exactly on
    # that edge puts the head in clear air, pointing back into the box.
    arrow(ax, (0.520, 0.085), (0.779, 0.690), colour=ACCENT, rad=0.78,
          gap=0.0, gap_a=0.010)
    arrow(ax, (0.220, 0.702), (0.590, 0.702))
    arrow(ax, (0.677, 0.620), (0.677, 0.465))
    ax.text(0.44, 0.925, 'a recovered CDS can MATCH the rule\nand extend the cluster core',
            ha='center', va='center', fontsize=7.5, color=ACCENT, linespacing=1.5)
    ax.set_title('Recovery feeds detection', pad=9, fontsize=9.5)


def mech_clustering(ax, m):
    """Clustering is a conditional sub-DAG, and the reference pass works on a copy.

    Two separate facts a flow chart usually hides. The partition branch exists because
    BiG-SCAPE is quadratic in memory; the copy exists because a reference inside the GCF
    cutoff would JOIN a family, and every downstream consumer reads the database without
    knowing which records are the dataset.
    """
    blank(ax)
    box(ax, 0.020, 0.595, 0.150, 0.17, 'region\nGenBanks', C_ARTEFACT, fs=8.0)
    ax.text(0.405, 0.912, 'bigscape_partition', ha='center', fontsize=7.6,
            color=FAINT, style='italic')
    box(ax, 0.245, 0.735, 0.340, 0.130, 'false:  BIGSCAPE\n(one job)', fs=8.0)
    box(ax, 0.245, 0.415, 0.340, 0.225,
        'true:\nPARTITION_BGCS →\nBIGSCAPE_PARTITION × N →\nMERGE_BIGSCAPE', fs=7.2)
    box(ax, 0.645, 0.595, 0.180, 0.17, 'BiG-SCAPE\nSQLite', C_ARTEFACT, fs=8.0)
    box(ax, 0.645, 0.290, 0.180, 0.135, 'copy', C_DISCARD, fs=8.0)
    box(ax, 0.645, 0.012, 0.180, 0.130, '+ 17 references', C_GATE, fs=7.9)

    arrow(ax, (0.170, 0.700), (0.245, 0.800), rad=0.12)
    arrow(ax, (0.170, 0.660), (0.245, 0.528), rad=-0.12)
    arrow(ax, (0.585, 0.800), (0.645, 0.700), rad=0.12)
    arrow(ax, (0.585, 0.528), (0.645, 0.660), rad=-0.12)
    arrow(ax, (0.735, 0.595), (0.735, 0.425), colour=ACCENT)
    arrow(ax, (0.735, 0.290), (0.735, 0.142), colour=ACCENT)
    ax.text(0.605, 0.215, 'the published database\nstays dataset-only',
            ha='right', va='center', fontsize=7.3, color=ACCENT, linespacing=1.5)
    ax.set_title('References never touch the original', pad=9, fontsize=9.5)


# ─── Main-body figure ───────────────────────────────────────────────────────

def figure_mechanisms(m, outdir):
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.6))
    for ax, fn, letter in zip(axes, (mech_screen, mech_recovery, mech_clustering),
                              'ABC'):
        fn(ax, m)
        panel_label(ax, letter, dx=-0.055, dy=1.12)
    fig.suptitle('Three mechanisms that decide what the pipeline finds',
                 fontsize=11.5, fontweight='bold', y=0.995)
    fig.text(0.5, 0.008,
             f"Measured on {m['genomes']:,} Enterobacterales genomes "
             f"({m['wall']}, {m['tasks']:,} tasks, 0 failures).",
             ha='center', fontsize=7.6, color=FAINT)
    fig.subplots_adjust(top=0.925, bottom=0.045, hspace=0.46)
    return save(fig, outdir, 'fig5_mechanisms')


# ─── Supplemental figure ───────────────────────────────────────────────────

STAGES = [
    ('DOWNLOAD_GENOMES', 'accessions → screened genomes', 'FETCH_RENAME_SCREEN'),
    ('ANTISMASH_ANALYSIS', 'recovery, then detection',
     'BUILD_PROTEIN_POOL · RECOVER_ORFS · ANTISMASH'),
    ('CLUSTERING  ∥  PHYLOGENY', 'independent; run concurrently',
     'BIGSCAPE · GTDBTK_CLASSIFY'),
    ('BGC_ANALYSIS', 'coupling, branch point, consensus',
     'GCF_CHARACTERISATION · …'),
    ('VISUALIZE_RESULTS', 'the report', 'bgc_report.html'),
]

ARTEFACTS = ['accessions.txt', 'renamed .gbff', 'region .gbk', 'SQLite + .tsv',
             'bgc_report.html']


# Measured CPU per stage, rolled up from the per-process figures in
# enterobacterales_refseq/summary.json: FETCH_RENAME_SCREEN 42.1 CPU-h,
# ANTISMASH 30.5 + RECOVER_ORFS 20.9, BIGSCAPE 0.7 + GTDBTK_CLASSIFY 9.1, and a
# 0.9 CPU-h remainder across every analysis and report process. 104.2 CPU-h total.
#
# One segment per stage box, so the strip can be read straight down from the stage it
# belongs to. The widths are proportional and the boxes are not, so correspondence is
# carried by COLOUR -- each stage box takes a tinted band in its segment's hue -- and by
# a leader dropping from the box to its segment. Aligning the widths instead would mean
# either distorting the shares or drawing a 0.9% stage box nothing could be written in.
CPU_STAGES = [
    ('DOWNLOAD_GENOMES',        40.4, '42.1 CPU-h', '6,028 tasks', '#1b5e7e'),
    ('ANTISMASH_ANALYSIS',      49.3, '51.4 CPU-h', 'recovery 20.0% + detection 29.3%', '#2f7d5d'),
    ('CLUSTERING ∥ PHYLOGENY',   9.5, ' 9.8 CPU-h', 'BiG-SCAPE 0.7% of it', '#7fa8b8'),
    ('BGC_ANALYSIS + report',    0.9, ' 0.9 CPU-h', '', ACCENT),
]


def spine_with_cpu(ax, m):
    """The stage flow, with the CPU strip directly beneath it on the same axes.

    Three channels read downward for each stage: what data survives it, what it is, and
    what artefact it emits -- then the strip, so the cost of a stage sits under the
    stage itself rather than in a separate panel the reader has to hold in mind.
    """
    blank(ax)
    n = len(STAGES)
    w, gap = 0.150, 0.052
    x0 = (1 - (n * w + (n - 1) * gap)) / 2
    survivors = [f"{m['genomes']:,}\ngenomes", f"{m['passed']:,}\npass (0.87%)",
                 f"{m['regions']:,}\nregions", f"{m['gcfs']}\nfamilies", '']
    # stage box -> CPU segment. CLUSTERING and PHYLOGENY share one box and one segment;
    # BGC_ANALYSIS and VISUALIZE_RESULTS share the 0.9% remainder.
    seg_of = [0, 1, 2, 3, 3]
    for i, ((name, sub, inner), surv) in enumerate(zip(STAGES, survivors)):
        x = x0 + i * (w + gap)
        box(ax, x, 0.555, w, 0.225, '')
        # A tinted band keyed to the stage's CPU segment, so the strip below can be
        # read back to its stage by colour where the widths cannot line up.
        ax.add_patch(Rectangle((x, 0.555), w, 0.028,
                               facecolor=CPU_STAGES[seg_of[i]][4], edgecolor='none',
                               zorder=5))
        # Three slots inside the box, each measured to fit, so a long process list
        # cannot run into the neighbouring stage the way it did in the first draft.
        fit_text(ax, x, 0.730, w, 0.045, sub, fs=7.2, colour=FAINT)
        fit_text(ax, x, 0.665, w, 0.058, name, fs=8.4, weight='bold')
        fit_text(ax, x, 0.588, w, 0.070, inner, fs=6.8, colour=INK)
        if surv:
            box(ax, x, 0.835, w, 0.135, surv, C_ARTEFACT, fs=8.0)
        ax.text(x + w / 2, 0.505, ARTEFACTS[i], ha='center', va='top', fontsize=7.2,
                color=FAINT)
        if i:
            arrow(ax, (x - gap, 0.668), (x, 0.668), lw=1.3, colour=AFTER)
    ax.text(x0 - 0.014, 0.900, 'surviving\ndata', ha='right', va='center',
            fontsize=7.4, color=FAINT, linespacing=1.4)
    ax.text(x0 - 0.014, 0.478, 'artefact', ha='right', va='center',
            fontsize=7.4, color=FAINT)

    # ── the strip, directly beneath ────────────────────────────────────────
    span, ys, hs = x0 + n * w + (n - 1) * gap - x0, 0.305, 0.105
    x = x0
    for j, (name, share, cpuh, note, col) in enumerate(CPU_STAGES):
        seg_w = span * share / 100.0
        ax.add_patch(Rectangle((x, ys), seg_w, hs, facecolor=col,
                               edgecolor='white', linewidth=1.0, zorder=3))
        if seg_w > 0.10:
            ax.text(x + seg_w / 2, ys + hs / 2, f'{share}%  ·  {cpuh.strip()}',
                    ha='center', va='center', fontsize=7.6, color='white',
                    fontweight='bold', zorder=4)
            if note:
                ax.text(x + seg_w / 2, ys - 0.048, note, ha='center', fontsize=6.9,
                        color=FAINT, zorder=4)
        else:
            # Two narrow segments sit side by side at the right-hand end and their
            # labels overlapped into an unreadable run of digits. Staggered, and
            # anchored by the end they are nearer so neither runs off the axes.
            drop = 0.100 if j % 2 == 0 else 0.172
            ax.annotate(f'{name.split()[0]}  {share}%',
                        xy=(x + seg_w / 2, ys), xytext=(x + seg_w / 2, ys - drop),
                        fontsize=7.2, color=col, fontweight='bold', ha='right',
                        arrowprops=dict(arrowstyle='-', lw=0.7, color=col,
                                        shrinkA=1, shrinkB=1))
        x += seg_w
    ax.text(x0 - 0.014, ys + hs / 2, 'CPU', ha='right', va='center', fontsize=7.8,
            color=INK, fontweight='bold')
    ax.text(x0 + span + 0.012, ys + hs / 2, '104.2\nCPU-h', ha='left', va='center',
            fontsize=7.2, color=FAINT, linespacing=1.4)
    ax.text(0.5, 0.055,
            'Cost runs opposite to the data: the stage that discards 99.1% of the input '
            'is the most expensive, and BiG-SCAPE —\nquadratic in both time and memory — '
            'is 0.7% of CPU. Segment widths are measured shares; stage boxes are not to '
            'scale.',
            ha='center', va='center', fontsize=7.7, color=INK, linespacing=1.6)


def figure_architecture(m, outdir):
    fig = plt.figure(figsize=(12.6, 9.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.32, 1.0], hspace=0.30, wspace=0.26)

    ax = fig.add_subplot(gs[0, :])
    spine_with_cpu(ax, m)
    panel_label(ax, 'A', dx=-0.012, dy=1.04)

    for col, (fn, letter) in enumerate(zip((mech_screen, mech_recovery,
                                            mech_clustering), 'BCD')):
        axm = fig.add_subplot(gs[1, col])
        fn(axm, m)
        panel_label(axm, letter, dx=-0.055, dy=1.14)

    fig.suptitle('ClusterQuest architecture: stages, artefacts and mechanisms',
                 fontsize=12, fontweight='bold', y=0.975)
    fig.text(0.5, 0.012,
             f"All figures measured on the Enterobacterales RefSeq run — "
             f"{m['genomes']:,} genomes, {m['wall']}, {m['tasks']:,} tasks, 0 failures.",
             ha='center', fontsize=7.6, color=FAINT)
    fig.subplots_adjust(top=0.93, bottom=0.055, left=0.055, right=0.975)
    return save(fig, outdir, 'figS1_architecture')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', default='docs/figures')
    a = ap.parse_args()
    m = measured()
    print(f"funnel: {m['genomes']:,} -> {m['passed']:,} -> {m['regions']:,} "
          f"-> {m['gcfs']} families")
    figure_mechanisms(m, a.outdir)
    figure_architecture(m, a.outdir)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
