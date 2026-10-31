#!/usr/bin/env python3
"""Figures 5 and S1 — what the pipeline does, and the three mechanisms that make it work.

Two figures from one generator, because they share the same facts and must not drift:

  fig5_mechanisms   main body. Three panels, one per mechanism, no stage overview.
  figS1_architecture  supplemental. Stages, data artefacts and the funnel in one view.

Every number is measured on the Enterobacterales GenBank run (2026-10-16, 1,349,504
genomes, 2 d 15 m, 0 failures) and sourced from
docs/comparisons/enterobacterales_genbank/summary.json. Nothing here is illustrative.

The architecture has one defining fact and one surprise, and both figures carry them:

  the funnel    1,349,504 genomes -> 1,530 pass the pepM screen (0.11%) -> 1,476
                regions -> 112 families. 99.89% is discarded before anything
                expensive runs.
  the inversion cost runs opposite to intuition. Fetch-and-screen is 91.2% of CPU over
                53,981 tasks; BiG-SCAPE, the one that is quadratic in both time and
                memory, is 0.07%. Moving and filtering data costs more than analysing
                it at this prevalence, which is why the screen is the architecture
                rather than an optimisation.

These figures were on the RefSeq run (150,690 genomes) until 2026-10-28. The switch
strengthens the thesis rather than restating it: at 8.9x the input and an eighth the
retention, fetch-and-screen goes 40.4% -> 91.2% and antiSMASH 29.3% -> 4.0%. The one
thing lost is peak disk and peak memory, which only the RefSeq run sampled -- the
panel that cites them now says so. See `not_measured` in the summary.

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
MEASURED = ROOT / 'docs' / 'comparisons' / 'enterobacterales_genbank' / 'summary.json'

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
        # Derived, never typed. These two appear in four places between the figures
        # and were hardcoded as 0.87% / 99.1% from the RefSeq run; when the basis
        # moved they were the only numbers that did not follow it.
        'pass_pct': f'{100 * passed / f["screened"]:.2f}%',
        'drop_pct': f'{100 * (1 - passed / f["screened"]):.2f}%',
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
    is never written rather than written and deleted.

    The disk figures here are the RefSeq run's -- a 150 GB peak of which the output
    directory is 30 GB -- and are labelled as such in the panel. They are the only
    measured ones: the GenBank run had no disk sampler, and peak occupancy cannot be
    recovered from a trace, which records bytes moved rather than bytes held. The
    mechanism is identical in both runs; only the scale it was sampled at differs.
    """
    blank(ax)
    ax.add_patch(Rectangle((0.05, 0.18), 0.90, 0.60, facecolor='none',
                           edgecolor=FAINT, linewidth=0.9, linestyle=(0, (4, 3)),
                           zorder=1))
    ax.text(0.50, 0.825, 'one FETCH_RENAME_SCREEN task  ·  × 53,981',
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
    ax.text(0.868, 0.678, m['pass_pct'], ha='center', fontsize=8.6, color=GOOD,
            fontweight='bold')
    ax.text(0.868, 0.232, m['drop_pct'], ha='center', fontsize=8.6, color=FAINT)
    ax.text(0.50, 0.075,
            f'The {m["drop_pct"]} that fail the screen are never written rather than written\n'
            'and deleted. Sampled on the RefSeq run: a 150 GB peak, outdir 30 GB.',
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
# enterobacterales_genbank/summary.json: FETCH_RENAME_SCREEN 1,410.4 CPU-h +
# NCBI_FETCH_METADATA 1.1, ANTISMASH 61.2 + RECOVER_ORFS 33.3, GTDBTK_CLASSIFY 38.2 +
# BIGSCAPE 1.0, and a 1.3 CPU-h remainder across every report process. 1,547.2 total.
#
# Those are COMPLETED and CACHED tasks together. The run was resumed four times, so a
# cached entry carries the realtime of the session that executed it and dropping them
# undercounts the fetch stage by more than half (705.2 CPU-h against 1,547.2).
#
# One segment per stage box, so the strip can be read straight down from the stage it
# belongs to. The widths are proportional and the boxes are not, so correspondence is
# carried by COLOUR -- each stage box takes a tinted band in its segment's hue -- and by
# a leader dropping from the box to its segment. Aligning the widths instead would mean
# either distorting the shares or drawing a 0.9% stage box nothing could be written in.
CPU_STAGES = [
    ('DOWNLOAD_GENOMES',       91.23, '1,411.5 CPU-h', '53,985 tasks', '#1b5e7e'),
    ('ANTISMASH_ANALYSIS',      6.11, '   94.6 CPU-h', 'recovery 2.2% + detection 4.0%', '#2f7d5d'),
    ('CLUSTERING ∥ PHYLOGENY',  2.57, '   39.8 CPU-h', 'BiG-SCAPE 0.07% of it', '#7fa8b8'),
    ('BGC_ANALYSIS + report',   0.09, '    1.3 CPU-h', '', ACCENT),
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
    survivors = [f"{m['genomes']:,}\ngenomes", f"{m['passed']:,}\npass ({m['pass_pct']})",
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
    # Every narrow label is right-aligned on ONE column, set at the centre of the
    # first segment too small to carry its own label. Leaders then fan out to the
    # RIGHT of that column, so a lower label's leader cannot cross the text of the
    # label above it -- which is what happened with each label over its own centre.
    label_x, _x = None, x0
    for _s in CPU_STAGES:
        _w = span * _s[1] / 100.0
        if _w <= 0.10 and label_x is None:
            label_x = _x + _w / 2
        _x += _w
    narrow, labels = 0, []
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
            # Narrow segments crowd into the right-hand end and their labels overlap
            # into an unreadable run of digits. Stagger by HOW MANY narrow segments
            # have been drawn, not by the segment index: on the RefSeq run exactly two
            # were narrow and `j % 2` happened to separate them, but at GenBank shares
            # three are, and that formula put the first and third on the same line --
            # rendering 'ANTISMASH_ANALYBGC_ANALYSIS 0.09%'. Counting cannot collide
            # however many fall below the threshold.
            drop = 0.100 + 0.072 * narrow
            narrow += 1
            # Leader now, label after the loop. One annotate carrying both would
            # draw its own arrow at its own zorder, so the LAST segment's leader
            # painted over the opaque background of a label already placed -- the
            # arrow and the text it must not cross belonged to the same artist.
            # Drawing every line first and every label second orders them properly
            # whatever the share distribution turns out to be.
            ax.annotate('', xy=(x + seg_w / 2, ys), xytext=(label_x, ys - drop),
                        arrowprops=dict(arrowstyle='-', lw=0.7, color=col,
                                        shrinkA=3, shrinkB=1), zorder=4)
            labels.append((label_x, ys - drop, f'{name.split()[0]}  {share}%', col))
        x += seg_w
    for lx, ly, text, col in labels:
        ax.text(lx, ly, text, ha='right', va='center', fontsize=7.2, color=col,
                fontweight='bold', zorder=8,
                bbox=dict(boxstyle='square,pad=0.22', facecolor='white',
                          edgecolor='none'))
    ax.text(x0 - 0.014, ys + hs / 2, 'CPU', ha='right', va='center', fontsize=7.8,
            color=INK, fontweight='bold')
    ax.text(x0 + span + 0.012, ys + hs / 2, '1,547\nCPU-h', ha='left', va='center',
            fontsize=7.2, color=FAINT, linespacing=1.4)
    ax.text(0.5, 0.055,
            f'Cost runs opposite to the data: the stage that discards {m["drop_pct"]} of the '
            'input is the most expensive, and BiG-SCAPE —\nquadratic in both time and '
            'memory — is 0.07% of CPU. Segment widths are measured shares; stage boxes are '
            'not to scale.',
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
             f"All figures measured on the Enterobacterales GenBank run — "
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
