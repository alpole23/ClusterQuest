#!/usr/bin/env python3
"""Figures 5 and S1 — what the pipeline does, and the three mechanisms that make it work.

Two figures from one generator, because they share the same facts and must not drift:

  fig5_mechanisms   main body. Three panels, one per mechanism, no stage overview.
  figS1_architecture  supplemental. Stages, data artefacts and the funnel in one view.

Every number is measured on the Enterobacterales RefSeq run (2026-09-27, 150,690
genomes, 43.1 h, 6,098 tasks, 0 failures) and read at render time from that run's
own records -- summary.json for the funnel, pipeline_trace.tsv.gz for every per-stage
task count, elapsed span and CPU figure, disk_samples.tsv for every peak. Nothing
here is illustrative and nothing is typed in by hand.

WHY THIS RUN AND NOT THE LARGER ONE. These figures sat on the GenBank run (1,349,504
genomes) from 2026-10-28 until 2026-11-02, which is 8.9x the input. That run has no
committed trace and no disk or memory sampler, so it can support a funnel and a CPU
share and nothing else -- and the brief for this figure is real data for EVERY step.
Taking per-stage wall and memory from RefSeq while taking the funnel from GenBank
would have put two runs in one spine, which is the same error class as mixing two CPU
bases. So the whole figure is RefSeq and the footer gives GenBank's scale in words.

The architecture has one defining fact and one surprise, and both figures carry them:

  the funnel    150,690 genomes -> 1,309 pass the pepM screen (0.87%) -> 1,303
                regions -> 72 families. 99.13% is discarded before anything
                expensive runs.
  the three     "Which stage is expensive" has three different right answers, so
  answers       panel A draws three strips rather than picking one. Fetch-and-screen
                is 75.4% of ELAPSED time and 68.2% of RESERVED cores but only 40.4%
                of CPU USED; antiSMASH_ANALYSIS is 12.8% / 20.8% / 49.5%. The fetch
                is network-bound on four allocated cores, so reserving them is most
                of what it does. BiG-SCAPE -- quadratic in both time and memory -- is
                0.7% of CPU used, because the funnel hands it 1,303 of 150,690.

A caveat the figure prints rather than hides: 72 families was measured at
`bigscape_min_bgc_length = 0`, which was the default until 2026-10-30. The shipped
default is 2,000 bp, this run contains exactly one region below it (1,001 bp, one
CDS), and it was never re-clustered. On the GenBank run removing the equivalent
record moved the count 112 -> 126, so this is not a rounding difference.

Usage:
    python scripts/figures/fig5_architecture.py --outdir docs/figures
"""
import argparse
import csv
import gzip
import json
import re
from datetime import datetime
from pathlib import Path

from figure_style import (ACCENT, AFTER, BEFORE, FAINT, GOOD, INK,
                          panel_label, plt, save)
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parent.parent.parent
RUN = ROOT / 'docs' / 'comparisons' / 'enterobacterales_refseq'
MEASURED = RUN / 'summary.json'
TRACE = RUN / 'pipeline_trace.tsv.gz'
SAMPLES = RUN / 'disk_samples.tsv'

# Which processes constitute each displayed stage. Written out rather than derived
# from the trace's subworkflow prefix, for one reason: GET_ANTISMASH_VERSION carries
# the ANTISMASH_ANALYSIS prefix and runs at t+0 beside the fetch, so a prefix rollup
# reports that stage's elapsed span as 38.03 h rather than its true 5.51 h. A stage's
# span has to be taken over the processes that do its work.
STAGE_PROCS = {
    'DOWNLOAD_GENOMES': ['NCBI_FETCH_METADATA', 'PEPM_MAKEDB', 'CREATE_NAME_MAP',
                         'FETCH_RENAME_SCREEN', 'EXTRACT_TAXONOMY'],
    'ANTISMASH_ANALYSIS': ['BUILD_PROTEIN_POOL', 'RECOVER_ORFS', 'ANTISMASH'],
    'CLUSTERING  ∥  PHYLOGENY': ['BIGSCAPE', 'BIGSCAPE_REFERENCES', 'CLUSTERING_STATS',
                                 'EXTRACT_GCF_REPRESENTATIVES', 'GENBANK_TO_FASTA',
                                 'GTDBTK_CLASSIFY', 'MERGE_GTDBTK'],
}
# Everything not named above -- the small analysis and report processes -- rolls into
# a single remainder segment, because a 0.3% stage cannot be drawn as a labelled box.
REMAINDER = 'BGC_ANALYSIS + report'

# Role colours. Deliberately few: a reader should not have to learn a key to follow
# the flow, so hue marks the KIND of thing (compute / artefact / discard) and nothing
# is encoded by colour alone -- every box is also labelled.
C_COMPUTE = '#dce6ee'
C_ARTEFACT = '#f2ead9'
C_DISCARD = '#ebedef'
C_GATE = '#f7ddcb'


def _duration_s(text):
    """Nextflow duration -> seconds. Spellings present: 'Nm Ns', 'Ns', 'Nms', 'Nm'."""
    total = 0.0
    for value, unit in re.findall(r'([\d.]+)\s*(ms|[dhms])', text or ''):
        total += float(value) * {'ms': 0.001, 's': 1, 'm': 60,
                                 'h': 3600, 'd': 86400}[unit]
    return total


def _stamp(text):
    return datetime.strptime(text.strip(), '%Y-%m-%d %H:%M:%S.%f')


def per_stage():
    """Tasks, elapsed span and CPU on BOTH bases, per stage, from the raw trace.

    Returns the stages in execution order plus the run totals. CACHED rows are summed
    alongside COMPLETED ones: a cached entry carries the realtime of the session that
    actually executed it, and dropping them undercounts every resumed stage.
    """
    rows = []
    with gzip.open(TRACE, 'rt') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['status'] not in ('COMPLETED', 'CACHED'):
                continue
            rows.append({
                'proc': r['name'].split(' (')[0].split(':')[-1],
                'rt': _duration_s(r['realtime']),
                'pct': float(r['%cpu'].rstrip('%')) / 100 if r['%cpu'].endswith('%') else 0.0,
                'cpus': int(r['cpus']),
                'start': _stamp(r['start']), 'end': _stamp(r['complete']),
            })
    t0 = min(r['start'] for r in rows)
    named = {p for procs in STAGE_PROCS.values() for p in procs}

    def roll(label, sel):
        return {
            'label': label, 'tasks': len(sel),
            'cpu_used': sum(r['rt'] * r['pct'] for r in sel) / 3600,
            'cpu_reserved': sum(r['rt'] * r['cpus'] for r in sel) / 3600,
            'start_h': (min(r['start'] for r in sel) - t0).total_seconds() / 3600,
            'end_h': (max(r['end'] for r in sel) - t0).total_seconds() / 3600,
        }

    stages = [roll(label, [r for r in rows if r['proc'] in procs])
              for label, procs in STAGE_PROCS.items()]
    stages.append(roll(REMAINDER, [r for r in rows if r['proc'] not in named]))
    for s in stages:
        s['wall_h'] = s['end_h'] - s['start_h']
        s['scattered'] = False
    # The remainder is a grab-bag, not a phase: GET_ANTISMASH_VERSION runs at t+0 and
    # COLLECT_VERSIONS at t+42.93, so its "span" is the whole run and any peak sampled
    # inside that window is the run's peak, not this work's. Only its CPU is meaningful.
    stages[-1].update(scattered=True, wall_h=None)

    return stages, {
        'cpu_used': sum(s['cpu_used'] for s in stages),
        'cpu_reserved': sum(s['cpu_reserved'] for s in stages),
        'tasks': sum(s['tasks'] for s in stages),
        'wall_h': max(s['end_h'] for s in stages),
    }


def peak_by_stage(stages):
    """Peak system memory and work-directory occupancy inside each stage's window.

    These are SYSTEM figures from a 60 s sampler, not per-task figures, and that is
    the only instrument here that can be trusted. Nextflow's own `peak_rss` sums RSS
    across the process tree: it reads 97 GB for GTDB-Tk and 18.7 GB for a BiG-SCAPE
    job that ru_maxrss puts near 2 GB, on a machine with 56 GB of RAM. A figure that
    printed those would be printing an artefact of the instrument.
    """
    samples = []
    with open(SAMPLES) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            samples.append({'h': float(r['elapsed_s']) / 3600,
                            'mem': float(r['mem_used_gb']),
                            'work': float(r['work_gb']),
                            'out': float(r['out_gb'])})
    for s in stages:
        if s.get('scattered'):
            s['peak_mem_gb'] = s['peak_work_gb'] = None
            continue
        win = [x for x in samples if s['start_h'] <= x['h'] < s['end_h']]
        s['peak_mem_gb'] = max((x['mem'] for x in win), default=None)
        s['peak_work_gb'] = max((x['work'] for x in win), default=None)
    return {'mem': max(x['mem'] for x in samples),
            'work': max(x['work'] for x in samples),
            'out': max(x['out'] for x in samples)}


def measured():
    """The run's own numbers, so the figure cannot drift from the record."""
    d = json.loads(MEASURED.read_text())
    f, r = d['found'], d['run']
    passed = int(str(f['passed_screen']).split()[0].replace(',', ''))
    stages, totals = per_stage()
    peaks = peak_by_stage(stages)
    return {
        'genomes': r['genomes'], 'screened': f['screened'], 'passed': passed,
        'regions': f['regions'], 'positive': f['bgc_positive_genomes'],
        'gcfs': f['gcfs'], 'wall': r['wall'], 'tasks': r['tasks'],
        'stages': stages, 'totals': totals, 'peaks': peaks,
        # Panel B frames one FETCH_RENAME_SCREEN task and multiplies by how many ran.
        # It was hardcoded to 53,981, the GenBank count, and did not follow the basis
        # back when the figures moved.
        'fetch_tasks': d['download']['batches'],
        # Derived, never typed. These two appear in four places between the figures,
        # and when the basis last moved they were the only numbers that did not
        # follow it.
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
    ax.text(0.50, 0.825,
            f"one FETCH_RENAME_SCREEN task  ·  × {m['fetch_tasks']:,}",
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
            f"and deleted. Sampled across this run: a {m['peaks']['work']:.0f} GB "
            f"work peak, outdir {m['peaks']['out']:.0f} GB.",
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

# One display box per MEASURED stage, so a strip segment and the box above it are the
# same thing. The previous version had five boxes against four CPU segments and needed
# a seg_of indirection table to map between them; two of the boxes then shared a
# segment, which is exactly the kind of near-correspondence a reader mis-reads.
STAGES = [
    ('DOWNLOAD_GENOMES', 'fetch, rename, screen, discard',
     'FETCH_RENAME_SCREEN', 'renamed .gbff', '#1b5e7e'),
    ('ANTISMASH_ANALYSIS', 'recovery, then detection',
     'BUILD_PROTEIN_POOL · RECOVER_ORFS · ANTISMASH', 'region .gbk', '#2f7d5d'),
    ('CLUSTERING  ∥  PHYLOGENY', 'independent in the DAG',
     'BIGSCAPE · GTDBTK_CLASSIFY', 'SQLite + GTDB .tsv', '#7fa8b8'),
    ('BGC_ANALYSIS + report', 'coupling, branch point, consensus',
     'GCF_CHARACTERISATION · VISUALIZE_RESULTS', 'bgc_report.html', ACCENT),
]


def _strip(ax, stages, x0, span, y, h, key, total, caption, fmt):
    """One proportional strip. Segments too narrow to label get a leader below.

    Three of these are drawn, for elapsed time, CPU used and CPU reserved, because
    they do NOT agree about which stage is the expensive one and a single strip would
    have to pick a winner. Segment order is stage order in all three, so the eye can
    track one stage down the three rows.
    """
    x, narrow = x0, 0
    fig = ax.get_figure()
    fig.canvas.draw()                       # a renderer must exist to measure
    ax_px = ax.get_window_extent().width

    def fits(text, seg, fs):
        """Does `text` fit inside a segment `seg` wide, at font size `fs`?

        Measured, not guessed at from a width threshold. A fixed threshold got both
        ends wrong on one strip: 9.8% was dropped as too narrow while 10.7% was
        labelled and overflowed, because the two differ by a percent and the label
        lengths differ by three characters.
        """
        probe = ax.text(0, -1, text, fontsize=fs, fontweight='bold')
        w = probe.get_window_extent(fig.canvas.get_renderer()).width
        probe.remove()
        return w <= seg * ax_px * 0.88

    for s in stages:
        value = s[key] or 0.0
        seg = span * value / total
        ax.add_patch(Rectangle((x, y), seg, h, facecolor=s['colour'],
                               edgecolor='white', linewidth=1.0, zorder=3))
        pct = f'{100 * value / total:.1f}%'
        full = f'{pct}  ·  {fmt(value)}'
        # Prefer share-and-value; fall back to the share alone, which is what the
        # strip is for; drop the label only when even that will not fit.
        text = full if fits(full, seg, 7.4) else (pct if fits(pct, seg, 7.0) else None)
        if text:
            ax.text(x + seg / 2, y + h / 2, text, ha='center', va='center',
                    fontsize=7.4 if text is full else 7.0, color='white',
                    fontweight='bold', zorder=4)
        else:
            narrow += 1
        x += seg
    if narrow:
        # Named in the caption instead of leadered: with three strips stacked there is
        # no room below a segment for a label that does not land in the strip beneath
        # it, and the bottom strip's landed in the caption text.
        ax.text(x0 + span, y - 0.014,
                f'{narrow} segment too small to label' if narrow == 1
                else f'{narrow} segments too small to label',
                ha='right', va='top', fontsize=6.3, color=FAINT, zorder=8)
    ax.text(x0 - 0.012, y + h / 2, caption, ha='right', va='center', fontsize=7.4,
            color=INK, fontweight='bold', linespacing=1.3)


def spine_with_cpu(ax, m):
    """The stage flow, with everything that was measured about each stage beneath it.

    Six channels read downward for one stage: the data entering it, what it is, the
    artefact it emits, its task count and sampled peaks, then three proportional
    strips -- elapsed, CPU used, CPU reserved.

    Segments under 1% of their strip are drawn but not labelled -- see _strip.

    THE REASON THERE ARE THREE STRIPS. They disagree, and the disagreement is the
    architectural fact. Fetch-and-screen is 75% of elapsed time and 68% of reserved
    cores but only 40% of CPU actually used, because it is network-bound and sits on
    four allocated cores achieving a fraction of one. antiSMASH_ANALYSIS is the
    reverse: 13% of elapsed, 21% of reserved, 50% of CPU used. "Which stage is the
    expensive one" has three different right answers, and an earlier version of this
    figure drew the single strip that happened to support its caption.
    """
    blank(ax)
    stages, totals = m['stages'], m['totals']
    # Pair each measured stage with its display row. Both are in execution order and
    # the assert keeps them that way: the measured labels come from STAGE_PROCS, and
    # a rename there must not silently re-colour the strips.
    assert len(stages) == len(STAGES), 'measured stages and display boxes diverged'
    for s, (name, *_rest, colour) in zip(stages, STAGES):
        assert s['label'].split()[0] == name.split()[0], (s['label'], name)
        s['colour'] = colour

    n = len(STAGES)
    w, gap = 0.170, 0.055
    x0 = (1 - (n * w + (n - 1) * gap)) / 2
    entering = [f"{m['genomes']:,}\ngenomes",
                f"{m['passed']:,} screened\n({m['pass_pct']} of input)",
                f"{m['regions']:,}\nBGC regions",
                f"{m['gcfs']} families *"]

    for i, ((name, sub, inner, artefact, colour), s) in enumerate(zip(STAGES, stages)):
        x = x0 + i * (w + gap)
        box(ax, x, 0.620, w, 0.215, '')
        ax.add_patch(Rectangle((x, 0.620), w, 0.026, facecolor=colour,
                               edgecolor='none', zorder=5))
        fit_text(ax, x, 0.788, w, 0.042, sub, fs=7.2, colour=FAINT)
        fit_text(ax, x, 0.728, w, 0.055, name, fs=8.4, weight='bold')
        fit_text(ax, x, 0.652, w, 0.068, inner, fs=6.8, colour=INK)
        box(ax, x, 0.880, w, 0.115, entering[i], C_ARTEFACT, fs=8.0)
        ax.text(x + w / 2, 0.600, artefact, ha='center', va='top', fontsize=7.2,
                color=FAINT)
        # Task count and sampled peaks. The remainder stage reports neither: its
        # processes are scattered across the whole run, so a window peak would be the
        # run's peak rather than this work's.
        if s.get('scattered'):
            detail = f"{s['tasks']} tasks\nscattered across the run"
        else:
            detail = (f"{s['tasks']:,} tasks · {s['wall_h']:.1f} h\n"
                      f"peak {s['peak_mem_gb']:.0f} GB RAM · "
                      f"{s['peak_work_gb']:.0f} GB work")
        ax.text(x + w / 2, 0.556, detail, ha='center', va='top', fontsize=6.9,
                color=INK, linespacing=1.5)
        if i:
            arrow(ax, (x - gap, 0.728), (x, 0.728), lw=1.3, colour=AFTER)

    ax.text(x0 - 0.012, 0.938, 'data entering', ha='right', va='center',
            fontsize=7.4, color=FAINT)
    ax.text(x0 - 0.012, 0.590, 'artefact', ha='right', va='center',
            fontsize=7.4, color=FAINT)
    ax.text(x0 - 0.012, 0.533, 'measured', ha='right', va='center',
            fontsize=7.4, color=FAINT)

    span, hs = n * w + (n - 1) * gap, 0.070
    rows = [(0.400, 'wall_h', 'elapsed', lambda v: f'{v:.1f} h',
             f"{totals['wall_h']:.1f} h"),
            (0.268, 'cpu_used', 'CPU used  (realtime × %cpu)',
             lambda v: f'{v:,.1f} CPU-h', f"{totals['cpu_used']:,.0f} CPU-h"),
            (0.136, 'cpu_reserved', 'CPU reserved  (realtime × cpus)',
             lambda v: f'{v:,.1f} CPU-h', f"{totals['cpu_reserved']:,.0f} CPU-h")]
    for y, key, caption, fmt, total_text in rows:
        _strip(ax, stages, x0, span, y, hs, key, totals[key], caption, fmt)
        ax.text(x0 + span + 0.010, y + hs / 2, total_text, ha='left', va='center',
                fontsize=7.2, color=FAINT)

    ax.text(0.5, 0.048,
            'The three strips disagree, and that is the point: fetch-and-screen is '
            '75% of ELAPSED time and 68% of RESERVED cores but 40% of CPU USED — it '
            'is network-bound on four\nallocated cores. Peak work disk is 150 GB and '
            '95 GB of it is pplacer scratch inside the single GTDB-Tk task, set by '
            "GTDB's reference data rather than by how many genomes reach it.\n"
            'The report stages are under 1% of all three strips and are drawn '
            'unlabelled.   * 72 families was measured at bigscape_min_bgc_length = 0, '
            'the default at the time; the shipped default is now 2,000 bp, this run '
            'contains one 1,001 bp single-CDS region, and it was never re-clustered.',
            ha='center', va='center', fontsize=7.4, color=INK, linespacing=1.7)


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
             f"{m['genomes']:,} genomes, {m['wall']}, {m['tasks']:,} tasks, 0 failures. "
             f"The same pipeline ran 1,349,504 GenBank genomes in 2 d 15 m; that run "
             f"had no disk or memory sampler and no committed trace, so it cannot "
             f"carry per-stage figures.",
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
