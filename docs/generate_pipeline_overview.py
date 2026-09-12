#!/usr/bin/env python3
"""Generate the pipeline overview figure for ClusterQuest.

Run from the repository root:  python docs/generate_pipeline_overview.py

Layout notes, since they encode real facts about the workflow:
  * PEPM_PRESCREEN sits between download and antiSMASH. It is the stage that
    makes order-scale runs tractable, so it gets its own box rather than being
    folded into ANTISMASH_ANALYSIS.
  * CLUSTERING and PHYLOGENY are drawn side by side because they are
    independent — both consume antiSMASH output and neither feeds the other, so
    Nextflow runs them concurrently. The old figure drew them as a chain.
  * Reference databases are named in each stage's caption rather than given a
    column of cylinders; the pinned version is part of the caption because two
    runs on different database versions are not comparable.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

# ─── Palette ────────────────────────────────────────────────────────────────
# Text colours are chosen against their own fill: the output boxes previously
# used near-white on gold (~1.8:1), which is below any usable contrast floor.
C = {
    'genome':    '#6baed6',   # blue    — genome retrieval
    'prescreen': '#f0a878',   # light orange — pepM screen (gates BGC detection)
    'bgc':       '#e6854a',   # orange  — BGC detection
    'cluster':   '#74c476',   # green   — GCF clustering
    'phylo':     '#9e7fc0',   # purple  — phylogenetics
    'score':     '#4fa8a0',   # teal    — analysis & scoring
    'viz':       '#7a9eb8',   # slate   — visualization
    'output':    '#f0c75e',   # gold    — published output
    'ink':       '#2b2b2b',
    'muted':     '#4a4a4a',
    'flow':      '#2d8a2d',
    'reuse':     '#c0392b',
}

fig, ax = plt.subplots(1, 1, figsize=(12, 10))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

ax.text(50, 99, 'ClusterQuest', fontsize=24, fontweight='bold',
        ha='center', va='top', color=C['ink'])
ax.text(50, 94, 'Phosphonate BGC discovery at order scale',
        fontsize=12, ha='center', va='top', style='italic', color=C['muted'])


def stage(x, y, w, h, color, title, caption, note=None, title_size=13):
    """A workflow stage. `caption` says what it does; `note` is the scale lever."""
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle='round,pad=0.02,rounding_size=0.5',
        facecolor=color, edgecolor='#333333', linewidth=1.8))
    ax.text(x, y + h * 0.18, title, fontsize=title_size, fontweight='bold',
            ha='center', va='center', color=C['ink'])
    ax.text(x, y - h * 0.16, caption, fontsize=9.5, ha='center', va='center',
            style='italic', color=C['muted'])
    if note:
        ax.text(x, y - h * 0.38, note, fontsize=8, ha='center', va='center',
                color=C['muted'])


def output(x, y, w, h, label, sub=None):
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle='round,pad=0.02,rounding_size=0.3',
        facecolor=C['output'], edgecolor='#333333', linewidth=1.4))
    dy = h * 0.17 if sub else 0
    ax.text(x, y + dy, label, fontsize=9.5, ha='center', va='center',
            fontweight='bold', color=C['ink'])
    if sub:
        ax.text(x, y - h * 0.22, sub, fontsize=8, ha='center', va='center',
                color=C['muted'])


def arrow(start, end, color=None, style='-', lw=2.0, rad=None):
    kw = dict(arrowstyle='->', color=color or C['flow'], lw=lw, linestyle=style,
              shrinkA=0, shrinkB=0)
    if rad is not None:
        kw['connectionstyle'] = f'arc3,rad={rad}'
    ax.annotate('', xy=end, xytext=start, arrowprops=kw)


# ─── Main spine ─────────────────────────────────────────────────────────────
SX, SW, SH = 45, 40, 10          # spine x, width, height
Y = dict(dl=85, screen=71, as_=57, split=41, score=25, viz=10)

stage(SX, Y['dl'], SW, SH, C['genome'], 'DOWNLOAD_GENOMES',
      'NCBI Datasets  ·  taxdump 2026-08-01')

stage(SX, Y['screen'], SW, SH, C['prescreen'], 'PEPM_PRESCREEN',
      'DIAMOND vs 7 curated PEP mutase references',
      '2,771 → 306 genomes at 298/298 sensitivity (Erwiniaceae)')

stage(SX, Y['as_'], SW, SH, C['bgc'], 'ANTISMASH_ANALYSIS',
      'phosphonate rule  ·  KnownClusterBlast vs MIBiG',
      'batched — 50 genomes per task')

# Parallel branch: independent of one another, so drawn side by side.
BW, BH = 19, 10
PHYLO_X, CLUST_X = 34, 56       # pair centred on SX so the fan-out stays symmetric
stage(PHYLO_X, Y['split'], BW, BH, C['phylo'], 'PHYLOGENY',
      'GTDB-Tk  ·  GTDB r226', 'sharded — 5,000 genomes', title_size=12)
stage(CLUST_X, Y['split'], BW, BH, C['cluster'], 'CLUSTERING',
      'BiG-SCAPE  ·  Pfam 38.2', 'optional pepM partitioning', title_size=12)

stage(SX, Y['score'], SW, SH, C['score'], 'GCF ANALYSIS  +  NOVELTY_SCORE',
      'coupling enzyme classes  ·  pepM all-by-all',
      'priority = distance × evidence', title_size=12)

stage(SX, Y['viz'], SW, SH, C['viz'], 'VISUALIZE_RESULTS',
      'interactive report, deterministic output')

# ─── Published output ───────────────────────────────────────────────────────
OX, OW, OH = 84, 25, 7          # right-hand output column
LX, LW = 11, 20                  # left-hand output (the parallel row's other branch)
output(OX, Y['dl'], OW, OH, 'Renamed genomes', 'ncbi_genomes/')
output(OX, Y['screen'], OW, OH, 'Screening verdicts', 'prescreen_results/')
output(OX, Y['as_'], OW, OH, 'BGC regions', 'antismash_results/')
output(OX, Y['split'], OW, OH, 'Gene cluster families', 'bigscape_results/')
output(LX, Y['split'], LW, OH, 'Placement', 'gtdbtk_results/')
output(OX, Y['score'], OW, OH, 'Ranked families', 'novelty_ranking.tsv')
output(OX, Y['viz'], OW, OH, 'Report', 'bgc_report.html')

# ─── Flow ───────────────────────────────────────────────────────────────────
spine_top = SX
for a, b in [('dl', 'screen'), ('screen', 'as_')]:
    arrow((spine_top, Y[a] - SH / 2), (spine_top, Y[b] + SH / 2))

# fan out to the two independent branches, then back in
arrow((SX, Y['as_'] - SH / 2), (PHYLO_X, Y['split'] + BH / 2), rad=0.15)
arrow((SX, Y['as_'] - SH / 2), (CLUST_X, Y['split'] + BH / 2), rad=-0.15)
arrow((PHYLO_X, Y['split'] - BH / 2), (SX, Y['score'] + SH / 2), rad=-0.15)
arrow((CLUST_X, Y['split'] - BH / 2), (SX, Y['score'] + SH / 2), rad=0.15)
arrow((SX, Y['score'] - SH / 2), (SX, Y['viz'] + SH / 2))

# stage → published output
for k in Y:
    if k == 'split':
        arrow((CLUST_X + BW / 2, Y[k]), (OX - OW / 2, Y[k]), lw=1.6)
        arrow((PHYLO_X - BW / 2, Y[k]), (LX + LW / 2, Y[k]), lw=1.6)
        continue
    arrow((SX + SW / 2, Y[k]), (OX - OW / 2, Y[k]), lw=1.6)

# ─── Cross-taxon reuse (published results feeding a later run) ───────────────
arrow((OX - OW / 2, Y['as_'] - 3.2), (SX + SW / 2, Y['as_'] - 3.2),
      color=C['reuse'], style='--', lw=2.2, rad=-0.3)
arrow((LX + LW / 2, Y['split'] - 3.2), (PHYLO_X - BW / 2, Y['split'] - 3.2),
      color=C['reuse'], style='--', lw=2.2, rad=0.3)

# ─── Legend ─────────────────────────────────────────────────────────────────
handles = [
    mpatches.Patch(facecolor=C['genome'], edgecolor='#333', label='Genome retrieval'),
    mpatches.Patch(facecolor=C['prescreen'], edgecolor='#333', label='pepM pre-screen'),
    mpatches.Patch(facecolor=C['bgc'], edgecolor='#333', label='BGC detection'),
    mpatches.Patch(facecolor=C['cluster'], edgecolor='#333', label='GCF clustering'),
    mpatches.Patch(facecolor=C['phylo'], edgecolor='#333', label='Phylogenetics'),
    mpatches.Patch(facecolor=C['score'], edgecolor='#333', label='Analysis & scoring'),
    mpatches.Patch(facecolor=C['viz'], edgecolor='#333', label='Visualization'),
    Line2D([0], [0], color=C['reuse'], linestyle='--', linewidth=2,
           label='Cross-taxon reuse'),
]
ax.legend(handles=handles, loc='lower center', ncol=4, fontsize=10,
          frameon=False, bbox_to_anchor=(0.5, -0.035),
          handlelength=2.2, handleheight=1.3, columnspacing=2.0)

plt.savefig('docs/pipeline_overview.png', dpi=160, facecolor='white',
            bbox_inches='tight', pad_inches=0.35)
print('Pipeline overview figure saved to docs/pipeline_overview.png')
