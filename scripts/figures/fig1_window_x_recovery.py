#!/usr/bin/env python3
"""Figure 1 — two blind spots in BGC detection, and why neither alone is enough.

antiSMASH reports a rule core plus a fixed flank, over whatever genes the
deposit happened to annotate. Two things therefore go missing, for unrelated
reasons:

  * the flank is too narrow for this chemistry -- the phosphonate rule ships
    NEIGHBOURHOOD 5, and the HiVir cluster runs past it
  * genes the submitter never annotated are invisible, because antiSMASH runs
    gene finding ONLY on records with zero CDS features

These are drawn as a 2x2 rather than a sequence, because **they are not
independent**. Recovered genes can match the detection rule, which extends the
core, which moves the flank: region 2's 3' boundary shifts 1,566 bp from
recovery alone, at BOTH neighbourhood settings. A figure ordering one change
"before" the other would assert an independence the data does not have.

The two clusters in this one genome also need opposite fixes, which is the
argument for doing both: pantaphos gains more from the window, the
phosphonolipid gains more from recovery.

Usage:
    python scripts/figures/fig1_window_x_recovery.py --outdir docs/figures
"""
import argparse
from pathlib import Path

from figure_style import ACCENT, FAINT, GOOD, INK, panel_label, plt, save
from fig1_orf_recovery import (CATEGORY_COLOR, CATEGORY_LABEL, gene_category,
                               is_recovered, label_for, legend_handles,
                               load_transferred)

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent.parent

# Genome coordinates read from each run's region_tabulation.tsv, so the four
# cells can be drawn on one axis. The GBKs are re-based to 0 individually and
# cannot be overlaid without this.
CELLS = [
    ('5 kb  ·  deposited only', 'results/antismash_results/Erwiniaceae_pre_recovery',
     {'region001': (796909, 810246), 'region002': (2688754, 2712709)}),
    ('5 kb  ·  + ORF recovery', 'results/antismash_results/Erwiniaceae',
     {'region001': (796909, 810246), 'region002': (2688754, 2714275)}),
    ('10 kb  ·  deposited only',
     'results_fig1_pan10_norec/antismash_results/Pantoea_10kb',
     {'region001': (791909, 815246), 'region002': (2683754, 2717709)}),
    ('10 kb  ·  + ORF recovery',
     'results_fig1_pan10_rec/antismash_results/Pantoea_10kb',
     {'region001': (791909, 815246), 'region002': (2683754, 2719275)}),
]

GENOME = 'Pantoea_ananatis_LMG_5342'
RECORD = 'HE617160.1'

REGIONS = [
    ('region001', 'pantaphos / HiVir  —  gains most from the WIDER WINDOW',
     (801909, 814435)),      # pepM start to the documented cluster end
    ('region002', 'phosphonolipid  —  gains most from ORF RECOVERY', None),
]


def read_cells(region, transferred):
    """One entry per 2x2 cell: (label, start, end, genes in genome coordinates)."""
    out = []
    for label, base, coords in CELLS:
        path = ROOT / base / GENOME / f'{RECORD}.{region}.gbk'
        start, end = coords[region]
        if not path.exists():
            print(f'  missing: {path}')
            continue
        rec = next(SeqIO.parse(str(path), 'genbank'))
        genes = []
        for f in rec.features:
            if f.type != 'CDS':
                continue
            tag = f.qualifiers.get('locus_tag', [''])[0]
            name = transferred.get(tag)
            cat = gene_category(f)
            lab = label_for(f)
            if name:
                lab = f'{name[:30]}  †'
                if cat == 'other':
                    cat = 'transferred'
            genes.append({
                'start': start + int(f.location.start),
                'end': start + int(f.location.end),
                'strand': 1 if f.location.strand in (None, 1) else -1,
                'recovered': is_recovered(f),
                'category': cat, 'label': lab, 'tag': tag,
            })
        out.append((label, start, end, sorted(genes, key=lambda g: g['start'])))
    return out


def draw_panel(ax, region, title, cluster_extent, transferred):
    cells = read_cells(region, transferred)
    lo = min(c[1] for c in cells)
    hi = max(c[2] for c in cells)
    span = hi - lo

    # Genes present only in the richest cell are the ones the other three miss.
    richest = {g['tag'] for g in cells[-1][3]}
    baseline = {g['tag'] for g in cells[0][3]}

    for row, (label, start, end, genes) in enumerate(cells):
        y = len(cells) - row - 1
        # The captured interval, so a narrower window reads as a shorter band.
        ax.add_patch(plt.Rectangle((start, y - 0.30), end - start, 0.60,
                                   facecolor='#f2f4f7', edgecolor='#dfe4ea',
                                   lw=0.7, zorder=0))
        for g in genes:
            w = g['end'] - g['start']
            head = min(w * 0.30, span * 0.006)
            body = w - head
            col = CATEGORY_COLOR[g['category']]
            edge = INK if g['recovered'] else '#ffffff'
            lw = 1.4 if g['recovered'] else 0.5
            if g['strand'] >= 0:
                v = [(g['start'], y - 0.16), (g['start'] + body, y - 0.16),
                     (g['end'], y), (g['start'] + body, y + 0.16),
                     (g['start'], y + 0.16)]
            else:
                v = [(g['end'], y - 0.16), (g['end'] - body, y - 0.16),
                     (g['start'], y), (g['end'] - body, y + 0.16),
                     (g['end'], y + 0.16)]
            ax.add_patch(plt.Polygon(v, closed=True, facecolor=col,
                                     edgecolor=edge, lw=lw, zorder=3))
        gained = len(genes) - len(cells[0][3])
        ax.text(lo - span * 0.012, y, label, ha='right', va='center',
                fontsize=8.2,
                fontweight='bold' if row == len(cells) - 1 else 'normal',
                color=INK if row == len(cells) - 1 else FAINT)
        ax.text(hi + span * 0.012, y, f'{len(genes)} genes'
                + (f'   +{gained}' if gained else '   baseline'),
                ha='left', va='center', fontsize=8.2,
                color=GOOD if gained else FAINT,
                fontweight='bold' if gained else 'normal')

    if cluster_extent:
        cs, ce = cluster_extent
        ax.plot([cs, ce], [-0.75, -0.75], color=ACCENT, lw=2.2,
                solid_capstyle='butt')
        ax.text((cs + ce) / 2, -0.98, 'HiVir cluster', ha='center', va='top',
                fontsize=7.6, color=ACCENT, fontweight='bold')
        # The 5 kb window stops short of the cluster; say so where it happens.
        five_end = CELLS[0][2][region][1]
        if ce > five_end:
            # Above the top track, not on it: at y=3.0 this drew straight
            # through the "5 kb, deposited only" row and its count.
            ax.annotate('', xy=(ce, 3.48), xytext=(five_end, 3.48),
                        arrowprops=dict(arrowstyle='<->', color=ACCENT, lw=1.0))
            ax.text((five_end + ce) / 2, 3.60,
                    f'{(ce - five_end) / 1000:.1f} kb of cluster outside the 5 kb window',
                    ha='center', fontsize=7.4, color=ACCENT, fontweight='bold')

    missed = len(richest - baseline)
    ax.set_title(f'{title}      ({missed} genes invisible at 5 kb without recovery)',
                 loc='left', fontsize=9.5, pad=8)
    ax.set_xlim(lo - span * 0.20, hi + span * 0.16)
    ax.set_ylim(-1.35, 3.95)
    ax.set_yticks([])
    ticks = [lo, (lo + hi) / 2, hi]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t / 1000:,.0f} kb' for t in ticks], fontsize=7.6)
    for side in ('left', 'top', 'right'):
        ax.spines[side].set_visible(False)
    return cells


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='docs/figures')
    args = ap.parse_args()

    transferred = load_transferred(
        ROOT / 'results_fig1_pan10_rec/main_analysis_results/Pantoea_10kb/'
        'annotation_transfer/gcf_annotation_transfer.tsv')

    fig, axes = plt.subplots(2, 1, figsize=(12.2, 7.4))
    used = []
    for ax, (region, title, extent) in zip(axes, REGIONS):
        cells = draw_panel(ax, region, title, extent, transferred)
        counts = ' -> '.join(str(len(c[3])) for c in cells)
        print(f'{region}: {counts} genes across the 2x2')
        for c in cells:
            for g in c[3]:
                if g['category'] not in used:
                    used.append(g['category'])
    for ax, letter in zip(axes, 'AB'):
        panel_label(ax, letter, dx=-0.155, dy=1.20)

    order = [c for c in CATEGORY_COLOR if c in used]
    axes[-1].legend(handles=legend_handles(order), loc='upper center',
                    bbox_to_anchor=(0.5, -0.20), ncol=4, fontsize=8)
    fig.suptitle('Two independent blind spots, and they compound: '
                 'a too-narrow window and unannotated genes',
                 fontsize=11.5, fontweight='bold', y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, args.outdir, 'fig1_window_x_recovery')


if __name__ == '__main__':
    main()
