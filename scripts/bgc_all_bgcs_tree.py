#!/usr/bin/env python3
"""
Standalone NJ tree of ALL phosphonate BGCs using BiG-SCAPE pairwise distances.

Leaves (one per region-level BGC) are annotated with coupling enzyme class color.

BiG-SCAPE stores all pairwise distances (not just those below the cutoff),
so the full distance matrix is available for all BGCs.

Usage:
    python scripts/bgc_all_bgcs_tree.py \\
        --db      work/57/.../Pantoea.db \\
        --coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \\
        --outdir  results/bgc_trees/Pantoea \\
        [--cutoff 0.3] \\
        [--layout circular|linear]
"""

import argparse
import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import COUPLING_COLORS, COUPLING_ORDER, load_coupling_classes
from utils.tree_building import build_nj_tree as _build_nj_tree


# ─── Data loading ─────────────────────────────────────────────────────────────


def load_bgc_data(conn, cutoff):
    """
    Return:
      bgc_ids   : sorted list of region-level bgc_record ids
      id_to_meta: {record_id: {genome, gbk_basename, gcf_id, organism}}
      distances : {(id_a, id_b): distance}  id_a < id_b
    """
    cur = conn.cursor()

    # All region-level records with GCF assignment
    cur.execute("""
        SELECT br.id, g.path, g.organism, rf.family_id
        FROM bgc_record br
        JOIN gbk g ON g.id = br.gbk_id
        LEFT JOIN bgc_record_family rf ON rf.record_id = br.id
        WHERE br.record_type = 'region'
        ORDER BY br.id
    """)
    rows = cur.fetchall()

    id_to_meta = {}
    for rec_id, path, organism, family_id in rows:
        genome      = os.path.basename(os.path.dirname(path))
        gbk_base    = os.path.splitext(os.path.basename(path))[0]
        id_to_meta[rec_id] = {
            'genome':      genome,
            'gbk_basename': gbk_base,
            'gcf_id':      family_id,   # None for singletons
            'organism':    organism or '',
        }

    bgc_ids = sorted(id_to_meta)

    # All pairwise distances between region-level BGCs
    id_set = set(bgc_ids)
    cur.execute("""
        SELECT d.record_a_id, d.record_b_id, d.distance
        FROM distance d
        JOIN bgc_record br_a ON br_a.id = d.record_a_id
        JOIN bgc_record br_b ON br_b.id = d.record_b_id
        WHERE br_a.record_type = 'region' AND br_b.record_type = 'region'
    """)
    distances = {}
    for a, b, dist in cur.fetchall():
        key = (min(a, b), max(a, b))
        distances[key] = dist

    return bgc_ids, id_to_meta, distances


# ─── Tree building ─────────────────────────────────────────────────────────────

def build_nj_tree(bgc_ids, distances):
    """Build NJ tree from full pairwise distance matrix."""
    n      = len(bgc_ids)
    labels = [str(i) for i in bgc_ids]
    dm_rows = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                a, b = bgc_ids[j], bgc_ids[i]
                key  = (min(a, b), max(a, b))
                row.append(distances.get(key, 1.0))
        dm_rows.append(row)
    print(f'  Building NJ tree for {n} BGCs...')
    return _build_nj_tree(labels, dm_rows)


# ─── Layout helpers (linear) ──────────────────────────────────────────────────

def assign_layout(clade, counter, depth=0):
    """Assign ._x (leaf y-position) and ._depth to every node."""
    clade._depth = depth
    if clade.is_terminal():
        clade._x = counter[0]
        counter[0] += 1
        return
    for child in clade.clades:
        assign_layout(child, counter, depth + 1)
    clade._x = sum(c._x for c in clade.clades) / len(clade.clades)


def max_depth(clade):
    if clade.is_terminal():
        return clade._depth
    return max(max_depth(c) for c in clade.clades)


def draw_cladogram(ax, clade, color='#333333', lw=0.6):
    """Root on left, leaves on right."""
    if clade.is_terminal():
        return
    x_node   = clade._depth
    child_ys = [c._x for c in clade.clades]
    ax.plot([x_node, x_node], [min(child_ys), max(child_ys)],
            color=color, lw=lw, solid_capstyle='round')
    for child in clade.clades:
        ax.plot([x_node, child._depth], [child._x, child._x],
                color=color, lw=lw, solid_capstyle='round')
        draw_cladogram(ax, child, color, lw)


# ─── Layout helpers (circular) ────────────────────────────────────────────────

def assign_circular_layout(clade, leaf_angles, counter, depth=0):
    """Assign ._angle and ._depth to every node for a fan/circular tree."""
    clade._depth = depth
    if clade.is_terminal():
        clade._angle = leaf_angles[counter[0]]
        counter[0] += 1
        return
    for child in clade.clades:
        assign_circular_layout(child, leaf_angles, counter, depth + 1)
    child_angles = [c._angle for c in clade.clades]
    clade._angle = (min(child_angles) + max(child_angles)) / 2


def max_depth_circ(clade):
    if clade.is_terminal():
        return clade._depth
    return max(max_depth_circ(c) for c in clade.clades)


def draw_circular_cladogram(ax, clade, md, color='#888888', lw=0.35):
    """Draw fan-tree branches: arcs at parent radius + radial arms to children."""
    if clade.is_terminal():
        return
    r_n = clade._depth / md
    child_angles = [c._angle for c in clade.clades]
    a_min, a_max = min(child_angles), max(child_angles)

    # Arc at parent radius spanning all children
    n_pts = max(3, int((a_max - a_min) * 100) + 2)
    arc_θ = np.linspace(a_min, a_max, n_pts)
    ax.plot(arc_θ, np.full(n_pts, r_n), color=color, lw=lw, solid_capstyle='butt')

    # Radial arm from parent radius to each child radius
    for child in clade.clades:
        r_c = child._depth / md
        ax.plot([child._angle, child._angle], [r_n, r_c],
                color=color, lw=lw, solid_capstyle='butt')
        draw_circular_cladogram(ax, child, md, color, lw)


# ─── Figure ───────────────────────────────────────────────────────────────────

def plot_tree(tree, bgc_ids, id_to_meta, coupling_classes, outdir):
    # Assign layout
    assign_layout(tree.root, [0])
    md = max_depth(tree.root)

    # Leaf order from tree traversal
    leaf_order = [int(clade.name) for clade in tree.get_terminals()]
    n_leaves   = len(leaf_order)

    # Unique GCF IDs (sorted by size desc, singletons last)
    gcf_counts  = Counter(id_to_meta[i]['gcf_id'] for i in bgc_ids
                          if id_to_meta[i]['gcf_id'] is not None)
    sorted_gcfs = [gid for gid, _ in gcf_counts.most_common()]

    # ── Figure dimensions ───────────────────────────────────────────────────────
    cell_h    = 0.12    # inches per leaf
    tree_w    = 3.0     # tree panel
    strip_w   = 0.18    # width of coupling enzyme strip
    gcf_lbl_w = 0.55    # width of GCF text column
    gap       = 0.05    # gap between panels
    label_w   = 3.2     # genome label area
    left_pad  = 0.2
    right_pad = 1.6     # legend (no GCF color legend needed)

    fig_h = n_leaves * cell_h + 0.5
    fig_w = left_pad + tree_w + gap + strip_w + gap + gcf_lbl_w + gap + label_w + right_pad

    fig = plt.figure(figsize=(fig_w, fig_h))

    f = lambda x: x / fig_w
    heat_bot  = 0.1 / fig_h
    heat_h    = (n_leaves * cell_h) / fig_h
    tree_left = left_pad / fig_w

    ax_tree = fig.add_axes([tree_left,                                        heat_bot, f(tree_w),     heat_h])
    ax_coup = fig.add_axes([tree_left + f(tree_w + gap),                      heat_bot, f(strip_w),    heat_h])
    ax_gcf  = fig.add_axes([tree_left + f(tree_w + gap + strip_w + gap),      heat_bot, f(gcf_lbl_w),  heat_h])
    ax_lbl  = fig.add_axes([tree_left + f(tree_w + gap + strip_w + gap + gcf_lbl_w + gap), heat_bot, f(label_w), heat_h])

    # ── Draw cladogram ──────────────────────────────────────────────────────────
    draw_cladogram(ax_tree, tree.root)
    ax_tree.set_xlim(-0.2, md + 0.5)
    ax_tree.set_ylim(-0.5, n_leaves - 0.5)
    ax_tree.axis('off')

    # ── Draw strips and labels ──────────────────────────────────────────────────
    for plot_idx, rec_id in enumerate(leaf_order):
        meta     = id_to_meta[rec_id]
        gbk_base = meta['gbk_basename']
        gcf_id   = meta['gcf_id']
        genome   = meta['genome'].replace('_', ' ')
        organism = meta['organism']

        # Coupling enzyme class
        cls        = coupling_classes.get(gbk_base, 'Unknown')
        coup_color = COUPLING_COLORS.get(cls, '#aaaaaa')

        y = plot_idx

        # Coupling strip
        ax_coup.barh(y, 1, height=0.85, color=coup_color, left=0)

        # GCF text label
        gcf_text = f'GCF-{gcf_id}' if gcf_id is not None else '—'
        ax_gcf.text(0.5, y, gcf_text, va='center', ha='center',
                    fontsize=5, family='monospace', color='#222222')

        # Genome label: species name + antiSMASH region ID
        label = f'{genome.replace("_", " ")}  {gbk_base}'
        ax_lbl.text(0.02, y, label, va='center', ha='left',
                    fontsize=5.5, family='monospace', color='#222222')

    for ax in (ax_coup, ax_gcf, ax_lbl):
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, n_leaves - 0.5)
        ax.axis('off')

    # ── Strip headers ───────────────────────────────────────────────────────────
    header_y = 1.0 - (0.15 / fig_h)
    for ax, label in [(ax_coup, 'CE'), (ax_gcf, 'GCF')]:
        ax.text(0.5, 1.005, label, transform=ax.transAxes,
                ha='center', va='bottom', fontsize=6, fontweight='bold')

    # ── Legends ────────────────────────────────────────────────────────────────
    coup_patches = [mpatches.Patch(color=COUPLING_COLORS[c], label=c)
                    for c in COUPLING_ORDER if c in coupling_classes.values()]

    legend_x = tree_left + f(tree_w + gap + strip_w + gap + gcf_lbl_w + gap + label_w + 0.15)
    fig.legend(handles=coup_patches, bbox_to_anchor=(legend_x, 0.65),
               loc='upper left', bbox_transform=fig.transFigure,
               fontsize=6.5, title='Coupling enzyme', title_fontsize=7,
               framealpha=0.9)

    fig.suptitle('Phosphonate BGC biosynthetic tree — all BGCs\n'
                 f'NJ tree from BiG-SCAPE pairwise distances  |  '
                 f'{n_leaves} BGCs  |  {len(sorted_gcfs)} GCFs',
                 fontsize=9, y=0.995, va='top')

    # ── Save ───────────────────────────────────────────────────────────────────
    os.makedirs(outdir, exist_ok=True)
    out_png = os.path.join(outdir, 'all_bgcs_biosynthetic_tree.png')
    out_svg = os.path.join(outdir, 'all_bgcs_biosynthetic_tree.svg')
    fig.savefig(out_png, dpi=180, bbox_inches='tight')
    fig.savefig(out_svg,           bbox_inches='tight')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_svg}')
    plt.close(fig)


# ─── Circular figure ──────────────────────────────────────────────────────────

def plot_circular_tree(tree, bgc_ids, id_to_meta, coupling_classes, outdir):
    """Draw a circular (fan) cladogram colored by coupling enzyme class."""
    leaf_order = [int(clade.name) for clade in tree.get_terminals()]
    n = len(leaf_order)

    # Leaf angles: 0 = top (North), clockwise, with a small gap so first and
    # last leaves don't overlap. theta_direction=-1 means CW in polar axes.
    gap_frac = 0.015
    leaf_angles = np.linspace(0, 2 * np.pi * (1 - gap_frac), n)

    assign_circular_layout(tree.root, leaf_angles, [0])
    md = max_depth_circ(tree.root)

    fig_size = 18
    fig, ax = plt.subplots(figsize=(fig_size, fig_size),
                           subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)   # clockwise

    # Tree branches
    draw_circular_cladogram(ax, tree.root, md)

    # Leaf dots (coupling class) and labels
    for i, rec_id in enumerate(leaf_order):
        meta     = id_to_meta[rec_id]
        gbk_base = meta['gbk_basename']
        gcf_id   = meta['gcf_id']
        genome   = meta['genome']
        θ = leaf_angles[i]

        cls   = coupling_classes.get(gbk_base, 'Unknown')
        color = COUPLING_COLORS.get(cls, '#aaaaaa')

        # Colored dot at leaf tip
        ax.scatter([θ], [1.0], s=10, c=[color], zorder=5, linewidths=0)

        # Label: compute readable rotation
        deg     = np.degrees(θ) % 360
        mpl_rot = 90 - deg          # convert CW-from-N to CCW-from-E
        if 90 < deg < 270:          # left half: flip for readability
            mpl_rot += 180
            ha = 'right'
        else:
            ha = 'left'

        gcf_text = f'GCF-{gcf_id}' if gcf_id is not None else '—'
        label = f'{gcf_text}  {genome.replace("_", " ")}'
        ax.text(θ, 1.035, label,
                ha=ha, va='center',
                fontsize=3.8, rotation=mpl_rot, rotation_mode='anchor',
                family='monospace', color='#333333')

    ax.set_ylim(0, 1.6)
    ax.set_rmax(1.6)
    ax.axis('off')

    # Legend
    present = set(coupling_classes.values())
    patches = [mpatches.Patch(color=COUPLING_COLORS[c], label=c)
               for c in COUPLING_ORDER if c in present]
    ax.legend(handles=patches, title='Coupling enzyme class',
              title_fontsize=8, fontsize=8, framealpha=0.92,
              loc='upper left', bbox_to_anchor=(0.0, 1.0),
              bbox_transform=ax.transAxes)

    gcf_count = len(set(m['gcf_id'] for m in id_to_meta.values()
                        if m['gcf_id'] is not None))
    fig.suptitle(
        f'Phosphonate BGC biosynthetic tree — all BGCs\n'
        f'NJ circular tree · BiG-SCAPE pairwise distances · {n} BGCs · {gcf_count} GCFs',
        fontsize=10, y=0.99, va='top')

    os.makedirs(outdir, exist_ok=True)
    out_png = os.path.join(outdir, 'all_bgcs_biosynthetic_tree_circular.png')
    out_svg = os.path.join(outdir, 'all_bgcs_biosynthetic_tree_circular.svg')
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    fig.savefig(out_svg,           bbox_inches='tight')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_svg}')
    plt.close(fig)


# ─── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='NJ tree of all BGCs from BiG-SCAPE distances',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--db',                  required=True)
    parser.add_argument('--coupling_annotation', required=True)
    parser.add_argument('--outdir',              required=True)
    parser.add_argument('--cutoff', type=float,  default=0.3)
    parser.add_argument('--layout', choices=['circular', 'linear'],
                        default='circular',
                        help='Tree layout (default: circular)')
    args = parser.parse_args()

    print('Loading coupling annotations...')
    coupling_classes = load_coupling_classes(args.coupling_annotation, region_only=True)
    print(f'  {len(coupling_classes)} BGC coupling labels loaded')

    print('Loading BGC data from database...')
    conn = sqlite3.connect(args.db)
    bgc_ids, id_to_meta, distances = load_bgc_data(conn, args.cutoff)
    conn.close()
    print(f'  {len(bgc_ids)} BGCs  |  {len(distances)} pairwise distances')

    tree = build_nj_tree(bgc_ids, distances)

    print('Plotting...')
    if args.layout == 'circular':
        plot_circular_tree(tree, bgc_ids, id_to_meta, coupling_classes, args.outdir)
    else:
        plot_tree(tree, bgc_ids, id_to_meta, coupling_classes, args.outdir)
    print('Done.')


if __name__ == '__main__':
    main()
