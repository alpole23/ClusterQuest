"""Cladogram layout and drawing helpers for Bio.Phylo trees.

Layout is done by annotating clades in place:
  linear   → ._depth (x, in node depth) and ._x   (leaf order position)
  circular → ._depth (radius, in node depth) and ._angle (radians)

Both bgc_all_bgcs_tree.py and bgc_gcf_tree.py draw cladograms this way.
"""

import numpy as np


# ─── Linear layout ────────────────────────────────────────────────────────────

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


# ─── Circular (fan) layout ────────────────────────────────────────────────────

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
