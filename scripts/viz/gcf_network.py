"""The GCF network, drawn as a distance matrix over family representatives.

BiG-SCAPE writes the complete all-pairs distance matrix and the report has never
shown it. A BGC-level view cannot help -- 1,476 regions is 1.09 million pairs, a
hairball at any threshold -- so each family collapses to its **representative**,
`family.center_id`, BiG-SCAPE's own medoid. 112 cells a side instead of 1,476.

This replaces the `isolation` column, which asked the same question as one number per
family and was read as a novelty rank however it was labelled (see
gcf_characterisation.py). A shape cannot be sorted.

**It is the same matrix the family-centre tree is built on** (`bgc_gcf_tree.py`,
Biosynthetic phylogeny pane), so the two panes are one object in two renderings: the
tree forces a hierarchy and shows the backbone, this shows every pair the hierarchy
had to discard. Before 2026-10-29 the pane was built on a different matrix -- single
linkage over all cross-family pairs -- and a reader who found two families adjacent
here and far apart in the tree had no way to reconcile them.

## Why a matrix and not a node-link diagram

Three node-link layouts were drawn against this data before the matrix was tried, and
the matrix wins on grounds that have nothing to do with taste:

  **It needs no threshold.** A node-link diagram must decide which edges to draw, so
  it needs an `EDGE_CUTOFF`, and no value of it is defensible. The matrix shows all
  6,216 pairs. The constant survives only to count neighbours and pick tooltips.

  **A hub reads as a stripe.** GCF-31 is within 0.50 of two thirds of all other
  families -- the grab-bag signature CLAUDE.md records for GCF-27 -- and as a full
  length line across every block that is unmissable, where 75 chords from one node
  are just a dense patch.

  **Sub-block structure survives.** Ordering by class then size puts near-diagonal
  blocks where sub-families are, and they are plainly visible. Neither the ring nor
  the arc showed them at all.

What it gives up: family size is not in the cells, so it rides a strip down the left
edge; and half the square is redundant by symmetry, which is accepted because scanning
a row OR a column is what makes a hub obvious.

The layouts it replaced, for the record. **Fruchterman-Reingold failed twice** -- plain
FR walks weakly-connected nodes into the clamp and draws them in rows along the frame,
and adding gravity made it worse, 89 of 112 pinned, because the repulsion of 111 other
nodes outweighs any centering term that does not crush the middle. The cause is the
hub, so *any* force layout resolves to a sunburst and the tuning only picks how ugly.
**A ring** fixed that and was shipped briefly, but it orders by class, so most edges
are within-class chords bundled unreadably against the rim while the interesting
cross-class ones have to be picked out of the crossings. **An arc diagram** -- the ring
unrolled -- is better on the long-range edges and fits a page better, but bundles the
same dense local block and still has to threshold.

## The ramp is clipped at 0.9, because that is where the data is

Representative distances on GenBank Enterobacterales: median **0.945**, 74% at 0.9 or
above, only 5% below 0.556. A linear ramp over [0, 1] spends almost all of its range
on pairs that are simply unrelated. One shade per tenth up to 0.9 and **blank at 0.9
and beyond** matches the decile histogram exactly, is one sentence to explain, and
keeps the drawn cells to about a quarter of the square -- which is also what keeps the
SVG small enough to inline.

A cell below the 0.30 family cutoff is legitimate, not a bug: BiG-SCAPE 2 clusters by
affinity propagation, NOT single linkage, so near-identical regions can land in
different families.

## Representatives are stable under the pipeline's pin

Worth stating because the opposite is written down elsewhere and is only true
unpinned. `docs/comparisons/bigscape_determinism/`: between identical runs **2 of 5
family centres changed unpinned, 0 with PYTHONHASHSEED=0**, which all four BiG-SCAPE
processes set. Every difference that survey saw was at distance >= 0.66, and **0 of
6,319 comparisons below 0.5 changed**. `bigscape_reproducibility/` adds that at 659
regions at most 2 members change family between runs, both boundary members of a
78-member family, which cannot move that family's medoid.
"""
import collections
import html as _html
import math
import os
import sqlite3

from utils.bgc_labels import label_from_path

# No longer decides what is drawn -- the matrix draws every pair. It is what "near"
# means when the caption counts a family's neighbours, and which cells are worth the
# bytes of an individual tooltip. Families are at least 0.30 apart by construction.
EDGE_CUTOFF = 0.50

# Shades run one per tenth up to BLANK_AT; at or beyond it a cell is left empty.
# See the module docstring: 74% of pairs sit at 0.9 or above on the run this was
# written against, so a ramp over the full [0, 1] would be almost entirely blank ink.
BLANK_AT = 0.9
BINS = 9

# One hue, per the report's sequential convention. **Index 0 is the DARKEST**,
# because index 0 is the closest bin: the ramp is ordered the way `shade()` indexes
# it, so that dark means similar. Writing it light-first and indexing by distance is
# the obvious mistake and silently draws the figure inside out -- the dense 0.8-0.9
# shoulder goes black and the close blocks fade to nothing.
RAMP = ['#17455a', '#2c6173', '#42788d', '#5b8fa4', '#78a5b7',
        '#97bac8', '#b4cdd7', '#cfdfe5', '#e7eef1']
_DIAGONAL = '#d2dcdf'      # a reading guide, not a measurement: a family to itself
_UNMEASURED = '#f3ece2'

# cell stays at 8 deliberately. A bigger cell makes the SVG wider than the report's
# content column (~1040 px at a 1440 window), and `width:100%` then scales the whole
# figure DOWN -- so a 9 px cell renders the 7 px row labels at 6.5 px and reads worse
# than an 8 px cell at 1:1. 8 is the largest cell that keeps every label legible.
CELL_MAX, CELL_MIN, BODY_MAX = 8, 3, 900
PAD_L, PAD_T = 92, 76          # row labels left; class names, strip, column labels above

_CLASS_FALLBACK = '#9aa7ad'
_SENTINELS = ('', '-', 'Unknown', 'unknown')


def family_network(db_path, cutoff=0.3):
    """(sizes, representatives, distances, measured_pairs, total_pairs) for a cutoff.

    sizes            {family_id: n_members}
    representatives  {family_id: (region label, organism)} for `family.center_id`
    distances        {(fam_a, fam_b): distance between the two representatives},
                     EVERY measured pair -- the matrix thresholds nothing

    The last two count how many of the n(n-1)/2 representative pairs the database
    actually holds a comparison for. On an unpartitioned run that is all of them; a
    merged partitioned database holds only within-partition comparisons, and a
    missing pair means NOT COMPARED, never distance zero. Those cells get their own
    fill and the caption says how many -- which is the trap the all-BGCs tree fell
    into, substituting a constant for 42.9% of its matrix with nothing printed.
    """
    db = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        rows = db.execute(
            """SELECT f.id, f.center_id, g.path, g.organism
                 FROM family f
                 JOIN bgc_record r ON r.id = f.center_id
                 LEFT JOIN gbk g ON g.id = r.gbk_id
                WHERE f.cutoff = ?""", (cutoff,)).fetchall()
        if len(rows) < 2:
            return {}, {}, {}, 0, 0

        fam_of_centre = {int(cid): fid for fid, cid, _, _ in rows}
        reps = {fid: (label_from_path(path or ''),
                      (organism or '').strip() or _organism_from_path(path))
                for fid, _, path, organism in rows}
        sizes = collections.Counter(
            f for (f,) in db.execute(
                """SELECT bf.family_id FROM bgc_record_family bf
                     JOIN family f ON f.id = bf.family_id
                    WHERE f.cutoff = ?""", (cutoff,)))

        ids = ','.join(str(c) for c in fam_of_centre)
        measured = {}
        for a, b, d in db.execute(
                'SELECT record_a_id, record_b_id, distance FROM distance '
                f'WHERE record_a_id IN ({ids}) AND record_b_id IN ({ids})'):
            fa, fb = fam_of_centre.get(a), fam_of_centre.get(b)
            if fa is None or fb is None or fa == fb:
                continue
            key = (fa, fb) if fa < fb else (fb, fa)
            # One pair can carry rows under several edge_params; keep the closest.
            if d < measured.get(key, 2.0):
                measured[key] = d
    finally:
        db.close()

    n = len(reps)
    return dict(sizes), reps, measured, len(measured), n * (n - 1) // 2


def _organism_from_path(path):
    """Fall back to the directory a region GBK sits in, which is the genome name."""
    return os.path.basename(os.path.dirname(path or '')).replace('_', ' ')


def _clean(cls):
    c = (cls or '').strip()
    return '' if c in _SENTINELS else c


def _matrix_order(nodes, size, classes, class_order):
    """Rows and columns: grouped by coupling class, largest first inside.

    The order is the whole reason the near-diagonal blocks mean anything. Grouping by
    class makes an off-diagonal block two classes whose representatives resemble each
    other, which is worth more attention than a block inside one class; sorting by
    size inside puts the families a reader has heard of at the top of their group.
    """
    known = [c for c in class_order if any(_clean(classes.get(f)) == c for f in nodes)]
    rest = sorted({_clean(classes.get(f)) for f in nodes} - set(known) - {''})
    buckets = known + rest + ['']                # unclassified last, always
    return [f for b in buckets
            for f in sorted((n for n in nodes if _clean(classes.get(n)) == b),
                            key=lambda n: (-size[n], str(n)))]


def _spans(order, classes):
    """[(class, first index, last index)] for the contiguous runs around the axis."""
    out = []
    for i, f in enumerate(order):
        c = _clean(classes.get(f))
        if out and out[-1][0] == c:
            out[-1][2] = i
        else:
            out.append([c, i, i])
    return out


def generate_gcf_network_html(db_path, cutoff=0.3, classes=None, colors=None,
                              class_order=(), hrefs=None, edge_cutoff=EDGE_CUTOFF):
    """Inline SVG of the family distance matrix, or '' when there is nothing to draw.

    classes      {family_id: coupling class} -- orders the axes and colours the
                 strips; anything missing or sentinel takes the neutral grey and
                 sits in the unclassified group, never a colour that reads as a class.
    colors       {class name: hex} from utils.constants.COUPLING_COLORS.
    class_order  the canonical class order, so the axes match the rest of the report
                 rather than whatever order the data arrived in.
    hrefs        {family_id: relative href} -- makes a row label a link to its page.
    edge_cutoff  what "near" means when counting a family's neighbours. It does NOT
                 decide what is drawn; every measured pair is.
    """
    size, reps, dist, measured, want = family_network(db_path, cutoff)
    if len(size) < 2:
        return ''
    classes, colors, hrefs = classes or {}, colors or {}, hrefs or {}

    order = _matrix_order(list(size), size, classes, list(class_order))
    n = len(order)
    at = {f: i for i, f in enumerate(order)}
    cell = max(CELL_MIN, min(CELL_MAX, BODY_MAX // n))
    body = cell * n
    W = PAD_L + body + 16
    H = PAD_T + body + 48
    x0 = lambda i: PAD_L + i * cell
    y0 = lambda i: PAD_T + i * cell

    def shade(d):
        """Bin index into RAMP, or None for a cell left blank."""
        if d >= BLANK_AT:
            return None
        return min(BINS - 1, int(d / BLANK_AT * BINS))

    parts = [f'<svg viewBox="0 0 {W} {H}" width="100%" '
             f'style="max-width:{W}px;height:auto;display:block;margin:0 auto;" '
             f'role="img" aria-label="Distance matrix of {n} gene cluster family '
             f'representatives, every family named on both axes and grouped by coupling '
             f'enzyme class. '
             f'Darker means more similar; a cell is blank at {BLANK_AT:g} or beyond.">',
             f'<rect x="{PAD_L}" y="{PAD_T}" width="{body}" height="{body}" '
             f'fill="#fbfcfc"/>']

    # Cells, collected per shade and emitted as ONE path each. 112 squared is 12,544
    # elements; as individual rects with their own fill that is over a megabyte of
    # SVG inlined into the report. Per-bin paths with integer coordinates bring the
    # whole figure under 100 kB, and the cells worth hovering get real elements below.
    buckets = collections.defaultdict(list)
    for (a, b), d in dist.items():
        k = shade(d)
        if k is None:
            continue
        for i, j in ((at[a], at[b]), (at[b], at[a])):
            buckets[k].append(f'M{x0(j)} {y0(i)}h{cell}v{cell}h-{cell}z')
    gaps = [f'M{x0(j)} {y0(i)}h{cell}v{cell}h-{cell}z'
            for ia, a in enumerate(order) for b in order[ia + 1:]
            if (min(a, b), max(a, b)) not in dist
            for i, j in ((at[a], at[b]), (at[b], at[a]))]
    if gaps:
        parts.append(f'<path d="{"".join(gaps)}" fill="{_UNMEASURED}"/>')
    for k in sorted(buckets):
        parts.append(f'<path d="{"".join(buckets[k])}" fill="{RAMP[k]}"/>')
    parts.append(f'<path d="{"".join(f"M{x0(i)} {y0(i)}h{cell}v{cell}h-{cell}z" for i in range(n))}" '
                 f'fill="{_DIAGONAL}"/>')

    # Individual elements only for the pairs a reader would chase, so the tooltips
    # exist where they are worth the bytes.
    near = sorted(((d, k) for k, d in dist.items() if d <= edge_cutoff))
    for d, (a, b) in near:
        tip = (f'GCF-{a} ({size[a]}) – GCF-{b} ({size[b]}): '
               f'representatives {d:.3f} apart')
        for i, j in ((at[a], at[b]), (at[b], at[a])):
            parts.append(f'<rect x="{x0(j)}" y="{y0(i)}" width="{cell}" '
                         f'height="{cell}" fill="{RAMP[shade(d)]}">'
                         f'<title>{_html.escape(tip)}</title></rect>')

    degree = collections.Counter(
        f for k, d in dist.items() if d <= edge_cutoff for f in k)

    # Class strips, left and top, plus a size bar down the left edge.
    for c, i0, i1 in _spans(order, classes):
        fill = colors.get(c, _CLASS_FALLBACK)
        run = (i1 - i0 + 1) * cell
        parts.append(f'<rect x="{PAD_L - 11}" y="{y0(i0)}" width="8" height="{run}" '
                     f'fill="{fill}"><title>{_html.escape(c or "no dominant class")}: '
                     f'{i1 - i0 + 1} families</title></rect>')
        parts.append(f'<rect x="{x0(i0)}" y="22" width="{run}" height="8" '
                     f'fill="{fill}"/>')
        # Above the matrix, not below: with every family labelled down the left edge
        # the eye enters at the top-left, and a class name under 896 px of cells is a
        # long way from the columns it names. Gated on measured width rather than on
        # family count -- a 4-family span is 32 px and cannot hold "Decarboxylase".
        name = c or 'unclassified'
        if i1 - i0 >= 1 and run >= 5.0 * len(name):
            parts.append(f'<text x="{x0(i0) + run / 2:.0f}" y="16" '
                         f'text-anchor="middle" font-size="10" '
                         f'font-family="system-ui,-apple-system,sans-serif" '
                         f'fill="#44555e">{_html.escape(name)}</text>')

    biggest = max(size.values())
    for f in order:
        w = max(1, round(20 * math.sqrt(size[f] / biggest)))
        parts.append(f'<rect x="{PAD_L - 35 + (20 - w)}" y="{y0(at[f])}" width="{w}" '
                     f'height="{cell}" fill="#8fa6b0"/>')

    # EVERY family is labelled on BOTH axes: down the left edge, and rotated above
    # the columns. 7 px text in an 8 px pitch is a dense table rather than a
    # collision, and rotating the column labels costs 29 px of height -- 3% of the
    # figure -- not the "roughly triple" an earlier version of this comment claimed.
    # That number came from reasoning about HORIZONTAL labels along the top, where
    # 28 px of text over an 8 px column really is impossible, and carrying the
    # conclusion across to rotated ones without measuring. Labelling both axes is
    # what lets a cell be read in from either direction, which is the whole point of
    # keeping the redundant lower triangle.
    #
    # The ten largest and anything near a quarter of the run are drawn darker, so the
    # families a reader came looking for still stand out from the other hundred.
    hub_at = max(4, n // 4)
    named = ({f for f in sorted(order, key=lambda g: -size[g])[:10]}
             | {f for f in order if degree.get(f, 0) >= hub_at})
    for f in order:
        label, organism = reps.get(f, ('', ''))
        deg = degree.get(f, 0)
        tip = (f'GCF-{f}: {size[f]} BGC{"s" if size[f] != 1 else ""}; '
               f'{deg} famil{"y" if deg == 1 else "ies"} within {edge_cutoff:g}')
        if label:
            tip += f'; representative {label}'
            if organism:
                tip += f' ({organism})'
        strong = f in named
        style = (f'dominant-baseline="middle" font-size="7" '
                 f'font-family="system-ui,-apple-system,sans-serif" '
                 f'fill="{"#1d3440" if strong else "#7b8a92"}"'
                 f'{" font-weight=\"600\"" if strong else ""}')
        name = f'GCF-{_html.escape(str(f))}'
        title = f'<title>{_html.escape(tip)}</title>'
        href = hrefs.get(f)

        ry = y0(at[f]) + cell / 2.0
        cx = x0(at[f]) + cell / 2.0
        cy = PAD_T - 5                      # column labels run UP from the matrix
        for text in (
                f'<text x="{PAD_L - 40}" y="{ry:.1f}" text-anchor="end" {style}>'
                f'{name}{title}</text>',
                f'<text x="{cx:.1f}" y="{cy}" transform="rotate(-90 {cx:.1f} {cy})" '
                f'text-anchor="start" {style}>{name}{title}</text>'):
            parts.append(f'<a href="{_html.escape(href)}">{text}</a>'
                         if href else text)

    # Ramp key along the bottom.
    key_y = PAD_T + body + 16
    for k in range(BINS):
        parts.append(f'<rect x="{PAD_L + k * 26}" y="{key_y}" width="26" height="9" '
                     f'fill="{RAMP[k]}"/>')
    parts.append(f'<rect x="{PAD_L + BINS * 26 + 6}" y="{key_y}" width="26" '
                 f'height="9" fill="#fbfcfc" stroke="#dfe4e6"/>')
    # Numbers alone are not a key. The first reader of this figure read the ramp
    # backwards -- reasonably, since "dark" reads as "more" and the quantity being
    # shown is DISsimilarity, so more distance is lighter. The ends say what they
    # mean in words.
    for lab, x, anchor in ((f'0.0 identical', PAD_L, 'start'),
                           (f'{BLANK_AT:g}', PAD_L + BINS * 26, 'end'),
                           (f'≥ {BLANK_AT:g} unrelated',
                            PAD_L + BINS * 26 + 38, 'start')):
        parts.append(f'<text x="{x}" y="{key_y + 20}" text-anchor="{anchor}" '
                     f'font-size="9" font-family="system-ui,-apple-system,sans-serif" '
                     f'fill="#66777e">{lab}</text>')
    parts.append(f'<text x="{PAD_L - 40}" y="{key_y + 7}" text-anchor="end" '
                 f'font-size="9" font-family="system-ui,-apple-system,sans-serif" '
                 f'fill="#66777e">distance</text>')
    parts.append(f'<text x="{PAD_L - 40}" y="{key_y + 20}" text-anchor="end" '
                 f'font-size="9" font-family="system-ui,-apple-system,sans-serif" '
                 f'fill="#8b979d">darker = closer</text>')
    parts.append('</svg>')

    lonely = sum(1 for f in order if degree.get(f, 0) == 0)
    hub = degree.most_common(1)[0] if degree else None
    used = [c for c in dict.fromkeys(_clean(classes.get(f)) for f in order) if c]
    legend = ''.join(
        f'<span style="display:inline-block;margin:0 14px 6px 0;white-space:nowrap;">'
        f'<span style="display:inline-block;width:11px;height:11px;border-radius:2px;'
        f'background:{colors.get(c, _CLASS_FALLBACK)};border:1px solid #17333f;'
        f'vertical-align:-1px;"></span>&nbsp;{_html.escape(c)}</span>'
        for c in used)

    hub_note = ''
    if hub and hub[1] >= hub_at:
        hub_note = (f' <strong>GCF-{_html.escape(str(hub[0]))}</strong> is the long '
                    f'stripe &mdash; it sits within {edge_cutoff:g} of {hub[1]} of '
                    f'them. A family near almost everything is usually a loose one '
                    f'rather than a central one, so check its cohesion before reading '
                    f'anything into that row.')

    gap_note = ''
    if want and measured < want:
        gap_note = (f'<p style="color:#8a5a2b;font-size:.88em;max-width:78ch;">'
                    f'<strong>{want - measured:,} of {want:,}</strong> representative '
                    f'pairs were never compared, and are drawn in a separate sand '
                    f'colour &mdash; this database holds only within-partition '
                    f'comparisons. Those cells mean <em>not measured</em>, not '
                    f'<em>far apart</em>.</p>')

    return f'''
    <h3>Gene Cluster Family Network</h3>
    <p style="color:#555;max-width:78ch;margin-bottom:4px;">
        Every family is one row and one column &mdash; its <strong>representative</strong>,
        the member BiG-SCAPE picks as the medoid, sitting closest to the rest of its
        family. A cell is the measured BiG-SCAPE distance between two
        representatives, <strong>darker the closer they are</strong>. Rows and columns
        are grouped by coupling enzyme class and sorted by family size inside each
        group; the bar down the left edge is that size. Hover a dark cell or a row
        label for the numbers and the representative&rsquo;s
        accession{', click a label for its family page' if hrefs else ''}.
    </p>
    <p style="color:#777;font-size:.88em;max-width:78ch;">
        <strong>Nothing is thresholded</strong> &mdash; all
        {len(dist):,} measured pairs are here, which is why this is a matrix and not a
        network of nodes and edges. Shades run one per tenth of distance and a cell is
        left blank at {BLANK_AT:g} or beyond, where the two families simply are not
        related; on a typical run that is most of the square. This is the same
        distance matrix as the family-centre tree under Biosynthetic phylogeny, so the
        two panes agree by construction &mdash; the tree shows the backbone, this
        shows every pair the tree had to discard. Families are at least 0.30 apart by
        definition, so a cell darker than that means BiG-SCAPE&rsquo;s affinity
        propagation split near-identical regions, which it legitimately does.
        <strong>{lonely}</strong> of <strong>{n}</strong> families have no other
        family within {edge_cutoff:g}.{hub_note}
    </p>
    <p style="color:#777;font-size:.88em;max-width:78ch;">
        <strong>A dark off-diagonal cell is worth taking seriously.</strong> Some
        families here carry an <em>identical</em> rule core and separate only on
        flanking chromosomal context, which differs by host species &mdash; so they
        are one chemistry reported as several. Measured on GenBank Enterobacterales,
        the two largest <em>Pantoea</em> pantaphos families sit 0.680 apart over the
        whole region and <strong>0.237</strong> apart over the biosynthetic cassette
        alone. Families are affinity propagation&rsquo;s, and it splits inside a
        connected component whatever the cutoff; this matrix shows the distances so
        that judgement stays with the reader. See
        <code>docs/comparisons/bigscape_cutoff_sweep/</code>.
    </p>
    {gap_note}
    <div style="margin:10px 0 14px;font-size:.88em;color:#444;">{legend}</div>
    <div style="border:1px solid #e3e6e4;border-radius:10px;padding:8px;background:#fff;">
        {''.join(parts)}
    </div>
    '''
