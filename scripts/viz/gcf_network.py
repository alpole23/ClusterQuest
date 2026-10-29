"""The GCF network: one node per family representative, from BiG-SCAPE's distances.

BiG-SCAPE writes the complete all-pairs distance matrix, and the report has never
shown it as a graph. A BGC-level network cannot help: 1,476 regions is 1.09 million
pairs, and the picture is a hairball whatever you threshold it at. **Collapsing each
family to its representative is what makes the relationships visible** -- 112 nodes
instead, where the question "is this family near anything else?" is answered by
looking at it.

This replaces the `isolation` column, which asked the same question as one number per
family and was read as a novelty rank however it was labelled (see
gcf_characterisation.py). A shape cannot be sorted.

## A node is a real BGC, and an edge is one measured comparison

The node is `family.center_id` -- the medoid, the member with the smallest summed
distance to the rest of its family, chosen by BiG-SCAPE itself. An edge is the
**stored distance between two representatives**: a single comparison between two real
regions, which the tooltip names and the family page links to.

**This is the same matrix the family-centre tree is built on** (`bgc_gcf_tree.py`,
Biosynthetic phylogeny pane). That consistency is the point of the change. The two
panes are now one object in two layouts -- the tree forces a hierarchy and shows the
backbone, the network keeps the cycles and the hub that a tree has to discard -- so a
reader who finds GCF-5 beside GCF-18 in one can find the same relationship in the
other. Before this they could not: the network was built on a different matrix.

## What the previous version did, and what switching cost

Until 2026-10-29 an edge was **single linkage**: the smallest distance between any
member of A and any member of B. Measured on GenBank Enterobacterales, switching to
representatives at the 0.50 cutoff:

    edges                327 -> 258        isolated families   9 -> 13
    largest degree        79 -> 75         representative-only edges: 0

Zero in that last cell at every threshold tried (0.30 through 0.70), and necessarily
so: the representative pair is one of the cross-family pairs single linkage minimises
over, so distance(rep_A, rep_B) >= min(A, B) always. **The representative network is a
strict subset of the closest-pair network.** It is the conservative reading.

The difference is small where it is measurable at all. 43 of 112 families are
singletons, where the representative *is* the only member and the two definitions
coincide by construction; over the 69 multi-member families the median
representative-minus-closest gap is +0.021, the mean +0.049, the worst +0.717.

So the trade is: an edge now means "these two families resemble each other" rather
than "these two families each hold one region that resembles the other's". The second
is a lead worth chasing and the first is a statement you can publish, and the caption
no longer has to talk a reader out of the obvious reading of a line between two nodes.
What is lost -- a family pair touching closely at one point while their bulk differs
-- is cohesion, which `docs/comparisons/family_cohesion/` measures directly and
better; that is where GCF-27's 0.10 across 21 sub-clusters is recorded.

## Representatives are stable under the pipeline's pin

Worth stating because the opposite is written down elsewhere and is only true
unpinned. `docs/comparisons/bigscape_determinism/`: between identical runs **2 of 5
family centres changed unpinned, 0 with PYTHONHASHSEED=0**, which the pipeline sets.
Every difference that survey saw was at distance >= 0.66, and **0 of 6,319 comparisons
below 0.5 changed** -- so the edges this figure draws sit in the regime that was
reproducible even without the pin. `bigscape_reproducibility/` adds that at 659
regions at most 2 members change family between runs, both boundary members of a
78-member family, which cannot move that family's medoid.

## Why a ring and not a force layout

The first version was Fruchterman-Reingold, and it failed on this data twice. Plain
FR walks weakly-connected nodes to the wall, where the clamp parks them in straight
rows; adding gravity to pull them back made it worse, 89 of 112 pinned, because the
repulsion of 111 other nodes outweighs any centering term you can add without
crushing the middle. The deeper problem is the data: GenBank Enterobacterales has one
family of degree 75 -- GCF-31 is within 0.50 of two thirds of all other families, the
grab-bag signature CLAUDE.md records for GCF-27 -- so ANY force layout resolves to a
sunburst around it, and the tuning only decides how ugly the sunburst is.

A ring sidesteps all of it. Position carries one honest meaning (coupling class, then
size), every edge is visible as a chord, a hub reads immediately as a node with
chords to everywhere, and the layout is deterministic by construction rather than by
seeding -- which the report's byte-reproducibility requirement needs, and which a
seeded RNG does not actually give you, since family ids shift between runs.
"""
import collections
import html as _html
import math
import os
import sqlite3

from utils.bgc_labels import label_from_path

# Families are at least 0.30 apart by construction, so an edge can only live above
# that. 0.50 keeps the pairs a reader would actually chase -- 258 edges over 112
# families on GenBank Enterobacterales, about 2 per family -- where 0.70 gives 502
# and the picture starts to fill in. Raising it does not reveal more structure; it
# reveals that almost everything is weakly related to almost everything.
EDGE_CUTOFF = 0.50

W, H = 1080, 780
CX, CY = W / 2.0, H / 2.0 - 8
RING = 286.0
# Chords bow toward the middle. 0.30 is enough to separate a short chord from the
# straight line between its ends without sending every long chord through one point.
BOW = 0.30

_CLASS_FALLBACK = '#9aa7ad'
_SENTINELS = ('', '-', 'Unknown', 'unknown')


def family_network(db_path, cutoff=0.3, edge_cutoff=EDGE_CUTOFF):
    """(sizes, representatives, edges, measured_pairs, total_pairs) for one cutoff.

    sizes            {family_id: n_members}
    representatives  {family_id: (region label, organism)} for `family.center_id`
    edges            {(fam_a, fam_b): distance between the two representatives},
                     thresholded at `edge_cutoff`

    The last two count how many of the n(n-1)/2 representative pairs the database
    actually holds a comparison for. On an unpartitioned run that is all of them; a
    merged partitioned database holds only within-partition comparisons, and a
    missing pair means NOT COMPARED, never distance zero. The caller says so rather
    than quietly drawing a sparser graph -- which is the exact trap the all-BGCs tree
    fell into, substituting a constant for 42.9% of its matrix with nothing printed.
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
    return (dict(sizes), reps,
            {k: v for k, v in measured.items() if v <= edge_cutoff},
            len(measured), n * (n - 1) // 2)


def _organism_from_path(path):
    """Fall back to the directory a region GBK sits in, which is the genome name."""
    return os.path.basename(os.path.dirname(path or '')).replace('_', ' ')


def _clean(cls):
    c = (cls or '').strip()
    return '' if c in _SENTINELS else c


def _ring_order(nodes, size, classes, class_order):
    """Families around the ring: grouped by coupling class, largest first inside.

    Grouping by class is what makes the picture answer a second question for free --
    a chord between two groups is two different coupling classes whose
    representatives are close, which is worth more attention than a chord inside one.
    """
    known = [c for c in class_order if any(_clean(classes.get(f)) == c for f in nodes)]
    rest = sorted({_clean(classes.get(f)) for f in nodes} - set(known) - {''})
    buckets = known + rest + ['']                # unclassified last, always
    return [f for b in buckets
            for f in sorted((n for n in nodes if _clean(classes.get(n)) == b),
                            key=lambda n: (-size[n], str(n)))]


def generate_gcf_network_html(db_path, cutoff=0.3, classes=None, colors=None,
                              class_order=(), hrefs=None, edge_cutoff=EDGE_CUTOFF):
    """Inline SVG of the family network, or '' when there is nothing to draw.

    classes      {family_id: coupling class} -- colours the nodes and groups the
                 ring; anything missing or sentinel takes the neutral grey and sits
                 in the unclassified group, never a colour that reads as a class.
    colors       {class name: hex} from utils.constants.COUPLING_COLORS.
    class_order  the canonical class order, so the ring matches the rest of the
                 report rather than whatever order the data arrived in.
    hrefs        {family_id: relative href} -- makes a node a link to its page.
    """
    size, reps, edges, measured, want = family_network(db_path, cutoff, edge_cutoff)
    if len(size) < 2:
        return ''
    classes, colors, hrefs = classes or {}, colors or {}, hrefs or {}

    order = _ring_order(list(size), size, classes, list(class_order))
    n = len(order)
    # -90deg start so the ring opens at the top; the first class group begins there.
    ang = {f: (2 * math.pi * i / n) - math.pi / 2 for i, f in enumerate(order)}
    pos = {f: (CX + RING * math.cos(a), CY + RING * math.sin(a)) for f, a in ang.items()}

    biggest = max(size.values())
    def radius(f):
        return 4.5 + 17.0 * math.sqrt(size[f] / biggest)

    degree = collections.Counter(f for k in edges for f in k)
    pull = {k: (edge_cutoff - d) / edge_cutoff for k, d in edges.items()}

    parts = [f'<svg viewBox="0 0 {W} {H}" width="100%" '
             f'style="max-width:{W}px;height:auto;display:block;margin:0 auto;" '
             f'role="img" aria-label="Ring of {n} gene cluster families grouped by '
             f'coupling enzyme class. Each node is the family representative; a chord '
             f'joins two representatives within {edge_cutoff:g} BiG-SCAPE distance.">']

    # Class arcs just outside the ring, so the grouping is visible without reading
    # the node colours one by one.
    groups = []
    for i, f in enumerate(order):
        c = _clean(classes.get(f))
        if groups and groups[-1][0] == c:
            groups[-1][2] = i
        else:
            groups.append([c, i, i])
    for c, i0, i1 in groups:
        if i1 - i0 < 1:
            continue
        a0, a1 = ang[order[i0]], ang[order[i1]]
        rr = RING + 26
        x0, y0 = CX + rr * math.cos(a0), CY + rr * math.sin(a0)
        x1, y1 = CX + rr * math.cos(a1), CY + rr * math.sin(a1)
        large = 1 if (a1 - a0) > math.pi else 0
        parts.append(f'<path d="M {x0:.1f} {y0:.1f} A {rr} {rr} 0 {large} 1 '
                     f'{x1:.1f} {y1:.1f}" fill="none" stroke-linecap="round" '
                     f'stroke="{colors.get(c, _CLASS_FALLBACK)}" stroke-width="5" '
                     f'stroke-opacity="0.75"><title>'
                     f'{_html.escape(c or "no dominant class")}: '
                     f'{i1 - i0 + 1} families</title></path>')

    # Chords, weakest first so the close pairs are not buried under them.
    for (a, b), d in sorted(edges.items(), key=lambda kv: -kv[1]):
        (x1, y1), (x2, y2) = pos[a], pos[b]
        mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
        qx, qy = CX + (mx - CX) * BOW, CY + (my - CY) * BOW
        w = pull[(a, b)]
        parts.append(f'<path d="M {x1:.1f} {y1:.1f} Q {qx:.1f} {qy:.1f} '
                     f'{x2:.1f} {y2:.1f}" fill="none" stroke="#43606d" '
                     f'stroke-opacity="{0.10 + 0.55 * w:.2f}" '
                     f'stroke-width="{0.6 + 2.2 * w:.2f}"><title>'
                     f'GCF-{_html.escape(str(a))} &#8211; GCF-{_html.escape(str(b))}: '
                     f'representatives {d:.3f} apart</title></path>')

    for f in order:
        x, y = pos[f]
        r = radius(f)
        cls = _clean(classes.get(f))
        deg = degree.get(f, 0)
        label, organism = reps.get(f, ('', ''))
        tip = (f'GCF-{f}: {size[f]} BGC{"s" if size[f] != 1 else ""}'
               f'{", " + cls if cls else ", no dominant class"}; '
               f'{deg} famil{"y" if deg == 1 else "ies"} within {edge_cutoff:g}')
        if label:
            tip += f'\nrepresentative {label}'
            if organism:
                tip += f' ({organism})'
        circle = (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                  f'fill="{colors.get(cls, _CLASS_FALLBACK)}" fill-opacity="0.92" '
                  f'stroke="#17333f" stroke-width="1.1">'
                  f'<title>{_html.escape(tip)}</title></circle>')
        href = hrefs.get(f)
        parts.append(f'<a href="{_html.escape(href)}">{circle}</a>' if href else circle)

    # Labels outside the ring, for the families a reader would look for: the ten
    # largest and anything joined to a quarter of the run. Labelling all 112 makes a
    # solid band of text with the network behind it.
    #
    # Each label runs RADIALLY, along its own spoke, rather than horizontally. The
    # ring is sorted by size inside each class group, so the ten largest are always
    # ADJACENT -- and horizontal labels on adjacent nodes near the top or bottom of
    # the ring differ by ~16 px of x against a ~44 px word, which overprinted them
    # into one blob. Radially they are separated by the angular gap instead: 11 px of
    # text height against 18 px of arc at 112 nodes, so they cannot collide at any
    # family count this figure is readable at.
    hub_at = max(4, n // 4)
    named = ({f for f in sorted(order, key=lambda g: -size[g])[:10]}
             | {f for f in order if degree.get(f, 0) >= hub_at})
    for f in sorted(named, key=lambda g: order.index(g)):
        a = ang[f]
        rr = RING + 34                       # outside the class arcs at RING + 26
        x, y = CX + rr * math.cos(a), CY + rr * math.sin(a)
        right = math.cos(a) >= 0
        deg = math.degrees(a) + (0 if right else 180)
        parts.append(f'<text x="{x:.1f}" y="{y:.1f}" '
                     f'transform="rotate({deg:.1f} {x:.1f} {y:.1f})" '
                     f'text-anchor="{"start" if right else "end"}" '
                     f'dominant-baseline="middle" font-size="11" '
                     f'font-family="system-ui,-apple-system,sans-serif" fill="#2b3f49">'
                     f'GCF-{_html.escape(str(f))}</text>')
    parts.append('</svg>')

    singles = sum(1 for f in order if degree.get(f, 0) == 0)
    hub = degree.most_common(1)[0] if degree else None
    used = [c for c in dict.fromkeys(_clean(classes.get(f)) for f in order) if c]
    legend = ''.join(
        f'<span style="display:inline-block;margin:0 14px 6px 0;white-space:nowrap;">'
        f'<span style="display:inline-block;width:11px;height:11px;border-radius:50%;'
        f'background:{colors.get(c, _CLASS_FALLBACK)};border:1px solid #17333f;'
        f'vertical-align:-1px;"></span>&nbsp;{_html.escape(c)}</span>'
        for c in used)

    hub_note = ''
    if hub and hub[1] >= hub_at:
        hub_note = (f' <strong>GCF-{_html.escape(str(hub[0]))}</strong> is joined to '
                    f'{hub[1]} of them &mdash; a family near almost everything is '
                    f'usually a loose one rather than a central one, so check its '
                    f'cohesion before reading anything into those chords.')

    # A merged partitioned database holds only within-partition comparisons, so some
    # representative pairs were never measured. Say so; do not let a sparser picture
    # pass for a more isolated one.
    gap_note = ''
    if want and measured < want:
        gap_note = (f'<p style="color:#8a5a2b;font-size:.88em;max-width:78ch;">'
                    f'<strong>{want - measured:,} of {want:,}</strong> representative '
                    f'pairs were never compared &mdash; this database holds only '
                    f'within-partition comparisons. A missing chord there means '
                    f"<em>not measured</em>, not <em>far apart</em>.</p>")

    return f'''
    <h3>Gene Cluster Family Network</h3>
    <p style="color:#555;max-width:78ch;margin-bottom:4px;">
        Every family is one node &mdash; its <strong>representative</strong>, the
        member BiG-SCAPE picks as the medoid, sitting closest to the rest of its
        family. Nodes are sized by how many BGCs the family holds, coloured by
        coupling enzyme class and grouped around the ring by that class. A chord is
        the measured BiG-SCAPE distance between two representatives, drawn when it is
        within <strong>{edge_cutoff:g}</strong>; the darker and thicker the chord, the
        closer the pair. Hover any node or chord for the numbers and the
        representative&rsquo;s accession{', click a node for its family page' if hrefs else ''}.
    </p>
    <p style="color:#777;font-size:.88em;max-width:78ch;">
        This is the same distance matrix as the family-centre tree under Biosynthetic
        phylogeny, so the two panes agree by construction &mdash; the tree shows the
        backbone, this shows the cycles a tree has to discard. A chord is a claim
        about the two families, not about their closest members: two families can hold
        a closer pair of individual regions than their representatives show, and this
        view deliberately does not draw that. Families are at least 0.30 apart by
        definition, so a chord shorter than that means BiG-SCAPE&rsquo;s affinity
        propagation split near-identical regions, which it legitimately does.
        <strong>{singles}</strong> of <strong>{n}</strong> families have no neighbour
        within {edge_cutoff:g}.{hub_note}
    </p>
    {gap_note}
    <div style="margin:10px 0 14px;font-size:.88em;color:#444;">{legend}</div>
    <div style="border:1px solid #e3e6e4;border-radius:10px;padding:8px;background:#fff;">
        {''.join(parts)}
    </div>
    '''
