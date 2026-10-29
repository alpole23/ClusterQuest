"""The GCF network: one node per family, drawn from BiG-SCAPE's own distances.

BiG-SCAPE writes the complete all-pairs distance matrix, and the report has never
shown it. A BGC-level network cannot help: 1,476 regions is 1.09 million pairs, and
the picture is a hairball whatever you threshold it at. **Collapsing each family to a
single node is what makes the relationships visible** -- 112 nodes instead, where the
question "is this family near anything else?" is answered by looking at it.

This replaces the `isolation` column, which asked the same question as one number per
family and was read as a novelty rank however it was labelled (see
gcf_characterisation.py). A shape cannot be sorted.

**An edge is single-linkage: the SMALLEST distance between any member of A and any
member of B.** Not the median, and deliberately not the centre-to-centre distance.
The question a reader brings is "could these two families be the same chemistry?",
and one close pair is enough to make that worth checking -- a median would bury it
under the bulk of two families that are mostly unalike. It also means an edge is a
claim about a pair of REGIONS, not about the families as wholes, which the caption
says.

Note BiG-SCAPE 2 clusters by affinity propagation inside each connected component,
NOT by single linkage, so two families can legitimately hold near-identical members:
on GenBank Enterobacterales the closest cross-family pair is 0.000. Those are the
edges most worth looking at and they would be impossible under single-linkage
clustering, so do not read one as a bug.

## Why a ring and not a force layout

The first version was Fruchterman-Reingold, and it failed on this data twice. Plain
FR walks weakly-connected nodes to the wall, where the clamp parks them in straight
rows; adding gravity to pull them back made it worse, 89 of 112 pinned, because the
repulsion of 111 other nodes outweighs any centering term you can add without
crushing the middle. The deeper problem is the data: GenBank Enterobacterales has one
family of degree 79 -- GCF-31 is within 0.50 of 70% of all other families, the
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
import sqlite3

# Families are at least 0.30 apart by construction, so an edge can only live above
# that. 0.50 keeps the pairs a reader would actually chase -- 327 edges over 112
# families on GenBank Enterobacterales, about 3 per family -- where 0.70 gives 691
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


def family_edges(db_path, cutoff=0.3, edge_cutoff=EDGE_CUTOFF):
    """({family_id: n_members}, {(fam_a, fam_b): min_distance}) for one cutoff.

    One pass over the distance table, keeping only the running minimum per family
    pair, so memory is O(pairs of families) rather than O(pairs of regions) -- 6,216
    against 1,088,550 on the run this was written against.
    """
    db = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        fam_of = dict(db.execute(
            """SELECT bf.record_id, bf.family_id
                 FROM bgc_record_family bf
                 JOIN family f ON f.id = bf.family_id
                WHERE f.cutoff = ?""", (cutoff,)))
        if not fam_of:
            return {}, {}
        size = collections.Counter(fam_of.values())
        best = {}
        for a, b, d in db.execute('SELECT record_a_id, record_b_id, distance FROM distance'):
            fa, fb = fam_of.get(a), fam_of.get(b)
            if fa is None or fb is None or fa == fb:
                continue
            key = (fa, fb) if fa < fb else (fb, fa)
            if d < best.get(key, 2.0):
                best[key] = d
    finally:
        db.close()
    return dict(size), {k: v for k, v in best.items() if v <= edge_cutoff}


def _clean(cls):
    c = (cls or '').strip()
    return '' if c in _SENTINELS else c


def _ring_order(nodes, size, classes, class_order):
    """Families around the ring: grouped by coupling class, largest first inside.

    Grouping by class is what makes the picture answer a second question for free --
    a chord between two groups is two different coupling classes holding a pair of
    near-identical regions, which is worth more attention than a chord inside one.
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
    size, edges = family_edges(db_path, cutoff, edge_cutoff)
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
             f'coupling enzyme class. A chord joins two families that hold at least '
             f'one pair of regions within {edge_cutoff:g} BiG-SCAPE distance.">']

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
                     f'closest pair {d:.3f}</title></path>')

    for f in order:
        x, y = pos[f]
        r = radius(f)
        cls = _clean(classes.get(f))
        deg = degree.get(f, 0)
        tip = (f'GCF-{f}: {size[f]} BGC{"s" if size[f] != 1 else ""}'
               f'{", " + cls if cls else ", no dominant class"}; '
               f'{deg} famil{"y" if deg == 1 else "ies"} within {edge_cutoff:g}')
        circle = (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                  f'fill="{colors.get(cls, _CLASS_FALLBACK)}" fill-opacity="0.92" '
                  f'stroke="#17333f" stroke-width="1.1">'
                  f'<title>{_html.escape(tip)}</title></circle>')
        href = hrefs.get(f)
        parts.append(f'<a href="{_html.escape(href)}">{circle}</a>' if href else circle)

    # Labels outside the ring, for the families a reader would look for: the ten
    # largest and anything joined to a quarter of the run. Labelling all 112 makes a
    # solid band of text with the network behind it.
    hub_at = max(4, n // 4)
    named = ({f for f in sorted(order, key=lambda g: -size[g])[:10]}
             | {f for f in order if degree.get(f, 0) >= hub_at})
    for f in sorted(named, key=lambda g: order.index(g)):
        a = ang[f]
        rr = RING + 40
        x, y = CX + rr * math.cos(a), CY + rr * math.sin(a)
        right = math.cos(a) >= 0
        parts.append(f'<text x="{x:.1f}" y="{y:.1f}" '
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

    return f'''
    <h3>Gene Cluster Family Network</h3>
    <p style="color:#555;max-width:78ch;margin-bottom:4px;">
        Every family is one node, sized by how many BGCs it holds, coloured by its
        coupling enzyme class and grouped around the ring by that class. A chord joins
        two families when at least one region in each sits within
        <strong>{edge_cutoff:g}</strong> BiG-SCAPE distance of the other; the darker
        and thicker the chord, the closer that pair. Hover any node or chord for the
        numbers{', click a node for its family page' if hrefs else ''}.
    </p>
    <p style="color:#777;font-size:.88em;max-width:78ch;">
        A chord is the <em>closest pair</em>, not an average, so it says two families
        have something in common &mdash; not that they are alike throughout. Families
        are at least 0.30 apart by definition, so a chord shorter than that means
        BiG-SCAPE&rsquo;s affinity propagation split near-identical regions, which it
        legitimately does. <strong>{singles}</strong> of <strong>{n}</strong> families
        have no neighbour within {edge_cutoff:g}.{hub_note}
    </p>
    <div style="margin:10px 0 14px;font-size:.88em;color:#444;">{legend}</div>
    <div style="border:1px solid #e3e6e4;border-radius:10px;padding:8px;background:#fff;">
        {''.join(parts)}
    </div>
    '''
