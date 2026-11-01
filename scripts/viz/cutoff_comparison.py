"""How the clustering changes across BiG-SCAPE cutoffs, in one section.

`bigscape_cutoffs` accepts a comma-separated list and BiG-SCAPE writes every one of
them into the same database, keyed by `family.cutoff`. Until now the report picked the
first and ignored the rest, so the only way to compare cutoffs was to run the pipeline
once per cutoff and diff two reports by eye.

The question a reader has is not "what does 0.40 look like" -- that is just another
clustering -- but **"does the cutoff I picked change the answer?"**. Two things decide
that, and the table reports both:

  shape      families, largest, singletons, and how many regions land in a family of
             more than one. A cutoff that moves the family count while leaving the
             largest family alone is splitting the tail, not restructuring.
  agreement  between CONSECUTIVE cutoffs, how much co-membership survives. Reported as
             pairs SPLIT and pairs MERGED separately, not just as an index, because
             they are different failures: a merge joins chemistry that a tighter cutoff
             kept apart, a split breaks a family you may already have written about.
             The Adjusted Rand Index is given alongside, since it is the number people
             quote, but it hides which direction the change went.

Adjusted Rand is computed from the contingency table rather than by enumerating pairs:
1,476 regions is 1.09 million pairs and the closed form is O(n + k1*k2).

**This section is not cutoff-aware and must not be.** Every other cutoff-dependent
block in the report shows one cutoff at a time and is switched by the selector at the
top; this one exists to show them together, and switching it would defeat it.
"""
import collections
import html as _html
import math
import sqlite3


def _comb2(n):
    return n * (n - 1) // 2


def memberships(db_path):
    """{cutoff: {record_id: family_id}}, region records only, ordered by cutoff.

    Region-only matters: BiG-SCAPE assigns families to region records, and on the
    runs seen so far nothing else appears in `bgc_record_family` -- but the filter is
    the same one every other query in this report applies, and a merged partitioned
    database is exactly where that assumption would stop holding quietly.
    """
    db = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        out = collections.defaultdict(dict)
        for cut, rec, fam in db.execute(
                """SELECT f.cutoff, bf.record_id, bf.family_id
                     FROM bgc_record_family bf
                     JOIN family f ON f.id = bf.family_id
                     JOIN bgc_record r ON r.id = bf.record_id
                    WHERE r.record_type = 'region'"""):
            out[cut][rec] = fam
    finally:
        db.close()
    return {c: out[c] for c in sorted(out)}


def shape(member):
    """Family-count statistics for one cutoff's {record: family} map."""
    sizes = collections.Counter(member.values())
    vals = sorted(sizes.values())
    return {
        'regions': len(member),
        'families': len(sizes),
        'largest': vals[-1] if vals else 0,
        'singletons': sum(1 for v in vals if v == 1),
        'median': vals[len(vals) // 2] if vals else 0,
        'in_multi': sum(v for v in vals if v > 1),
    }


def agreement(a, b):
    """Co-membership change from clustering `a` to clustering `b`.

    Over the records both clusterings assign. `split` is pairs together in A and apart
    in B; `merged` is the reverse. ARI is the standard contingency-table form.
    """
    shared = sorted(set(a) & set(b))
    if len(shared) < 2:
        return None
    ct = collections.Counter((a[r], b[r]) for r in shared)
    ra = collections.Counter(a[r] for r in shared)
    rb = collections.Counter(b[r] for r in shared)
    n = len(shared)
    both = sum(_comb2(v) for v in ct.values())
    sa = sum(_comb2(v) for v in ra.values())
    sb = sum(_comb2(v) for v in rb.values())
    total = _comb2(n)
    expected = sa * sb / total if total else 0.0
    denom = 0.5 * (sa + sb) - expected
    return {
        'n': n,
        'split': sa - both,          # together in A, apart in B
        'merged': sb - both,         # apart in A, together in B
        'kept': both,
        'ari': (both - expected) / denom if denom else 1.0,
    }


def build_cutoff_comparison_section(db_path, primary=None):
    """The section HTML, or '' when the database holds fewer than two cutoffs."""
    try:
        members = memberships(db_path)
    except Exception as exc:                      # a missing or unreadable db is not
        print(f'Warning: cutoff comparison unavailable: {exc}')   # worth failing on
        return ''
    if len(members) < 2:
        return ''

    cuts = list(members)
    shapes = {c: shape(members[c]) for c in cuts}

    def _row(c):
        # Same PEP 701 constraint as viz/gcf_network.py: no backslash inside an
        # f-string expression below Python 3.12. Build the fragments first.
        is_primary = primary is not None and abs(c - primary) < 1e-9
        tr = ' style="background:#eef4f8;"' if is_primary else ''
        tag = (' <span style="color:#2c5aa0;font-size:.82em;">primary</span>'
               if is_primary else '')
        num = 'padding:7px 12px;text-align:right;font-variant-numeric:tabular-nums;'
        s = shapes[c]
        return (
            f'<tr{tr}>'
            f'<td style="padding:7px 12px;font-variant-numeric:tabular-nums;">'
            f'<strong>{c:g}</strong>{tag}</td>'
            f'<td style="{num}">{s["families"]}</td>'
            f'<td style="{num}">{s["largest"]}</td>'
            f'<td style="{num}">{s["singletons"]}</td>'
            f'<td style="{num}">{s["median"]}</td>'
            f'<td style="{num}">{s["in_multi"]} '
            f'<span style="color:#888;font-size:.85em;">'
            f'({100.0 * s["in_multi"] / s["regions"]:.0f}%)</span></td>'
            f'</tr>')

    rows = ''.join(_row(c) for c in cuts)

    pair_rows = ''
    for lo, hi in zip(cuts, cuts[1:]):
        ag = agreement(members[lo], members[hi])
        if not ag:
            continue
        tone = '#2f7d5d' if ag['ari'] >= 0.95 else ('#b07d2b' if ag['ari'] >= 0.8 else '#a8452c')
        pair_rows += (
            f'<tr><td style="padding:7px 12px;font-variant-numeric:tabular-nums;">'
            f'{lo:g} &rarr; {hi:g}</td>'
            f'<td style="padding:7px 12px;text-align:right;font-variant-numeric:tabular-nums;'
            f'color:{tone};font-weight:600;">{ag["ari"]:.4f}</td>'
            f'<td style="padding:7px 12px;text-align:right;font-variant-numeric:tabular-nums;">'
            f'{ag["split"]:,}</td>'
            f'<td style="padding:7px 12px;text-align:right;font-variant-numeric:tabular-nums;">'
            f'{ag["merged"]:,}</td>'
            f'<td style="padding:7px 12px;text-align:right;font-variant-numeric:tabular-nums;">'
            f'{ag["kept"]:,}</td></tr>')

    th = ('padding:7px 12px;text-align:right;font-weight:600;color:#fff;'
          'background:#2c3e50;')
    th_l = th.replace('text-align:right', 'text-align:left')

    return f'''
    <h3>Cutoff Comparison</h3>
    <p style="color:#555;max-width:78ch;">
        This run clustered at <strong>{len(cuts)}</strong> cutoffs, all written to the
        same BiG-SCAPE database. Every other cutoff-dependent section of this report
        shows one at a time, chosen with the selector at the top of the page; this one
        shows them together, because the question it answers is whether the choice
        matters.
    </p>
    <table class="data-table" style="margin:14px 0 26px;border-collapse:collapse;">
        <thead><tr>
            <th style="{th_l}">Cutoff</th>
            <th style="{th}">Families</th>
            <th style="{th}">Largest</th>
            <th style="{th}">Singletons</th>
            <th style="{th}">Median size</th>
            <th style="{th}">Regions in a family &gt; 1</th>
        </tr></thead>
        <tbody>{rows}</tbody>
    </table>
    <h3 style="margin-top:4px;">Agreement between consecutive cutoffs</h3>
    <p style="color:#555;max-width:78ch;">
        Over the regions both cutoffs assign. <strong>Split</strong> counts pairs of
        regions in the same family at the looser cutoff and in different families at
        the tighter one; <strong>merged</strong> counts the reverse. They are reported
        separately because they are different problems &mdash; a merge joins chemistry
        the tighter cutoff kept apart, a split breaks a family you may already have
        written about. The Adjusted Rand Index is the number usually quoted, but it
        hides which direction the change went.
    </p>
    <table class="data-table" style="margin:14px 0 8px;border-collapse:collapse;">
        <thead><tr>
            <th style="{th_l}">Cutoffs</th>
            <th style="{th}">ARI</th>
            <th style="{th}">Pairs split</th>
            <th style="{th}">Pairs merged</th>
            <th style="{th}">Pairs kept</th>
        </tr></thead>
        <tbody>{pair_rows}</tbody>
    </table>
    <p style="color:#777;font-size:.88em;max-width:78ch;">
        A cutoff that moves the family count while leaving the largest family alone is
        splitting the tail rather than restructuring. Note BiG-SCAPE is not exactly
        reproducible at scale even at one cutoff &mdash; 3 of 1,302 members moved
        between two identical runs on Enterobacterales &mdash; so treat a handful of
        changed pairs as noise rather than signal.
    </p>
    '''
