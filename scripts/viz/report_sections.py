"""HTML section builders for the BGC report.

Each function returns a self-contained block of markup that generate_html_report()
drops into the page, keeping the report generator readable. The _build_* helpers are
pure string building; build_coupling_table_rows additionally reads the BiG-SCAPE DB so
the coupling table reflects the current run rather than hardcoded family IDs.
"""

import html as _html
import json as _json
import sqlite3

from utils.constants import load_coupling_classes, COUPLING_COLORS, KCB_THRESHOLDS
from utils.coupling_confidence import BACKGROUND_CEILING_PCT
from utils.gene_diagram import generate_gene_svg
from viz.report_assets import MAX_PANES


def _build_kcb_content(kcb_stats, taxon_clean, gcf_data, gcf_classes=None, gcf_hrefs=None):
    """Compute KCB tab contents.

    Returns dict with keys: kcb_mapping_section, novel_bgcs_tab_content.
    """
    kcb_mapping_section = ''''''
    novel_bgcs_tab_content = '''
            <h2>Detected BGC regions</h2>
            <p style="color: #666;">No KnownClusterBlast data available. Run antiSMASH with <code>--antismash_cb_knownclusters true</code> to record known-cluster matches.</p>'''

    if kcb_stats.get('total_regions', 0) > 0:
        # Every region, ranked or not, above the floor or below it. This used to be
        # two tables: one region-per-row for the below-floor regions and one
        # cluster-per-row for the rest. The second said nothing the first could not,
        # ordered it differently, and truncated its region lists at five.
        novel_bgcs = kcb_stats.get('all_regions') or kcb_stats.get('novel_bgcs', [])
        novel_count = kcb_stats.get('novel_bgc_count', 0)

        gcf_lookup = {}
        has_gcf_data = False
        if gcf_data and isinstance(gcf_data, dict):
            bgc_to_gcf = gcf_data.get('bgc_to_gcf', {})
            if bgc_to_gcf:
                for key_str, info in bgc_to_gcf.items():
                    parts = key_str.split('|')
                    if len(parts) == 2:
                        gcf_lookup[(parts[0], parts[1])] = info
                        has_gcf_data = True

        if novel_bgcs:
            detail_rows = ''
            for bgc in novel_bgcs:
                genome = bgc.get('genome', 'unknown')
                region = bgc.get('region', '?')
                region_name = bgc.get('region_name', region)
                record_index = bgc.get('record_index', 1)
                product = bgc.get('product', 'Unknown')
                contig_edge = bgc.get('contig_edge', '')
                # The best KnownClusterBlast hit, whatever its similarity. On this
                # chemistry every hit falls at or below the floor, so showing only
                # above-floor hits made a weak hit indistinguishable from none at all.
                # `kcb_*` is the hit that cleared the floor; `top_*` is the best hit
                # regardless. They are the same hit when one cleared.
                cleared = str(bgc.get('kcb_hit') or '')
                top_hit = str(bgc.get('top_hit') or '') or cleared
                top_sim = bgc.get('top_sim')
                acc = str(bgc.get('top_acc') or bgc.get('kcb_acc') or '')
                try:
                    sim_val = float(top_sim)
                except (TypeError, ValueError):
                    sim_val = None
                sim_txt = f'{sim_val:.0f}%' if sim_val is not None else ('?' if top_hit else '')
                if top_hit:
                    name = top_hit[:28] + ('…' if len(top_hit) > 28 else '')
                    tone = '#1d6fa5' if cleared else '#6c757d'
                    link = (f'<a href="https://mibig.secondarymetabolites.org/repository/{acc}" '
                            f'target="_blank" style="color:{tone};">{name}</a>' if acc else name)
                    if cleared:
                        note = (f'<span style="color:#1d6fa5;font-weight:600;" title="at or above '
                                f'the {KCB_THRESHOLDS["low"]}% floor">{sim_txt}</span>')
                    else:
                        note = (f'<span style="color:#999;" title="below the '
                                f'{KCB_THRESHOLDS["low"]}% floor — too weak to call this '
                                f'cluster known">{sim_txt}</span>')
                    kcb_cell = f'{link} {note}'
                else:
                    kcb_cell = '<span style="color:#ccc;">no hit</span>'
                # MIBiG accession in its own column, because it is the identifier a
                # reader takes elsewhere -- and it makes the column sortable, which
                # a name truncated to 28 characters inside a link is not.
                mibig_cell = (f'<a href="https://mibig.secondarymetabolites.org/repository/{acc}" '
                              f'target="_blank" style="color:#6c757d;font-variant-numeric:'
                              f'tabular-nums;">{acc}</a>' if acc
                              else '<span style="color:#ccc;">—</span>')
                sort_sim = f'{sim_val:.4f}' if sim_val is not None else '-1'
                edge_badge = '<span style="background: #e74c3c; color: white; padding: 1px 5px; border-radius: 3px; font-size: 0.75em;">edge</span>' if str(contig_edge).lower() == 'true' else ''
                antismash_link = f'../../antismash_results/{taxon_clean}/{genome}/index.html#r{record_index}c{region}'
                gcf_cell = ''
                if has_gcf_data:
                    gcf_info = gcf_lookup.get((genome, str(region_name)), {})
                    if gcf_info:
                        fid = gcf_info.get('family_id', '')
                        mc = gcf_info.get('member_count', 1)
                        # Badge colour = the GCF's dominant coupling enzyme class,
                        # the same palette the GCF tree and heatmap legends use, so a
                        # red GCF-1 here is the red GCF-1 in the tree.
                        cls = (gcf_classes or {}).get(fid)
                        badge_bg = COUPLING_COLORS.get(cls, '#3498db')
                        badge_title = f' title="Coupling class: {cls}"' if cls else ''
                        # The badge links to the family's own page, which carries the
                        # representative cluster and the consensus gene content.
                        href = (gcf_hrefs or {}).get(str(fid), '')
                        badge = (f'<span{badge_title} style="background: {badge_bg}; color: white; '
                                 f'padding: 2px 8px; border-radius: 4px; font-size: 0.85em; '
                                 f'display: inline-block; white-space: nowrap;">GCF-{fid}</span>')
                        if href:
                            badge = f'<a href="{href}" style="text-decoration: none;">{badge}</a>'
                        gcf_cell = (f'<td style="text-align: center;" data-sort="{int(fid):06d}">'
                                    f'{badge}</td>'
                                    f'<td style="text-align: center;" data-sort="{int(mc):06d}">{mc}</td>')
                    else:
                        gcf_cell = ('<td style="text-align: center; color: #999;" data-sort="~">-</td>'
                                    '<td style="text-align: center; color: #999;" data-sort="-1">-</td>')
                # data-sort carries the sort key for the columns where the rendered
                # text does not sort correctly: a % inside a span, an accession beside
                # a link, a GCF badge whose text is "GCF-7" and must order as 7.
                detail_rows += f'''
                <tr>
                    <td><a href="genomes/{genome}.html" title="{genome}" style="color: #2c5aa0; display: inline-block; max-width: 340px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; vertical-align: middle;">{genome}</a></td>
                    <td style="text-align: center;"><a href="{antismash_link}" target="_blank" style="color: #28a745; font-weight: bold;">Region {region_name}</a></td>
                    <td>{product}</td>
                    <td style="text-align: center;" data-sort="{'1' if str(contig_edge).lower() == 'true' else '0'}">{edge_badge}</td>
                    {gcf_cell}
                    <td style="font-size: .85em;" data-sort="{_html.escape(top_hit.lower(), quote=True) or '~'}">{kcb_cell}</td>
                    <td style="font-size: .85em;" data-sort="{acc or '~'}">{mibig_cell}</td>
                    <td style="text-align: right; font-size: .85em;" data-sort="{sort_sim}">{sim_txt}</td>
                </tr>'''

            kcb_mapping_section = ''''''
            th = ('onclick="sortRegions(this)" style="cursor:pointer;" '
                  'title="Click to sort; click again to reverse"')
            gcf_header = (f'<th {th}>GCF Family</th><th {th}>Members</th>'
                          if has_gcf_data else '')
            gcf_description = (' Follow a GCF badge to that family’s page.'
                               if has_gcf_data else '')
            kcb_floor = KCB_THRESHOLDS['low']
            n_cleared = sum(1 for b in novel_bgcs if b.get('kcb_hit'))
            novel_bgcs_tab_content = f'''
            <h2>Detected BGC regions</h2>
            <p style="color: #666; margin-bottom: 15px; max-width: 78ch;">
                <em>Every region found. <strong>Sort by any column</strong> — sorting on
                <strong>MIBiG ID</strong> or <strong>Best KCB hit</strong> groups the regions that
                matched the same characterised cluster. Of {len(novel_bgcs):,} regions,
                {n_cleared:,} have a KnownClusterBlast hit at or above the {kcb_floor}% floor
                (shown in blue); the rest are greyed, because a hit below the floor is too weak
                to call the cluster known, and a region with no hit at all says "no hit". For
                phosphonate chemistry a miss is weak evidence of novelty — MIBiG holds few
                characterised pathways. "edge" marks a region on a contig boundary, which may be
                incomplete.{gcf_description}</em>
            </p>
            <div class="search-box">
                <input type="text" id="novelSearch" placeholder="Search by genome, strain, region, GCF or MIBiG ID (e.g. BGC0000904)" onkeyup="filterNovelBGCs()">
            </div>
            <div class="table-container">
                <table id="novelTable">
                    <thead>
                        <tr>
                            <th {th}>Genome</th>
                            <th {th}>antiSMASH Region</th>
                            <th {th}>Product Type</th>
                            <th {th}>Contig Edge</th>
                            {gcf_header}
                            <th {th}>Best KCB hit</th>
                            <th {th}>MIBiG ID</th>
                            <th {th} style="cursor:pointer;text-align:right;">Similarity</th>
                        </tr>
                    </thead>
                    <tbody id="novelTableBody">
                        {detail_rows}
                    </tbody>
                </table>
            </div>'''
        elif kcb_stats.get('total_regions', 0) > 0:
            kcb_mapping_section = ''''''
            novel_bgcs_tab_content = '''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666;">All detected BGC regions matched characterized clusters in the MIBiG database. No potentially novel BGCs identified.</p>'''
        else:
            novel_bgcs_tab_content = '''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666;">No KnownClusterBlast data available. Run antiSMASH with <code>--antismash_cb_knownclusters true</code> to identify potentially novel BGCs.</p>'''

    return {
        'kcb_mapping_section':    kcb_mapping_section,
        'novel_bgcs_tab_content': novel_bgcs_tab_content,
    }


def _stat_tile(value, label, sub=''):
    """One overview tile. Uniform surface — the previous grid gave each of 13 tiles a
    different hue, which encoded nothing and competed with the figures below."""
    sub_html = (f'<div style="font-size: 0.78em; margin-top: 3px; color: #888;">{sub}</div>'
                if sub else '')
    return (f'<div style="background: #f6f7f8; border: 1px solid #e3e6e8; '
            f'padding: 11px 14px; border-radius: 8px;">'
            f'<div style="font-size: 1.5em; font-weight: 600; color: #2c3e50;">{value}</div>'
            f'<div style="color: #555; font-size: 0.85em;">{label}</div>'
            f'{sub_html}</div>')


def _stat_group(title, tiles):
    """A labelled band of tiles. Grouping replaces a flat 13-tile grid in which nothing
    signalled what mattered or how the numbers related."""
    if not tiles:
        return ''
    return (f'<div style="margin-top: 18px;">'
            f'<div style="font-size: 0.78em; text-transform: uppercase; letter-spacing: 0.06em; '
            f'color: #8a8a8a; margin-bottom: 7px;">{title}</div>'
            f'<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); '
            f'gap: 10px;">' + ''.join(tiles) + '</div></div>')


def build_overview_stats(stats, kcb_stats, gcf_data, rarefaction_stats=None):
    """The Overview tab's summary, grouped into sampling / BGCs / diversity.

    Replaces a flat grid of 13 tiles that carried roughly 7 distinct facts: genomes
    with and without a BGC are both derivable from the total, "% matching known
    clusters" and "potentially novel" were the same fact stated inversely, and "most
    common BGC type" is constant because detection is hardcoded to the phosphonate
    rule. Chao2 coverage — arguably the headline result — was prose beneath the chart
    rather than a figure anyone would see.
    """
    total = stats.get('total_genomes', 0) or 0
    with_bgc = stats.get('genomes_with_bgcs', 0) or 0
    total_bgcs = stats.get('total_bgcs', 0) or 0
    pct = f'{100.0 * with_bgc / total:.1f}%' if total else '—'

    sampling = [
        _stat_tile(f'{total:,}', 'Genomes analysed'),
        _stat_tile(f'{with_bgc:,}', 'Carry a phosphonate BGC', f'{pct} of those analysed'),
    ]

    # Per BGC-positive genome, not per genome analysed: averaging across genomes with
    # no BGC mostly restates how many lack one, which the tile above already says.
    per_positive = f'{total_bgcs / with_bgc:.2f}' if with_bgc else '—'
    bgcs = [
        _stat_tile(f'{total_bgcs:,}', 'BGC regions'),
        _stat_tile(per_positive, 'Per BGC-positive genome',
                   f"range 1–{stats.get('max_bgcs', 0)}"),
        _stat_tile(f"{kcb_stats.get('contig_edge_count', 0):,}", 'On a contig edge',
                   'possibly incomplete'),
        _stat_tile(f"{kcb_stats.get('novel_bgc_count', 0):,}", 'No MIBiG match',
                   'all appear under Novel BGCs'),
    ]

    diversity = []
    if gcf_data and isinstance(gcf_data, dict) and gcf_data.get('gcfs'):
        gcfs = gcf_data['gcfs']
        summary = gcf_data.get('summary', {})
        total_gcfs = summary.get('total', len(gcfs))
        singletons = summary.get('singletons',
                                 sum(1 for g in gcfs if g.get('is_singleton')))
        largest = max((g.get('member_count', 0) for g in gcfs), default=0)
        diversity.append(_stat_tile(total_gcfs, 'Gene cluster families',
                                    f'{singletons} seen in one genome'))
        diversity.append(_stat_tile(f'{largest:,}', 'Largest family'))
    if rarefaction_stats and rarefaction_stats.get('chao2'):
        c = rarefaction_stats['chao2']
        diversity.append(_stat_tile(f"{c.get('coverage', 0):.0f}%", 'Coverage (Chao2)',
                                    f"{c.get('s_est', 0):.0f} families estimated"))

    return (_stat_group('Sampling', sampling)
            + _stat_group('BGCs', bgcs)
            + _stat_group('Diversity', diversity))


def _build_gcf_support_section(gcf_support_rows):
    """GCF Analysis block pairing GCF membership with coupling-call support."""
    if not gcf_support_rows:
        return ''
    return f'''
            <h3 style="margin-top: 30px;">Coupling enzyme support by GCF</h3>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 12px; max-width: 74ch;">
                Membership comes from BiG-SCAPE (whole gene neighbourhood); the coupling
                class comes from antiSMASH domain markers. <strong>Support</strong> is the
                identity of the enzyme that drove the call, measured against the nearest
                characterised reference of its class — advisory, not a verdict. A low value
                may mean the call is wrong, or that the enzyme is a novel variant; both
                warrant a look.
            </p>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 12px; max-width: 74ch;">
                <strong>Read it against the reference named.</strong> Enzymes of
                <em>different</em> classes score 26.7–29.7% against each other, so at or
                below ~30% an identity carries no class information. Above that there is no
                threshold to quote: five of the seven references are <em>Streptomyces</em>,
                so a modest score may be the genus gap rather than a doubtful call, and the
                one class with a <em>Pantoea</em> reference scores near 100% for the same
                reason.
            </p>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>GCF</th><th>Dominant coupling class</th><th>Members</th>
                            <th>Median support</th><th>Range</th><th>Nearest reference</th>
                            <th>Refs</th><th></th>
                        </tr>
                    </thead>
                    <tbody>
{gcf_support_rows}
                    </tbody>
                </table>
            </div>'''


def build_bigscape_stats_section(bigscape_stats_html, taxon_clean):
    """Network-level clustering numbers, and where to open BiG-SCAPE's own output.

    The per-family tables that used to sit in this block are their own sections
    now; what is left is the network summary and the pointer to the interactive
    output, which is the one thing this report cannot reproduce.
    """
    if not bigscape_stats_html:
        return ''
    return f'''
            <div class="clustering-section">
                {bigscape_stats_html}
                <div class="info-box" style="margin-top: 20px;">
                    <p><strong>The interactive BiG-SCAPE output</strong> — the similarity
                    network, per-family alignments and the distance browser — is not
                    reproduced here. To open it:</p>
                    <ol style="margin: 10px 0 10px 20px; line-height: 1.8;">
                        <li>Open: <code>results/bigscape_results/{taxon_clean}/index.html</code></li>
                        <li>Select database: <code>results/bigscape_results/{taxon_clean}/{taxon_clean}.db</code></li>
                    </ol>
                </div>
            </div>'''


def _build_versions_html(versions_data):
    """Return HTML table for software versions."""
    if not versions_data:
        return ''
    tool_names = {
        'antismash': 'antiSMASH', 'bigscape': 'BiG-SCAPE', 'gtdbtk': 'GTDB-Tk',
        'taxonkit': 'TaxonKit', 'nextflow': 'Nextflow', 'pipeline_version': 'Pipeline Version',
    }
    version_rows = ''
    for tool, version in versions_data.items():
        display_name = tool_names.get(tool, tool.title())
        version_rows += f'''
                <tr>
                    <td style="padding: 8px 12px; border-bottom: 1px solid #eee; font-weight: 500;">{display_name}</td>
                    <td style="padding: 8px 12px; border-bottom: 1px solid #eee; font-family: monospace;">{version}</td>
                </tr>'''
    return f'''
            <div class="versions-section" style="margin-top: 30px;">
                <h3>Software Versions</h3>
                <p style="color: #666; margin-bottom: 15px;">
                    <em>Software versions used in this pipeline run for reproducibility.</em>
                </p>
                <table style="width: 100%; max-width: 500px; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">
                    <thead>
                        <tr style="background: #2c5aa0; color: white;">
                            <th style="padding: 10px 12px; text-align: left;">Software</th>
                            <th style="padding: 10px 12px; text-align: left;">Version</th>
                        </tr>
                    </thead>
                    <tbody>
                        {version_rows}
                    </tbody>
                </table>
            </div>'''


def _build_rarefaction_section(rarefaction_stats):
    """Return rarefaction curve HTML block for the Overview tab."""
    if not (rarefaction_stats and rarefaction_stats.get('generated')):
        return ''
    _c = rarefaction_stats.get('chao2') or {}
    _denominator = rarefaction_stats.get('denominator')
    _coverage_note = (
        f"Chao2 estimates <strong>{_c.get('s_est', 0):.0f}</strong> gene cluster families in this clade, of which "
        f"<strong>{rarefaction_stats.get('total_gcfs', 0)}</strong> were observed — "
        f"<strong>{_c.get('coverage', 0):.0f}% coverage</strong> "
        f"({_c.get('q1', 0)} GCFs seen in a single genome, {_c.get('q2', 0)} in exactly two)."
        if _c else ''
    )
    _axis_note = (
        f"Sampled over all {rarefaction_stats.get('n_genomes', 0)} analysed genomes, "
        f"{rarefaction_stats.get('n_bgc_positive', 0)} of which carry a BGC."
        if _denominator == 'analysed'
        else "Sampled over BGC-positive genomes only — genomes with no BGC are not on the axis, "
             "so this overstates coverage of the full genome set."
    )
    _rarefaction_svg_b64 = rarefaction_stats.get('svg_b64')
    _rarefaction_img = (
        f'<img src="data:image/svg+xml;base64,{_rarefaction_svg_b64}" '
        f'alt="GCF Rarefaction Curve" style="max-width: 100%; height: auto;">'
        if _rarefaction_svg_b64
        else '<img src="rarefaction_curve.png" alt="GCF Rarefaction Curve" style="max-width: 100%; height: auto;">'
    )
    return f'''
            <!-- GCF Rarefaction Curve -->
            <div class="plot" style="max-width: 800px; margin: 30px auto;">
                <h3 style="margin-top: 0;">GCF Rarefaction Curve</h3>
                <p style="color: #666; font-size: 0.9em; margin-bottom: 15px;">
                    Shows how the number of unique Gene Cluster Families (GCFs) increases as more genomes are sampled.
                    A plateau indicates sampling saturation — additional genomes would yield diminishing discovery of new GCFs.
                    {_axis_note}
                </p>
                {_rarefaction_img}
                <p style="color: #666; font-size: 0.9em; margin-top: 12px;">{_coverage_note}</p>
            </div>'''


# Coupling class display metadata (class_id → (display_name, marker, pathway, reference_genes))
# class_id -> (display name, marker, product, reference genes)
# The product column names what the enzyme makes from phosphonopyruvate. It used to be
# written with arrows ("→ phosphonomethylmalate → phosphinothricin-type"), which left a
# dangling arrow at the start of every cell; the immediate product is now named plainly
# with the downstream chemistry in parentheses.
_COUPLING_META = {
    'Synthase': (
        'Synthase', 'SMCOG1271 (HMGL-like)',
        'Phosphonomethylmalate (phosphinothricin-type)', 'FrbC, HvrC'),
    'Reductase': (
        'Reductase', 'Fe-ADH',
        'Phosphonolactate', 'VlpB'),
    'Decarboxylase': (
        'Decarboxylase', 'SMCOG1055 (ThDP)',
        '2-Phosphonoacetaldehyde (2-AEP)', 'DhpF, Fom2, Ppd'),
    'Transaminase': (
        # PnaA is the sequence the Support column is scored against; PalB names the
        # class in the literature but is not in assets/reference_sequences/, so
        # listing it bare implied an identity that is never computed.
        'Transaminase', 'SMCOG1019 (Aminotran_1_2 / PF00155)',
        'L-Phosphonoalanine', 'PnaA (PalB-like)'),
    'Unknown': (
        'Unknown', 'no marker matched',
        'not assignable', 'none'),
}
_COUPLING_ROW_ORDER = ['Synthase', 'Reductase', 'Decarboxylase', 'Transaminase', 'Unknown']

def gcf_coupling_classes(coupling_annotation_path, bigscape_db_path, cutoff=0.3):
    """Map each GCF id to its dominant coupling enzyme class.

    The single source of truth for "what class is this GCF" — used both by the
    coupling table and by the GCF badges in the Novel BGCs table, so the badge
    colour always agrees with the tree legend and the table.

    Returns {family_id: class_name}, or {} when the inputs are unavailable.
    """
    import os as _os
    from collections import Counter as _Counter
    try:
        coupling_classes = load_coupling_classes(str(coupling_annotation_path), region_only=True)
        conn = sqlite3.connect(str(bigscape_db_path))
        cur = conn.cursor()
        cur.execute("""
            SELECT f.id
            FROM family f
            JOIN bgc_record_family rf ON rf.family_id = f.id
            JOIN bgc_record br        ON br.id = rf.record_id
            WHERE f.cutoff = ? AND br.record_type = 'region'
            GROUP BY f.id
        """, (cutoff,))
        gcf_ids = [r[0] for r in cur.fetchall()]

        out = {}
        for gcf_id in gcf_ids:
            cur.execute("""
                SELECT g.path
                FROM bgc_record_family rf
                JOIN bgc_record br ON br.id = rf.record_id
                JOIN gbk g         ON g.id = br.gbk_id
                WHERE rf.family_id = ? AND br.record_type = 'region'
            """, (gcf_id,))
            counts = _Counter()
            for (path,) in cur.fetchall():
                gbk_base = _os.path.splitext(_os.path.basename(path))[0]
                counts[coupling_classes.get(gbk_base, 'Unknown')] += 1
            out[gcf_id] = counts.most_common(1)[0][0] if counts else 'Unknown'
        conn.close()
        return out
    except Exception as e:
        print(f"Warning: could not derive GCF coupling classes: {e}")
        return {}


def build_gcf_support_rows(coupling_support_path, coupling_annotation_path,
                           bigscape_db_path, cutoff=0.3):
    """Per-GCF table rows: dominant coupling class and how well evidenced it is.

    Pairs the two signals the classification now rests on — GCF membership from
    BiG-SCAPE (whole gene neighbourhood) and coupling class from SMCOG markers — and
    shows the reference support behind the second so a reader can see which calls are
    thinly evidenced.

    Deliberately a table rather than a chart: with ~13 GCFs the numbers matter more
    than the shape, and the report already carries several embedded images.

    Returns HTML <tr> rows, or None when the inputs are unavailable.
    """
    import csv
    import statistics
    from collections import defaultdict
    try:
        gcf_class = gcf_coupling_classes(coupling_annotation_path, bigscape_db_path, cutoff)
        if not gcf_class:
            return None

        # BGC label -> GCF, from the BiG-SCAPE record/family mapping
        conn = sqlite3.connect(str(bigscape_db_path))
        cur = conn.cursor()
        cur.execute("""
            SELECT g.path, rf.family_id
            FROM bgc_record_family rf
            JOIN bgc_record br ON br.id = rf.record_id
            JOIN gbk g         ON g.id = br.gbk_id
            JOIN family f      ON f.id = rf.family_id
            WHERE f.cutoff = ? AND br.record_type = 'region'
        """, (cutoff,))
        import os as _os
        bgc_to_gcf = {_os.path.splitext(_os.path.basename(p))[0]: fid
                      for p, fid in cur.fetchall()}
        conn.close()

        with open(coupling_support_path) as f:
            lines = [ln for ln in f if not ln.startswith('#')]
        per_gcf = defaultdict(list)
        n_refs = {}
        ref_of = {}
        for row in csv.DictReader(lines, delimiter='\t'):
            fid = bgc_to_gcf.get(row['bgc'])
            if fid is None:
                continue
            try:
                per_gcf[fid].append(float(row['assigned_pct_id']))
            except ValueError:
                continue
            n_refs[fid] = row.get('assigned_n_refs', '?')
            org = row.get('assigned_ref_organism', '') or ''
            ref = row.get('assigned_ref', '') or '—'
            # Binomials are italicised by convention; the gene name is not
            ref_of[fid] = (f'{ref} <em style="color:#666;">({org})</em>'
                           if org and org != '-' else ref)
        if not per_gcf:
            return None

        rows_html = []
        for i, fid in enumerate(sorted(per_gcf, key=lambda k: -len(per_gcf[k]))):
            vals = per_gcf[fid]
            cls = gcf_class.get(fid, 'Unknown')
            med = statistics.median(vals)
            colour = COUPLING_COLORS.get(cls, '#999999')
            # One boundary, and only one, is defensible from the current reference
            # set: characterised enzymes of *different* classes score 26.7-29.7%
            # against each other, so at or below ~30% an identity carries no class
            # information. Above it there are too few references, from too few genera,
            # to calibrate a verdict — so none is offered. The previous three tiers
            # ("well evidenced" >=90, "moderate" >=50, else "weak") were invented: no
            # GCF ever fell in the 50-90% band, and the 90% mark simply tracked
            # whether a class happened to have a same-genus reference.
            note = ('at superfamily background' if med <= BACKGROUND_CEILING_PCT else '')
            bg = ' style="background:#fafafa;"' if i % 2 else ''
            td = 'padding: 7px 12px; border-bottom: 1px solid #eee;'
            rows_html.append(
                f'<tr{bg}>'
                f'<td style="{td}"><span style="background: {colour}; color: white; '
                f'padding: 2px 8px; border-radius: 4px; font-size: 0.85em; '
                f'display: inline-block; white-space: nowrap;">GCF-{fid}</span></td>'
                f'<td style="{td}">{cls}</td>'
                f'<td style="{td} text-align: center;">{len(vals)}</td>'
                f'<td style="{td} text-align: center;">{med:.1f}%</td>'
                f'<td style="{td} text-align: center;">{min(vals):.1f}–{max(vals):.1f}%</td>'
                f'<td style="{td}">{ref_of.get(fid, "—")}</td>'
                f'<td style="{td} text-align: center;">{n_refs.get(fid, "?")}</td>'
                f'<td style="{td} color: #666;">{note}</td>'
                f'</tr>')
        return '\n'.join(rows_html)
    except Exception as e:
        print(f"Warning: could not build GCF support table: {e}")
        return None


def build_coupling_table_rows(coupling_annotation_path, bigscape_db_path, cutoff=0.3):
    """Return HTML <tr> rows for the coupling enzyme table, driven by live data.

    Falls back to the static hardcoded rows when inputs are unavailable.
    """
    try:
        from collections import defaultdict as _defaultdict
        gcf_class = gcf_coupling_classes(coupling_annotation_path, bigscape_db_path, cutoff)
        if not gcf_class:
            return None
        class_to_gcfs = _defaultdict(list)
        for gcf_id, dominant in gcf_class.items():
            class_to_gcfs[dominant].append(gcf_id)

        # Sort GCF IDs within each class and build HTML rows
        rows_html = []
        for i, cls_id in enumerate(_COUPLING_ROW_ORDER):
            gcf_ids = sorted(class_to_gcfs.get(cls_id, []))
            gcf_label = ', '.join(str(g) for g in gcf_ids) if gcf_ids else '—'
            display, marker, pathway, refs = _COUPLING_META[cls_id]
            bg = ' style="background:#fafafa;"' if i % 2 == 1 else ''
            is_last = (i == len(_COUPLING_ROW_ORDER) - 1)
            border = '' if is_last else 'border-bottom: 1px solid #eee; '
            td = f'padding: 7px 12px; {border}'
            swatch = (f'<span style="display: inline-block; width: 10px; height: 10px; '
                      f'border-radius: 2px; background: {COUPLING_COLORS.get(cls_id, "#999999")}; '
                      f'margin-right: 8px; vertical-align: middle;"></span>')
            rows_html.append(
                f'<tr{bg}>'
                f'<td style="{td} white-space: nowrap;">{swatch}{display}</td>'
                f'<td style="{td}">{marker}</td>'
                f'<td style="{td}">{pathway}</td>'
                f'<td style="{td}">{refs}</td>'
                f'<td style="{td}">{gcf_label}</td>'
                f'</tr>'
            )
        return '\n'.join(rows_html)

    except Exception as e:
        print(f"Warning: could not build dynamic coupling table ({e}); using static fallback")
        rows = []
        for i, cls_id in enumerate(_COUPLING_ROW_ORDER):
            display, marker, pathway, refs = _COUPLING_META[cls_id]
            bg = ' style="background:#fafafa;"' if i % 2 == 1 else ''
            is_last = (i == len(_COUPLING_ROW_ORDER) - 1)
            border = '' if is_last else 'border-bottom: 1px solid #eee; '
            td = f'padding: 7px 12px; {border}'
            rows.append(
                f'<tr{bg}>'
                f'<td style="{td}">{display}</td>'
                f'<td style="{td}">{marker}</td>'
                f'<td style="{td}">{pathway}</td>'
                f'<td style="{td}">{refs}</td>'
                f'<td style="{td}">—</td>'
                f'</tr>'
            )
        return '\n'.join(rows)


def build_partition_section(pepm_summary):
    """Whether pepM identity could partition BiG-SCAPE, for the pipeline info tab.

    Operational rather than biological: it reports how far the all-pairs
    clustering problem could be split without separating BGCs that belong in one
    family. Lossless means no pair BiG-SCAPE grouped would be cut apart.
    """
    parts = (pepm_summary or {}).get('partitioning')
    if not parts:
        return ''
    rows = ''.join(
        f'<tr><td style="padding:5px 10px;">{p["threshold"]:.2f}</td>'
        f'<td style="padding:5px 10px;text-align:right;">{p["components"]:,}</td>'
        f'<td style="padding:5px 10px;text-align:right;">{p["largest"]:,}</td>'
        f'<td style="padding:5px 10px;text-align:right;">{p["largest_fraction"]:.0%}</td>'
        f'<td style="padding:5px 10px;text-align:right;">{p["relative_work"]:.0%}</td>'
        f'<td style="padding:5px 10px;text-align:right;'
        f'{"color:#7a3;" if p["lossless"] else "color:#c33;font-weight:600;"}">'
        f'{"lossless" if p["lossless"] else str(p["same_gcf_pairs_lost"]) + " lost"}</td></tr>'
        for p in parts)
    return f'''
    <div class="section">
        <h3>BiG-SCAPE Partitioning Feasibility</h3>
        <p style="color:#555;max-width:70ch;">
            BiG-SCAPE compares every BGC against every other, so its memory grows
            quadratically. Splitting the input by pepM identity first can avoid that.
            This reports, for this dataset, how far it could be split and whether doing so
            would separate BGCs that BiG-SCAPE placed in the same family.
            <strong>Work</strong> is the resulting compute as a share of one unsplit job.
        </p>
        <table style="border-collapse:collapse;font-size:0.9em;">
            <thead><tr style="background:#f6f7f8;">
                <th style="padding:5px 10px;text-align:left;">pepM cut</th>
                <th style="padding:5px 10px;">partitions</th>
                <th style="padding:5px 10px;">largest</th>
                <th style="padding:5px 10px;">share</th>
                <th style="padding:5px 10px;">work</th>
                <th style="padding:5px 10px;">families</th>
            </tr></thead>
            <tbody>{rows}</tbody>
        </table>
        <p style="color:#777;font-size:0.85em;margin-top:10px;max-width:70ch;">
            &ldquo;Lossless&rdquo; is a lower bound rather than a guarantee: it counts pairs
            above the GCF similarity cutoff, while BiG-SCAPE families are transitively
            closed and include pairs below it. Treat any non-zero loss as disqualifying.
            Partitioning is enabled with <code>--bigscape_partition</code> and is off by
            default.
        </p>
    </div>
    '''


# ─── Sidebar navigation ────────────────────────────────────────────────────────

def build_tabbed_nav(groups):
    """Sidebar entries and content panes for the whole report.

    `groups` is [(group_label, [(entry_label, pane_html), ...]), ...].

    The report grew by merging tabs — the five that were left held between one and
    seven unrelated sections each, and Gene Cluster Families was a single pane about
    two metres long in which the only way to reach the family trees was to scroll
    past every consensus table. Each section is now its own pane, and the rail is
    what makes them reachable.

    A pane with no HTML is dropped and so is a group left empty by that, so a run
    without clustering has fewer entries rather than entries that open onto nothing.
    A group with one pane renders as a plain top-level entry; a group with several
    renders as a heading with indented entries beneath it.

    Numbering runs across the whole report, not per group, because the CSS rule that
    reveals a pane pairs `#tabN` with `#contentN`. Nothing may depend on a particular
    N: which sections a run emits decides them.
    """
    nav_parts, pane_parts, n = [], [], 0
    for group_label, entries in groups:
        live = [(label, html) for label, html in entries if html and html.strip()]
        if not live:
            continue
        if len(live) > 1:
            nav_parts.append(f'        <div class="nav-group">{_html.escape(group_label)}</div>')
        for label, html in live:
            n += 1
            if n > MAX_PANES:
                raise ValueError(
                    f'report has more than {MAX_PANES} panes; raise MAX_PANES in '
                    'viz/report_assets.py, which generates one CSS rule per pane')
            checked = ' checked' if n == 1 else ''
            cls = ' class="sub"' if len(live) > 1 else ''
            nav_parts.append(
                f'        <input type="radio" id="tab{n}" name="tabs"{checked}>\n'
                f'        <label for="tab{n}"{cls}>{_html.escape(label)}</label>')
            pane_parts.append(
                f'        <div class="tab-content" id="content{n}">\n{html}\n        </div>')
    return '\n'.join(nav_parts) + '\n\n' + '\n'.join(pane_parts)


# ─── Tab bodies ────────────────────────────────────────────────────────────────
# These were inline in a single 300-line f-string inside generate_html_report,
# which made every tab edit a careful string match into a wall of markup. Each
# returns the inner HTML of one `.tab-content` div; the surrounding div and the
# tab nav stay in visualize_results.py, where the tab numbering lives.

_MISSING = ('<div style="color: #999; padding: 20px; background: #f8f9fa; '
            'border-radius: 8px; text-align: center; font-size: 0.9em;">{}</div>')


def build_biosynthetic_phylogeny_section(coupling_table_rows):
    """The coupling-enzyme classification reference table.

    This is the key to every coupling class named elsewhere in the report — the
    tree colours, the support table, the branch points — so it is a reference
    page rather than a result. The trees themselves are under Family trees, and
    carry their own colour legends.
    """
    th = 'text-align: left; padding: 8px 12px; border-bottom: 2px solid #dee2e6;'
    return f'''
            <h3>GCF Biosynthetic Phylogeny</h3>
            <p style="color: #666; margin-bottom: 20px; max-width: 74ch;">
                <em>Phosphonate BGCs classified by the coupling enzyme acting on
                phosphonopyruvate — the branching step immediately downstream of PEP
                mutase, which fixes the rest of the pathway. This table is the key to
                every coupling class named elsewhere in the report.</em>
            </p>

            <div style="margin-top: 24px; background: #f8f9fa; padding: 20px; border-radius: 10px;">
                <h4 style="margin-top: 0;">Coupling Enzyme Classes</h4>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.9em;">
                    <thead>
                        <tr style="background: #e9ecef;">
                            <th style="{th}">Class</th>
                            <th style="{th}">Marker</th>
                            <th style="{th}">Product</th>
                            <th style="{th}">Reference genes</th>
                            <th style="{th}">GCFs</th>
                        </tr>
                    </thead>
                    <tbody>
                        {coupling_table_rows}
                    </tbody>
                </table>
            </div>'''


def build_no_clustering_notice():
    """Shown in place of every BiG-SCAPE-derived section when there was no run."""
    return ('<div class="info-box" style="background-color: #f8f9fa; '
            'border-left: 4px solid #6c757d;"><p style="color: #666;">'
            'No clustering analysis was performed. To enable clustering, run the '
            'pipeline with <code>--clustering bigscape</code>.</p></div>')


def build_gcf_trees_tab(gcf_tree_b64, gcf_tree_mime):
    """The family-centre tree.

    Every centre-to-centre distance is measured — via BIGSCAPE_CENTERS on a
    partitioned run — so nothing here is substituted, and the figure means the
    same thing in both run modes.

    Per-partition trees were tried here and removed. A partition is a pepM
    identity component sized to bound BiG-SCAPE's memory, not a biological
    grouping: on Erwiniaceae partition 0 held 236 BGCs but only 2 families, 215
    of them one family, so 91% of that figure was within-family variation drawn
    at a leaf count nobody can read. If per-BGC detail is wanted, the unit to
    draw is a family, not a partition.
    """
    centre = (f'<img src="data:{gcf_tree_mime};base64,{gcf_tree_b64}" '
              f'alt="GCF family-centre tree" '
              f'style="max-width: 100%; height: auto; display: block;">'
              if gcf_tree_b64 else
              _MISSING.format('Family-centre tree not generated.<br>'
                              'Run with <code>--clustering bigscape</code> to enable.'))
    return f'''
            <hr class="tab-section-divider">
            <h3>Gene Cluster Family Trees</h3>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Branch colours are the coupling enzyme classes tabulated above.</em>
            </p>

            <div class="plot" style="margin-top: 20px;">
                <h4 style="margin: 0 0 8px 0; color: #333;">Family-Centre Tree</h4>
                <p style="color: #666; font-size: 0.85em; margin: 0 0 12px 0;">
                    One representative (medoid) per Gene Cluster Family, circle size &prop; GCF membership.
                    Every centre-to-centre distance is measured, including across partitions,
                    so the backbone is real rather than substituted.
                </p>
                {centre}
            </div>'''


def build_pipeline_tab(resource_usage_html, partition_section_html, versions_html):
    """Nextflow resource usage, the partitioning feasibility table and versions."""
    no_trace = ('<div class="info-box warning"><p>No resource usage data available. '
                'Trace data will appear here after running the pipeline.</p></div>')
    return f'''
            <h2>Pipeline Information</h2>
            <h3 style="margin-top: 20px;">Resource Usage</h3>
            <p style="color: #666; margin-bottom: 20px; font-size: 0.9em;">
                <em>Resource consumption metrics from Nextflow trace data, showing CPU, memory, and runtime for each pipeline process.</em>
            </p>
            {resource_usage_html or no_trace}
            {partition_section_html}
            {versions_html}'''


def _twin_cell(other, differing):
    """The "Same chemistry as" cell: `=` only when the profiles really are equal.

    `differing` empty means the filtered profiles are identical and the split
    between the two families is not biosynthetic. Anything else means they are
    merely close, and the cell says so and names what differs -- a family that
    differs by one amide-bond ligase is a bench question, not a duplicate.
    """
    if not differing:
        return (f'<span style="color:#8a5a0c;" title="Same core/tailoring/lipid/'
                f'transport domains as GCF-{other} — the split between them is not '
                f'biosynthetic. Not counted against the rank.">= GCF-{other}</span>')
    return (f'<span style="color:#1d6fa5;" title="Nearly the same chemistry as '
            f'GCF-{other}, but NOT identical — differs by {differing}. '
            f'Worth checking before treating either as a duplicate. '
            f'Not counted against the rank.">~ GCF-{other}</span>')


def build_priority_section(ranking_path, bioprofile_path=None, gcf_hrefs=None,
                           branch_point_path=None, gcf_data=None):
    """Every family on one row: what it is, and how novel.

    This was two questions answered in one table of eleven columns, most of them
    inputs to the score rather than facts about the family. It is now one table of
    what a reader wants per family -- organism, known-cluster hit, coupling class,
    branch point, whether another family shares its chemistry -- and a second,
    collapsed, holding the score's components for anyone who wants to argue with
    the weighting.

    The "nearest reference" columns are gone. Every number in them was a non-hit:
    the closest characterised cluster to anything here sits at 0.37, against a
    family cutoff of 0.30, so the column reported the reference set's coverage
    rather than anything about the family, and read as a similarity when it was
    not one. Reference identity still zeroes a family whose chemistry is
    characterised; it is simply not shown as though it graded the rest.

    Each GCF links to its own page -- representative cluster, gene diagram, and the
    consensus gene content across every member. showGCF() is the fallback for a run
    whose per-family pages were not written.

    The weights are reasoned, not fitted: there is no set of leads that panned out
    to fit against. The components stay visible, one click away, so a reader can
    disagree with the weighting and re-order by eye.
    """
    import csv as _csv
    from pathlib import Path as _Path
    if not ranking_path or not _Path(ranking_path).exists():
        return ''

    # Families that are chemically indistinguishable from another, from
    # BIOSYNTHETIC_PROFILE. Shown as an ANNOTATION, never folded into the score.
    # Demoting on filtered-domain similarity would penalise exactly the families
    # this list exists to surface: the domain map covers ~93% of observed hits
    # and everything else is `other`, meaning unmapped rather than absent, so a
    # family whose chemistry nobody has curated would look like every other
    # empty profile. The flag tells the reader; it must not move a rank.
    # `=` (identical filtered profile) and `~` (near-identical) are kept apart.
    # Collapsing them put "same chemistry" on the pantaphos pair, which differs
    # by an ATP-grasp amide-bond ligase -- the one difference in this run worth
    # a bench experiment. A near miss names the domain instead.
    twin = {}
    if bioprofile_path and _Path(bioprofile_path).exists():
        for row in _csv.DictReader(_Path(bioprofile_path).open(), delimiter='\t'):
            if row.get('verdict'):
                twin[row['gcf']] = (row['nearest_gcf'],
                                    row.get('differing_domains', ''))
    # Branch point (2-AEP / 2-HEP / other) and the representative's organism and
    # known-cluster hit: all per-family facts that used to be spread over three
    # tables in two panes.
    branch = {}
    if branch_point_path and _Path(branch_point_path).exists():
        for row in _csv.DictReader(_Path(branch_point_path).open(), delimiter='\t'):
            branch[row['gcf']] = (row.get('branch_point', ''), row.get('support', ''))
    rep = {}
    for g in ((gcf_data or {}).get('gcfs') or []):
        rep[str(g.get('family_id'))] = (g.get('organism') or '',
                                        g.get('kcb_hit') or '', g.get('kcb_acc') or '')

    rows = list(_csv.DictReader(_Path(ranking_path).open(), delimiter='\t'))
    if not rows:
        return ''

    unc = [r for r in rows if r['status'] != 'ranked']
    ranked = [r for r in rows if r['status'] == 'ranked']

    # The case against ranking on reference identity is that the reference set is
    # small and taxonomically narrow, so identities clump rather than grading
    # smoothly. That was asserted here as "bimodal, nothing between 45% and 94%" --
    # an Erwiniaceae measurement that is simply untrue of Enterobacterales, where
    # the same column runs 18.5-100% with its widest gap at 75-94%. Measure it.
    pcts = sorted(float(r['reference_pct_id']) for r in rows
                  if (r.get('reference_pct_id') or '').strip())
    spread = ''
    if len(pcts) >= 4:
        lo, hi = pcts[0], pcts[-1]
        gap_lo, gap_hi = max(zip(pcts, pcts[1:]), key=lambda ab: ab[1] - ab[0])
        spread = f' Here they run {lo:.1f}–{hi:.1f}% across {len(pcts)} families'
        spread += (f', with nothing between {gap_lo:.1f}% and {gap_hi:.1f}%.'
                   if gap_hi - gap_lo >= 10 else '.')

    def bar(frac, tone):
        pct = max(0.0, min(1.0, frac)) * 100
        return (f'<div style="display:flex;align-items:center;gap:.45rem">'
                f'<div style="flex:0 0 46px;height:6px;background:#e9ecef;border-radius:3px;'
                f'overflow:hidden"><div style="width:{pct:.0f}%;height:100%;'
                f'background:{tone}"></div></div>'
                f'<span style="font-variant-numeric:tabular-nums">{frac:.2f}</span></div>')

    def gcf_link(fid, tone, title):
        """Link to the family's own page, or fall back to the in-report jump."""
        href = (gcf_hrefs or {}).get(str(fid))
        target = (f'href="{href}"' if href
                  else f'href="javascript:void(0)" onclick="showGCF(\'{fid}\')"')
        return (f'<a {target} style="color:{tone};font-weight:600;text-decoration:none;'
                f'border-bottom:1px dotted {tone};" title="{title}">GCF-{fid}</a>')

    unc_html = ''
    if unc:
        items = ''.join(
            f'<tr><td style="padding:6px 10px;white-space:nowrap;">'
            + gcf_link(r["gcf"], '#8a5a0c', 'Open this family\'s page')
            + '</td>'
            f'<td style="padding:6px 10px;text-align:right;">{r["members"]}</td>'
            f'<td style="padding:6px 10px;text-align:right;">{r["genomes"]}</td>'
            f'<td style="padding:6px 10px;text-align:right;">{float(r["intact"]):.0%}</td></tr>'
            for r in unc)
        unc_html = f'''
        <div style="background:#fdf6e3;border-left:3px solid #9a6b0f;padding:14px 18px;margin:0 0 22px;">
            <h4 style="margin:0 0 6px;">No clustering distances — {len(unc)} famil{'y' if len(unc)==1 else 'ies'}</h4>
            <p style="margin:0 0 10px;color:#555;font-size:.9em;max-width:66ch;">
                BiG-SCAPE produced no all-pairs distances for these, so there is nothing
                honest to rank them on. Worth a look by hand.
            </p>
            <table style="border-collapse:collapse;font-size:.9em;">
                <thead><tr style="background:#f4ead3;">
                    <th style="text-align:left;padding:5px 10px;">Family</th>
                    <th style="text-align:right;padding:5px 10px;">BGCs</th>
                    <th style="text-align:right;padding:5px 10px;">Genomes</th>
                    <th style="text-align:right;padding:5px 10px;">Intact</th>
                </tr></thead>
                <tbody>{items}</tbody>
            </table>
        </div>'''

    # Branch point is the one chemistry fact a reader can act on, so it is named
    # plainly and the hedged variants ("unknown (Ppd, no third enzyme found)") are
    # shortened to the call with the detail on hover.
    def branch_cell(gcf):
        bp, sup = branch.get(gcf, ('', ''))
        if not bp:
            return '<span style="color:#bbb;">—</span>'
        short = bp.split('(')[0].strip().rstrip(',')
        tone = '#0e5c6b' if not short.lower().startswith('unknown') else '#8a7a55'
        sup_txt = f' · {sup}' if sup else ''
        return (f'<span style="color:{tone};" title="{_html.escape(bp)}{sup_txt}">'
                f'{_html.escape(short)}</span>')

    def known_cell(gcf):
        _, hit, acc = rep.get(gcf, ('', '', ''))
        if not hit:
            return '<span class="tag-novel">none</span>'
        txt = _html.escape(hit[:26] + ('…' if len(hit) > 26 else ''))
        return (f'<a href="https://mibig.secondarymetabolites.org/repository/{acc}" '
                f'target="_blank" style="color:#1d6fa5;">{txt}</a>' if acc else txt)

    def org_cell(gcf):
        org = rep.get(gcf, ('', '', ''))[0]
        return f'<em>{_html.escape(org)}</em>' if org else '<span style="color:#bbb;">—</span>'

    body = ''.join(
        f'<tr>'
        f'<td style="padding:7px 10px;white-space:nowrap;">'
        + gcf_link(r["gcf"], '#2c5aa0',
                   'Open this family\'s page: representative cluster, gene diagram '
                   'and consensus gene content')
        + '</td>'
        f'<td style="padding:7px 10px;font-weight:600;text-align:right;'
        f'font-variant-numeric:tabular-nums;">{float(r["priority"]):.3f}</td>'
        f'<td style="padding:7px 10px;text-align:right;">{r["members"]}</td>'
        f'<td style="padding:7px 10px;font-size:.9em;">{org_cell(r["gcf"])}</td>'
        f'<td style="padding:7px 10px;font-size:.85em;">{known_cell(r["gcf"])}</td>'
        f'<td style="padding:7px 10px;font-size:.9em;">{r["coupling_class"]}</td>'
        f'<td style="padding:7px 10px;font-size:.85em;">{branch_cell(r["gcf"])}</td>'
        f'<td style="padding:7px 10px;font-size:.85em;">'
        + (_twin_cell(*twin[r["gcf"]]) if r["gcf"] in twin
           else '<span style="color:#bbb;">—</span>')
        + '</td>'
        f'</tr>'
        for r in ranked)

    # The score's components, collapsed. They justify the order rather than
    # describe the family, so they are one click away instead of eight columns
    # wide in the table a reader actually reads.
    score_rows = ''.join(
        f'<tr>'
        f'<td style="padding:6px 10px;color:#888;text-align:right;">{r["rank"]}</td>'
        f'<td style="padding:6px 10px;white-space:nowrap;">GCF-{r["gcf"]}</td>'
        f'<td style="padding:6px 10px;font-weight:600;text-align:right;'
        f'font-variant-numeric:tabular-nums;">{float(r["priority"]):.3f}</td>'
        f'<td style="padding:6px 10px;">{bar(float(r["distance"]), "#0e5c6b")}</td>'
        f'<td style="padding:6px 10px;">{bar(float(r["evidence"]), "#166b47")}</td>'
        f'<td style="padding:6px 10px;text-align:right;">{r["members"]}</td>'
        f'<td style="padding:6px 10px;text-align:right;">{r["genomes"]}</td>'
        f'<td style="padding:6px 10px;text-align:right;">{r["genera"]}</td>'
        f'<td style="padding:6px 10px;text-align:right;">{float(r["intact"]):.0%}</td>'
        f'</tr>'
        for r in ranked)

    return f'''
    <div class="section">
        <h3>Novelty Assessment</h3>
        <p style="color:#555;max-width:74ch;">
            One row per family, ordered by <strong>isolation × evidence</strong>.
            <strong>Isolation</strong> is how far a family sits from every other family
            in this run; <strong>evidence</strong> is independent genomes, independent
            genera, and the share of regions not truncated at a contig edge. They
            multiply because both are necessary. The weights are reasoned, not fitted,
            so treat the order as a considered opinion &mdash; the components are under
            <em>How the score is built</em> below.
        </p>
        {unc_html}
        <div class="table-container">
        <table style="width:100%;border-collapse:collapse;font-size:.9em;">
            <thead><tr style="background:#e9ecef;">
                <th style="text-align:left;padding:6px 10px;">Family</th>
                <th style="text-align:right;padding:6px 10px;" title="isolation × evidence">Novelty</th>
                <th style="text-align:right;padding:6px 10px;">BGCs</th>
                <th style="text-align:left;padding:6px 10px;">Representative organism</th>
                <th style="text-align:left;padding:6px 10px;"
                    title="Best KnownClusterBlast hit for the representative. MIBiG holds
                    few characterised phosphonate pathways, so "none" is weak evidence of
                    novelty.">Known cluster</th>
                <th style="text-align:left;padding:6px 10px;">Coupling class</th>
                <th style="text-align:left;padding:6px 10px;"
                    title="The intermediate the pathway branches through, from the enzymes
                    present. Hover a cell for the full call and its support.">Branch point</th>
                <th style="text-align:left;padding:6px 10px;" title="Another family with the
                    same core/tailoring/lipid/transport domains. A split between two such
                    families is not biosynthetic. Annotation only — it does not affect the
                    rank.">Same chemistry as</th>
            </tr></thead>
            <tbody>{body}</tbody>
        </table>
        </div>

        <details style="margin-top:20px;border:1px solid #dee2e6;border-radius:8px;">
            <summary style="cursor:pointer;padding:11px 16px;font-weight:600;
                            background:#f8f9fa;border-radius:8px;">
                How the score is built
            </summary>
            <div style="padding:6px 16px 16px;">
                <p style="color:#555;max-width:74ch;font-size:.92em;">
                    <strong>Isolation</strong> is measured over BiG-SCAPE&rsquo;s all-pairs
                    matrix, not against a reference: the characterised set is seven
                    proteins, five of them <em>Streptomyces</em>, and nothing in this run
                    falls inside the family cutoff of any of them.{spread} Reference
                    identity is used only to zero a family whose chemistry is already
                    characterised.
                </p>
                <div class="table-container">
                <table style="width:100%;border-collapse:collapse;font-size:.88em;">
                    <thead><tr style="background:#eef1f2;">
                        <th style="text-align:right;padding:6px 10px;">#</th>
                        <th style="text-align:left;padding:6px 10px;">Family</th>
                        <th style="text-align:right;padding:6px 10px;">Novelty</th>
                        <th style="text-align:left;padding:6px 10px;">Isolation</th>
                        <th style="text-align:left;padding:6px 10px;">Evidence</th>
                        <th style="text-align:right;padding:6px 10px;">BGCs</th>
                        <th style="text-align:right;padding:6px 10px;">Genomes</th>
                        <th style="text-align:right;padding:6px 10px;">Genera</th>
                        <th style="text-align:right;padding:6px 10px;">Intact</th>
                    </tr></thead>
                    <tbody>{score_rows}</tbody>
                </table>
                </div>
            </div>
        </details>
    </div>'''


def build_novelty_intro():
    """One line saying what the ranking orders on, above the ranking itself."""
    return '''
            <h2>BGC Novelty</h2>
            <p style="color:#666;max-width:72ch;">
                <em>Which gene cluster families are worth taking into the laboratory,
                and why.</em>
            </p>'''


def build_all_regions_section(all_regions_html, n_regions=0):
    """Every region found, as its own pane rather than a fold in the ranking.

    The ranking has one row per family and this has one per region — 72 against
    1,303 on Enterobacterales. Presenting them as equals buries the actionable part
    under a table where every row says the same thing, which is what the old "Novel
    BGCs" tab did at 30.7% of the whole report. It used to be hidden in a
    `<details>` under the ranking for that reason; a pane of its own separates them
    without hiding it.
    """
    if not all_regions_html:
        return ''
    # Count the rows actually rendered here. `n_regions` counts every region in the
    # tabulation, but this table holds only those below the KnownClusterBlast floor
    # -- 1,196 of 1,309 on Enterobacterales -- so using it labelled the table with a
    # total that did not match the rows under it.
    shown = max(0, all_regions_html.count('<tr') - all_regions_html.count('<thead'))
    count = f' ({shown:,})' if shown else ''
    return f'''
            <h2>All detected regions{count}</h2>
            <p style="color:#666;font-size:.9em;max-width:70ch;">
                Every region whose best KnownClusterBlast hit falls below the similarity
                floor — which for phosphonate chemistry is nearly all of them. The rest are
                under <strong>Known-cluster matches</strong>. Use this to locate a specific
                contig or genome; use <strong>Novelty assessment</strong> to decide what
                to work on.
            </p>
            {all_regions_html}'''


# ─── Consensus gene content ──────────────────────────────────────────────────

_ROLE_STYLE = {
    'core':       ('#0e5c6b', '#d9eef2', 'phosphonate pathway'),
    'tailoring':  ('#7a4b12', '#fbeedd', 'modifies the product'),
    'lipid':      ('#6a1b63', '#f7e4f6', 'lipid handling'),
    'transport':  ('#1a4f8a', '#e2ecf9', 'moves the product'),
    'regulation': ('#4a4a10', '#f3f2dd', 'controls expression'),
    'mobile':     ('#7a1f1f', '#fbe4e4', 'how the cluster arrived'),
    'catabolism': ('#1d5c2e', '#e2f3e6',
                   'degrades phosphonates — the C-P lyase operon, which scavenges '
                   'phosphorus rather than making a product'),
    'primary metabolism': ('#666666', '#eeeeee',
                           'central metabolism — a chromosomal neighbour, not part of the cluster'),
    'other':      ('#888888', '#f4f4f4', 'not classified'),
}


# Below this share of members a gene is noise in a consensus: present in a handful
# of regions, usually unnamed, and swept in at a region boundary rather than part of
# what the family is. GCF-27's table was 505 rows at no cutoff and is 45 at this one;
# across the run 3,263 rows become 2,214. The full list stays in
# gcf_consensus_clusters.tsv -- this is a display cut, not a filter on the analysis.
CONSENSUS_MIN_PREVALENCE = 0.10


def _in_gene_order(fam_rows):
    """Consensus genes in cluster order, matching the diagram drawn above them.

    Prevalence order was the old default and is wrong for reading a cluster: a BGC
    is a sequence, and the first thing anyone wants is left-to-right.

    Two orders are available and they disagree. The diagram is drawn on one member,
    so its order is that member's coordinates; `median_rank` is the family-wide
    consensus order. Sorting on median_rank alone put the table out of step with
    the diagram in 37 of 72 families here -- the same genes listed in a different
    sequence from the picture above them, which is worse than either order alone.

    So: genes the scaffold carries take their scaffold order, and the rest are
    slotted between them by interpolating their median_rank against the scaffold
    genes' ranks. One order, and it is the one drawn.

    A run predating these columns falls back to the prevalence sort.
    """
    def rank_of(r):
        v = (r.get('median_rank') or '').strip()
        return int(v) if v else None

    if not any(rank_of(r) is not None for r in fam_rows):
        return sorted(fam_rows, key=lambda r: -float(r['prevalence'] or 0))

    on = sorted((r for r in fam_rows if (r.get('scaffold_start') or '').strip()),
                key=lambda r: int(r['scaffold_start']))
    if not on:
        return sorted(fam_rows, key=lambda r: (rank_of(r) if rank_of(r) is not None else 1 << 30,
                                               -float(r['prevalence'] or 0)))

    pos = {id(r): float(i) for i, r in enumerate(on)}
    # (median_rank, scaffold position) for the drawn genes, to interpolate against.
    anchors = sorted((rank_of(r), pos[id(r)]) for r in on if rank_of(r) is not None)

    def slot(r):
        k = rank_of(r)
        if k is None:
            return len(on) + 0.5          # no rank at all: after everything drawn
        lo = max((a for a in anchors if a[0] <= k), default=None)
        hi = min((a for a in anchors if a[0] >= k), default=None)
        if lo is None:
            return hi[1] - 0.5
        if hi is None:
            return lo[1] + 0.5
        if lo[0] == hi[0]:
            return lo[1] + 0.5            # same rank as a drawn gene: just after it
        frac = (k - lo[0]) / (hi[0] - lo[0])
        return lo[1] + frac * (hi[1] - lo[1])

    return sorted(fam_rows,
                  key=lambda r: (pos.get(id(r), slot(r)),
                                 0 if id(r) in pos else 1,
                                 -float(r['prevalence'] or 0)))


# A gene in one of these roles is never hidden, whatever its prevalence. Learned
# the hard way: the confirmed phosphonolipid's defining enzyme sits at prevalence
# 0.16 in GCF-27, not because it is accessory but because that family is three
# architectures averaged together. Low prevalence can mean "this family is
# chimeric" as easily as "this gene is rare", and the display must not decide.
NEVER_HIDE = {'core', 'tailoring', 'lipid', 'transport', 'catabolism'}


def _member_sets(per_cds_path):
    """{family: {group: {region, ...}}} from gcf_annotation_transfer.tsv.

    Which members carry each orthologue group, which is what tells one insertion
    apart from several independent rare genes. Exact, from the group id the
    transfer records -- matching on product text mis-assigned the GCF-22
    pathogenicity island, putting STY4528 with the ParABE genes when its member
    set is a different 19.
    """
    from pathlib import Path as _P
    if not per_cds_path or not _P(per_cds_path).exists():
        return {}
    out = {}
    import csv as _c
    with _P(per_cds_path).open() as fh:
        for r in _c.DictReader(fh, delimiter='\t'):
            if 'group' not in r:
                return {}          # a run predating the group column
            out.setdefault(r['family'], {}).setdefault(r['group'], set()).add(r['region'])
    return out


def _cassette_bounds(ordered):
    """(lo, hi) indices of the biosynthetic cassette within a family's gene order.

    Anchored on biosynthetic genes the family actually keeps -- role in NEVER_HIDE
    and prevalence >= 0.5. A bare prevalence-free anchor is not usable: GCF-22 has
    a stray MFS transporter at 0.10 sitting the far side of its chromosomal
    context, and anchoring on it drags the boundary over GGDEF, SiaB and SpoIIE.

    Then extended outward across adjacent UNNAMED genes at >= 0.9, which are
    unmapped rather than known-irrelevant -- GCF-18 has one at prevalence 1.00
    that the role map cannot place. Extending across any >= 0.9 gene was tried and
    overreaches: it pulled four NADH-quinone oxidoreductase subunits and EF-P
    hydroxylase into GCF-15.
    """
    anchors = [i for i, r in enumerate(ordered)
               if r['role'] in NEVER_HIDE and float(r['prevalence'] or 0) >= 0.5]
    if not anchors:
        return 0, len(ordered) - 1
    lo, hi = min(anchors), max(anchors)
    unnamed_core = lambda r: (float(r['prevalence'] or 0) >= 0.9
                              and r['consensus_product'] == '(unnamed)')
    while lo > 0 and unnamed_core(ordered[lo-1]):
        lo -= 1
    while hi < len(ordered)-1 and unnamed_core(ordered[hi+1]):
        hi += 1
    return lo, hi


def _cooccurring(rows, member_sets, min_j=0.80):
    """Group variable genes whose member sets nearly coincide, into events.

    Five rows at 0.17 that are five independent rare genes and five rows at 0.17
    that are one insertion look identical in a prevalence column, and they are not
    the same thing. On GCF-22 this resolves nine scattered rows into two events: a
    five-gene sulfur block in 19 of 186 members and the ParABE / SPI-7
    pathogenicity island in a different 32.

    Returns a list of lists, each inner list one event, largest first.
    """
    if not member_sets:
        return [[r] for r in rows]
    ids = [r['group'] for r in rows]
    parent = {i: i for i in ids}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, a in enumerate(ids):
        for b in ids[i+1:]:
            sa, sb = member_sets.get(a), member_sets.get(b)
            if not sa or not sb:
                continue
            if len(sa & sb) / len(sa | sb) >= min_j:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
    out = {}
    for r in rows:
        out.setdefault(find(r['group']), []).append(r)
    return sorted(out.values(),
                  key=lambda b: (-len(b), -float(b[0]['prevalence'] or 0)))


def _consensus_diagram(fam, fam_rows, s):
    """The consensus cluster drawn on the member that best represents it.

    A consensus cluster has no coordinates of its own: it is an abstraction over
    members sitting at different contig offsets, sometimes in different orders. So
    rather than synthesise a layout that exists nowhere, the arrows are one real
    member's -- the scaffold, chosen in gcf_annotation_transfer.py as the member
    carrying the greatest summed prevalence, which is the member that best shows
    what defines the family (1,807 of 1,816 core groups across this run, against
    1,728 when the scaffold was picked on gene count alone).

    What that costs: a group the scaffold happens not to carry has no arrow. The
    table below is still the full list, and the caption says how many of each.
    """
    placed = [r for r in fam_rows if (r.get('scaffold_start') or '').strip()]
    if not placed:
        return ''
    span = s.get('scaffold_span') or []
    try:
        lo, hi = int(span[0]), int(span[1])
    except (IndexError, ValueError, TypeError):
        lo = min(int(r['scaffold_start']) for r in placed)
        hi = max(int(r['scaffold_end']) for r in placed)

    genes = []
    for r in sorted(placed, key=lambda r: int(r['scaffold_start'])):
        role = r.get('role') or 'other'
        prev = float(r['prevalence'] or 0)
        name = r['consensus_product']
        genes.append({
            'start': int(r['scaffold_start']), 'end': int(r['scaffold_end']),
            'strand': int(r.get('scaffold_strand') or 1),
            'color': _ROLE_STYLE.get(role, _ROLE_STYLE['other'])[0],
            'locus_tag': r.get('scaffold_locus') or '',
            'product': f'{name} — in {prev:.0%} of members ({role})',
            'gene_name': '',
        })
    svg = generate_gene_svg(genes, lo, hi, width=900, height=92)

    scaffold = s.get('scaffold') or ''
    n_groups = s.get('groups') or len(fam_rows)
    missing = n_groups - len(placed)
    # Core coverage, not raw count. The scaffold is picked to maximise summed
    # prevalence, so it is typically a compact member carrying every gene that
    # defines the family and few of the accessory neighbours -- "21 of 79" reads
    # as unrepresentative when all 20 core genes are there.
    core = [r for r in fam_rows if float(r['prevalence'] or 0) >= 0.5]
    core_on = [r for r in core if (r.get('scaffold_start') or '').strip()]
    core_txt = (f', including {len(core_on)} of the {len(core)} present in at least '
                f'half the members' if core else '')
    note = (f' The remaining {missing} are accessory, and are in the table below.'
            if missing > 0 else '')
    return f'''
        <div style="margin:4px 0 14px;">
            <div style="overflow-x:auto;">{svg}</div>
            <p style="color:#777;font-size:.82em;margin:4px 0 0;max-width:80ch;">
                Drawn on <code>{_html.escape(scaffold)}</code>, the member carrying the
                most of this family&rsquo;s shared gene content &mdash;
                {len(placed)} of {n_groups} genes{core_txt}.{note}
                Arrows are coloured by role and show direction; hover for the gene
                name and how many members carry it.
            </p>
        </div>'''


def _consensus_block(fam, fam_rows, s, gcf_hrefs=None, member_sets=None):
    """One family's consensus gene table, collapsed behind a summary line.

    Split out of build_consensus_clusters_section so the per-GCF detail pages
    render the identical table rather than a second implementation of it.
    """
    # A rare gene is hidden only if its role says nothing. Anything biosynthetic,
    # transported or catabolic stays however rare -- see NEVER_HIDE.
    shown = [r for r in fam_rows
             if float(r['prevalence'] or 0) >= CONSENSUS_MIN_PREVALENCE
             or r['role'] in NEVER_HIDE]
    n_hidden = len(fam_rows) - len(shown)
    genes = _in_gene_order(shown or fam_rows)
    lo, hi = _cassette_bounds(genes)
    cassette, flank = genes[lo:hi+1], genes[:lo] + genes[hi+1:]
    core = [r for r in cassette if float(r['prevalence'] or 0) >= 0.9]
    variable = [r for r in cassette if float(r['prevalence'] or 0) < 0.9]
    events = _cooccurring(variable, (member_sets or {}).get(str(fam), {}))
    hidden_txt = (f'<span title="present in fewer than '
                  f'{CONSENSUS_MIN_PREVALENCE:.0%} of members — mostly unnamed, and '
                  f'swept in at a region boundary rather than part of the family. '
                  f'All of them are in gcf_consensus_clusters.tsv."> '
                  f'(+{n_hidden} rare)</span>') if n_hidden else ''
    # The representative cluster's own diagram and gene table live on the family
    # page; this table is the other half of the same subject, so it links across.
    href = (gcf_hrefs or {}).get(str(fam))
    page_link = (f'<a href="{href}" style="float:right;color:#2c5aa0;font-size:.85em;'
                 f'text-decoration:none;" title="Representative cluster diagram, its '
                 f'gene table, and this consensus">representative BGC &rarr;</a>'
                 if href else '')
    n_mem = s.get('members', '')
    head = f'GCF-{fam}'
    if n_mem:
        head += f' — {n_mem} member{"s" if n_mem != 1 else ""}'
    # The arrow is shown only where the transfer actually moved something. On a
    # RefSeq-only run every genome already carries PGAP annotation, so 64 of 72
    # families here move by under a point and "78.0% → 78.5%" reads as a result
    # when it is rounding. The coverage itself still matters — it says how much of
    # the family has no name anywhere — so it stays, without the arrow.
    before, after = s.get('pct_before'), s.get('pct_after')
    if before is not None and after is not None:
        if after - before >= 1.0:
            head += (f' · annotation {before}% → <strong>{after}%</strong>'
                     f' <span style="color:#166b47;" title="gained by transferring '
                     f'product names from better-annotated members of this family"'
                     f'>+{after - before:.1f}</span>')
        else:
            head += (f' · <span title="share of this family\'s CDS carrying an '
                     f'informative product. Transferring names from better-annotated '
                     f'members moved it by less than a point.">annotation '
                     f'{after}%</span>')

    def row(g, indent=False):
        prev = float(g['prevalence'] or 0)
        role = g.get('role') or 'other'
        fg, bg, _ = _ROLE_STYLE.get(role, _ROLE_STYLE['other'])
        ns = (g.get('n_sources') or '').strip()
        prov = ''
        if ns and ns.isdigit():
            n = int(ns)
            if n == 1:
                prov = ('<span title="named from a single genome — treat as one '
                        'opinion, not consensus" style="color:#9a6b0f;">1 source</span>')
            elif n > 1:
                dis = (g.get('n_disagree') or '0').strip() or '0'
                extra = f', {dis} disagreed' if dis not in ('0', '') else ''
                prov = f'<span style="color:#777;">{n} sources{extra}</span>'
        unnamed = g['consensus_product'] == '(unnamed)'
        name_html = ('<em style="color:#999;">unnamed</em>' if unnamed
                     else _html.escape(g['consensus_product']))
        # Members, not just a decimal. "6 of 38" and "0.16" are the same number and
        # read very differently when the denominator is mostly fragments.
        size, total = g.get('group_size') or '', n_mem or ''
        count = f'{size} of {total}' if size and total else f'{prev:.2f}'
        pad = 'padding:5px 9px 5px 26px;' if indent else 'padding:5px 9px;'
        return (
            f'<tr><td style="{pad}">{name_html}</td>'
            f'<td style="padding:5px 9px;white-space:nowrap;">'
            f'<span style="background:{bg};color:{fg};padding:1px 7px;'
            f'border-radius:9px;font-size:.82em;">{role}</span></td>'
            f'<td style="padding:5px 9px;font-size:.85em;color:#555;">'
            f'{_html.escape(g.get("domains") or "—")}</td>'
            f'<td style="padding:5px 9px;">'
            f'<div style="display:flex;align-items:center;gap:.4rem;">'
            f'<div style="flex:0 0 40px;height:5px;background:#e9ecef;'
            f'border-radius:3px;overflow:hidden;"><div style="width:{prev*100:.0f}%;'
            f'height:100%;background:#0e5c6b;"></div></div>'
            f'<span style="font-variant-numeric:tabular-nums;font-size:.85em;" '
            f'title="prevalence {prev:.2f}">{count}</span></div></td>'
            f'<td style="padding:5px 9px;font-size:.85em;">{prov}</td></tr>')

    def band(label, note=''):
        return (f'<tr><td colspan="5" style="padding:9px 9px 4px;background:#f3f5f6;'
                f'font-size:.82em;font-weight:600;letter-spacing:.04em;'
                f'text-transform:uppercase;color:#5a6a78;">{label}'
                f'<span style="font-weight:400;text-transform:none;letter-spacing:0;'
                f'color:#8a97a3;"> {note}</span></td></tr>')

    body = [band('Cluster — core', f'· {len(core)} genes in \u226590% of members')]
    body += [row(g) for g in core]
    if variable:
        body.append(band('Cluster — variable',
                         f'· {len(variable)} genes in {len(events)} independent event'
                         f'{"s" if len(events) != 1 else ""}'))
        for ev in events:
            if len(ev) > 1:
                ms = (member_sets or {}).get(str(fam), {})
                u = set().union(*(ms.get(g['group'], set()) for g in ev))
                body.append(
                    f'<tr><td colspan="5" style="padding:6px 9px 2px;font-size:.84em;'
                    f'color:#1d6fa5;">&#9492; co-occurring &mdash; {len(ev)} genes in the '
                    f'same {len(u)} member{"s" if len(u) != 1 else ""}, so one event</td></tr>')
                body += [row(g, indent=True) for g in ev]
            else:
                body.append(row(ev[0]))
    if flank:
        body.append(band('Neighbourhood',
                         f'· {len(flank)} genes outside the cassette'))
        body += [row(g) for g in flank]

    return (f'''
    <details style="margin:0 0 10px;border:1px solid #e3e6e8;border-radius:5px;">
        <summary style="padding:9px 13px;cursor:pointer;background:#f7f8f9;
                        border-radius:5px;">{head}
            <span style="color:#888;font-size:.9em;"> · {len(genes)} genes{hidden_txt}</span>
            {page_link}
        </summary>
        <div style="padding:4px 10px 0;">{_consensus_diagram(fam, fam_rows, s)}</div>
        <div class="table-container" style="padding:0 10px 10px;">
        <table style="width:100%;border-collapse:collapse;font-size:.9em;">
            <thead><tr style="background:#eef1f2;">
                <th style="text-align:left;padding:5px 9px;">Gene</th>
                <th style="text-align:left;padding:5px 9px;">Role</th>
                <th style="text-align:left;padding:5px 9px;">Domains</th>
                <th style="text-align:left;padding:5px 9px;">Prevalence</th>
                <th style="text-align:left;padding:5px 9px;">Naming</th>
            </tr></thead>
            <tbody>{"".join(body)}</tbody>
        </table>
        </div>
    </details>''')


def consensus_blocks_by_family(consensus_path, transfer_summary_path=None,
                               per_cds_path=None):
    """{family_id: consensus-table HTML}, for the per-GCF detail pages.

    Same renderer the report section uses, keyed by family instead of concatenated,
    so a family's consensus table is identical wherever it is read.
    """
    import csv as _csv
    from pathlib import Path as _Path
    if not consensus_path or not _Path(consensus_path).exists():
        return {}
    rows = list(_csv.DictReader(_Path(consensus_path).open(), delimiter='\t'))
    if not rows:
        return {}
    summary = {}
    if transfer_summary_path and _Path(transfer_summary_path).exists():
        try:
            summary = _json.loads(_Path(transfer_summary_path).read_text()).get('per_family', {})
        except Exception:
            summary = {}
    by_fam = {}
    for r in rows:
        by_fam.setdefault(r['family'], []).append(r)
    ms = _member_sets(per_cds_path)
    return {fam: _consensus_block(fam, fam_rows, summary.get(fam, {}), member_sets=ms)
            for fam, fam_rows in by_fam.items()}


def build_consensus_clusters_section(consensus_path, transfer_summary_path=None,
                                     gcf_hrefs=None, per_cds_path=None):
    """One consensus cluster per family, assembled from every member.

    A single representative BGC shows one genome's annotation, which on a mixed
    assembly set is usually a bad draw. Measured on Erwiniaceae, where GenBank-only
    deposits are common: 40.2% of CDS carried an informative product and 170 of 333
    regions carried none at all; pooling orthologues across a family and taking the
    majority name lifted that to 80.2%, and the consensus pantaphos cluster recovered
    a GNAT acetyltransferase and an ATP-grasp protein that LMG 5342's own annotation
    calls "hypothetical". The gain is much smaller on a RefSeq-only set — 87.8% to
    88.2% on Enterobacterales — because PGAP has already annotated every genome. The
    rendered section reports whichever applies to the run in hand.

    Prevalence is what makes this more than a longer gene list: a gene at 1.00 is in
    every member and is part of what defines the family; one at 0.24 is accessory and
    may be a neighbouring gene the region boundary caught. Roles come from Pfam
    accessions, not product text, so they are computed identically whether or not
    NCBI annotated the assembly.
    """
    import csv as _csv
    from pathlib import Path as _Path
    if not consensus_path or not _Path(consensus_path).exists():
        return ''
    rows = list(_csv.DictReader(_Path(consensus_path).open(), delimiter='\t'))
    if not rows:
        return ''

    summary = {}
    if transfer_summary_path and _Path(transfer_summary_path).exists():
        try:
            summary = _json.loads(_Path(transfer_summary_path).read_text()).get('per_family', {})
        except Exception:
            summary = {}

    by_fam = {}
    for r in rows:
        by_fam.setdefault(r['family'], []).append(r)
    ms = _member_sets(per_cds_path)

    def order(fam):
        return -int(summary.get(fam, {}).get('members', 0) or 0), int(fam)

    blocks = []
    for fam in sorted(by_fam, key=order):
        blocks.append(_consensus_block(fam, by_fam[fam], summary.get(fam, {}),
                                       gcf_hrefs, member_sets=ms))

    # Coverage from this run, not from the run this section was written against.
    # The prose used to quote 40.2% -> 80.2% over 333 regions as though it were a
    # property of the method; those are Erwiniaceae numbers, and on a RefSeq-only
    # order the same step moves 87.8% -> 88.2%. Quoting them as "this run" made the
    # section describe data the reader was not looking at.
    overall = {}
    if transfer_summary_path and _Path(transfer_summary_path).exists():
        try:
            overall = _json.loads(_Path(transfer_summary_path).read_text())
        except Exception:
            overall = {}
    b, a = overall.get('pct_before'), overall.get('pct_after')
    n_cds = overall.get('cds')
    if b is not None and a is not None and n_cds:
        gain = a - b
        measured = (
            f'Across the {n_cds:,} CDS in these families that took annotation '
            f'coverage from <strong>{b}%</strong> to <strong>{a}%</strong>.'
            if gain >= 1.0 else
            f'These families were already <strong>{b}%</strong> annotated across '
            f'{n_cds:,} CDS, so the transfer had little left to do here — it earns '
            f'its place on assembly sets that include GenBank-only deposits with no '
            f'functional annotation at all, where it has lifted coverage from 40.2% '
            f'to 80.2%.')
    else:
        measured = ('A product name found on any member is propagated to the '
                    'orthologues that lack one.')

    legend = ' '.join(
        f'<span style="background:{bg};color:{fg};padding:1px 7px;border-radius:9px;'
        f'font-size:.82em;margin-right:6px;" title="{tip}">{role}</span>'
        for role, (fg, bg, tip) in _ROLE_STYLE.items() if role != 'other')

    return f'''
    <div class="section">
        <h3>Consensus Gene Content</h3>
        <p style="color:#555;max-width:72ch;">
            One consensus cluster per family, assembled from <strong>every member</strong>
            rather than a single representative. Annotation quality is distributed very
            unevenly within a family, so any one representative is a lottery; pooling
            orthologues across the family and taking the majority name removes that draw.
            {measured}
        </p>
        <p style="color:#555;max-width:72ch;font-size:.92em;">
            <strong>Prevalence is the column to read.</strong> A gene present in every
            member (1.00) is part of what defines the family; one at 0.24 is accessory, and
            may simply be a neighbour the region boundary caught. <strong>Roles come from
            Pfam accessions, not product text</strong>, so they are computed identically
            whether or not NCBI annotated the assembly — which is why
            <em>serine hydroxymethyltransferase</em> lands in <code>primary metabolism</code>
            rather than being counted as a tailoring methyltransferase.
        </p>
        <p style="color:#555;max-width:72ch;font-size:.92em;">
            A transferred name is an <strong>inference from a homologue, not an
            observation</strong>. The Naming column says how many independent genomes
            supported it; a single source is one genome's opinion.
        </p>
        <p style="margin:12px 0 14px;">{legend}</p>
        {"".join(blocks)}
    </div>'''
