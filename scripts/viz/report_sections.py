"""HTML section builders for the BGC report.

Each function returns a self-contained block of markup that generate_html_report()
drops into the page, keeping the report generator readable. The _build_* helpers are
pure string building; build_coupling_table_rows additionally reads the BiG-SCAPE DB so
the coupling table reflects the current run rather than hardcoded family IDs.
"""

import sqlite3

from utils.constants import load_coupling_classes, COUPLING_COLORS
from utils.coupling_confidence import BACKGROUND_CEILING_PCT


def _build_kcb_content(kcb_stats, taxon_clean, gcf_data, gcf_classes=None):
    """Compute KCB tab contents.

    Returns dict with keys: kcb_mapping_section, novel_bgcs_tab_content,
    kcb_hits_tab_content.
    """
    kcb_mapping_section = ''''''
    novel_bgcs_tab_content = '''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666;">No KnownClusterBlast data available. Run antiSMASH with <code>--antismash_cb_knownclusters true</code> to identify potentially novel BGCs.</p>'''
    kcb_hits_tab_content = '''
            <h2>KnownClusterBlast Hits</h2>
            <p style="color: #666;">No KnownClusterBlast data available. Run antiSMASH with <code>--antismash_cb_knownclusters true</code> to identify known cluster matches.</p>'''

    if kcb_stats.get('total_regions', 0) > 0:
        novel_bgcs = kcb_stats.get('novel_bgcs', [])
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
                        gcf_cell = (f'<td style="text-align: center;">'
                                    f'<span{badge_title} style="background: {badge_bg}; color: white; '
                                    f'padding: 2px 8px; border-radius: 4px; font-size: 0.85em; '
                                    f'display: inline-block; white-space: nowrap;">GCF-{fid}</span>'
                                    f'</td><td style="text-align: center;">{mc}</td>')
                    else:
                        gcf_cell = '<td style="text-align: center; color: #999;">-</td><td style="text-align: center; color: #999;">-</td>'
                detail_rows += f'''
                <tr>
                    <td><a href="genomes/{genome}.html" title="{genome}" style="color: #2c5aa0; display: inline-block; max-width: 340px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; vertical-align: middle;">{genome}</a></td>
                    <td style="text-align: center;"><a href="{antismash_link}" target="_blank" style="color: #28a745; font-weight: bold;">Region {region_name}</a></td>
                    <td>{product}</td>
                    <td style="text-align: center;">{edge_badge}</td>
                    {gcf_cell}
                </tr>'''

            kcb_mapping_section = ''''''
            gcf_header = '<th>GCF Family</th><th>Members</th>' if has_gcf_data else ''
            gcf_description = ' When BiG-SCAPE clustering is enabled, the GCF (Gene Cluster Family) assignment shows how these novel BGCs group together.' if has_gcf_data else ''
            novel_bgcs_tab_content = f'''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666; margin-bottom: 15px;">
                <em>These BGC regions did not return any hits from KnownClusterBlast (KCB) analysis against the MIBiG database,
                suggesting they may encode novel or uncharacterized biosynthetic pathways. Regions marked as "edge" are on contig
                boundaries and may be incomplete.{gcf_description}</em>
            </p>
            <div class="search-box">
                <input type="text" id="novelSearch" placeholder="Search by genome, strain, region or GCF (e.g. 5342)" onkeyup="filterNovelBGCs()">
            </div>
            <div class="table-container">
                <table id="novelTable">
                    <thead>
                        <tr>
                            <th>Genome</th>
                            <th>antiSMASH Region</th>
                            <th>Product Type</th>
                            <th>Contig Edge</th>
                            {gcf_header}
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

        # Build KCB Hits tab content
        cluster_mapping = kcb_stats.get('cluster_mapping', [])
        unique_clusters = kcb_stats.get('unique_known_clusters', 0)
        sim_breakdown = kcb_stats.get('similarity_breakdown', {})
        sim_badges = ''
        for sim_level, count in sim_breakdown.items():
            color = '#27ae60' if sim_level == 'high' else '#f39c12' if sim_level == 'medium' else '#95a5a6'
            sim_badges += f'<span style="background: {color}; color: white; padding: 4px 12px; border-radius: 4px; margin-right: 8px;">{sim_level}: {count}</span>'
        kcb_table_rows = ''
        for item in cluster_mapping:
            known_cluster = item['known_cluster']
            mibig_acc = item['mibig_acc']
            hit_count = item['count']
            regions = item['regions']
            product_types = list(set(r['product'] for r in regions))
            products_display = ', '.join(product_types[:3])
            if len(product_types) > 3:
                products_display += f' (+{len(product_types) - 3} more)'
            bgc_links = []
            for r in regions[:5]:
                genome = r['genome']; region = r['region']
                region_name = r.get('region_name', region); record_index = r.get('record_index', 1)
                antismash_link = f'../../antismash_results/{taxon_clean}/{genome}/index.html#r{record_index}c{region}'
                bgc_links.append(f'<a href="{antismash_link}" target="_blank" style="color: #28a745;">{genome} Region {region_name}</a>')
            bgc_display = ', '.join(bgc_links)
            if len(regions) > 5:
                bgc_display += f' (+{len(regions) - 5} more)'
            sim_counts = {}
            for r in regions:
                sim = r.get('similarity', 'unknown'); sim_counts[sim] = sim_counts.get(sim, 0) + 1
            sim_display = ' / '.join([f'{k}: {v}' for k, v in sim_counts.items()])
            if 'high' in sim_counts:
                row_bg = 'background: rgba(39, 174, 96, 0.15);'
            elif 'medium' in sim_counts:
                row_bg = 'background: rgba(243, 156, 18, 0.15);'
            elif 'low' in sim_counts:
                row_bg = 'background: rgba(149, 165, 166, 0.15);'
            else:
                row_bg = ''
            mibig_link = f'https://mibig.secondarymetabolites.org/repository/{mibig_acc}/' if mibig_acc else '#'
            kcb_table_rows += f'''
                <tr style="{row_bg}">
                    <td><a href="{mibig_link}" target="_blank" style="color: #2c5aa0; font-weight: bold;">{known_cluster}</a></td>
                    <td><a href="{mibig_link}" target="_blank" style="color: #666;">{mibig_acc}</a></td>
                    <td style="text-align: center; font-weight: bold;">{hit_count}</td>
                    <td>{products_display}</td>
                    <td>{sim_display}</td>
                    <td>{bgc_display}</td>
                </tr>'''
        if not cluster_mapping:
            # Zero hits is a result, not a failure: it means every region is
            # potentially novel. Rendering bare table headers reads like a bug.
            total_regions = kcb_stats.get('total_regions', 0)
            kcb_hits_tab_content = f'''
            <h2>KnownClusterBlast Hits</h2>
            <div style="background: #eef6ec; border: 1px solid #cfe3ca; border-radius: 6px; padding: 18px 20px; margin-top: 10px;">
                <strong>No KnownClusterBlast hits.</strong>
                <p style="margin: 8px 0 0; color: #555;">
                    None of the {total_regions} detected region{'s' if total_regions != 1 else ''}
                    matched a characterised cluster in the MIBiG database, so all of them appear in
                    the <em>Novel BGCs</em> tab. For phosphonate BGCs this is common — MIBiG holds
                    relatively few characterised phosphonate pathways — and it is a finding rather
                    than an error.
                </p>
            </div>'''
        else:
            kcb_hits_tab_content = f'''
            <h2>KnownClusterBlast Hits</h2>
            <div style="margin-bottom: 20px;">
                <strong>Total Hits by Similarity:</strong> {sim_badges}
            </div>
            <div style="background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 6px; padding: 15px; margin-bottom: 20px;">
                <strong>Similarity Legend:</strong>
                <div style="margin-top: 10px; display: flex; gap: 20px; flex-wrap: wrap;">
                    <div><span style="background: #27ae60; color: white; padding: 2px 8px; border-radius: 3px;">high</span> &gt;75% sequence similarity to MIBiG reference</div>
                    <div><span style="background: #f39c12; color: white; padding: 2px 8px; border-radius: 3px;">medium</span> 50-75% sequence similarity</div>
                    <div><span style="background: #95a5a6; color: white; padding: 2px 8px; border-radius: 3px;">low</span> 15-50% sequence similarity</div>
                </div>
                <p style="margin: 10px 0 0 0; font-size: 0.9em; color: #666;">
                    The "Similarity" column shows how many regions matched each known cluster at different similarity levels (e.g., "high: 1 / low: 2" means 1 region matched with &gt;75% similarity and 2 regions matched with 15-50% similarity).
                </p>
            </div>
            <div class="search-box">
                <input type="text" id="kcbSearch" placeholder="Search known clusters..." onkeyup="filterKCBHits()">
            </div>
            <div class="table-container">
                <table id="kcbTable">
                    <thead>
                        <tr>
                            <th>Known Cluster</th>
                            <th>MIBiG ID</th>
                            <th>Hits</th>
                            <th>Product Types</th>
                            <th>Similarity</th>
                            <th>BGC Regions</th>
                        </tr>
                    </thead>
                    <tbody id="kcbTableBody">
                        {kcb_table_rows}
                    </tbody>
                </table>
            </div>
            <p style="color: #666; font-style: italic; margin-top: 15px;">
                Table shows {unique_clusters} unique known clusters from MIBiG database. Each row represents a characterized BGC that matched one or more regions in your dataset.
            </p>'''

    return {
        'kcb_mapping_section':    kcb_mapping_section,
        'novel_bgcs_tab_content': novel_bgcs_tab_content,
        'kcb_hits_tab_content':   kcb_hits_tab_content,
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
            <p style="color: #666; font-size: 0.9em; margin-bottom: 12px;">
                GCF membership comes from BiG-SCAPE, which compares the whole gene
                neighbourhood. The coupling class comes from antiSMASH SMCOG/domain markers,
                which are deliberately broad — characterised phosphonate coupling enzymes are
                scarce, and a narrow reference-driven classifier would only recover chemistry
                already known. <em>Support</em> is the identity of the enzyme that drove each
                call to the nearest characterised reference of its class. It is advisory:
                a low value may mean the assignment is wrong, or that the enzyme is a novel
                variant unlike the one characterised example — both warrant a look.
                <strong>Read it against the reference it was scored on.</strong> Most
                characterised phosphonate enzymes come from <em>Streptomyces</em>, so a
                modest identity may simply reflect the genus gap rather than a doubtful
                call; the one class with a <em>Pantoea</em> reference (Synthase, HvrC)
                scores near 100% partly for that reason. Only one boundary is measurable
                here: characterised enzymes of <em>different</em> classes score 26.7–29.7%
                against each other, so at or below ~30% an identity carries no class
                information. Above that the reference set is too small, and drawn from too
                few genera, to support a verdict — so none is given.
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


def _build_bigscape_section_html(bigscape_stats_html, gcf_visualization_html, taxon_clean,
                                 gcf_support_rows=None):
    """Return the BiG-SCAPE GCF section HTML for the GCF Analysis tab."""
    if not (bigscape_stats_html or gcf_visualization_html):
        return ''
    return f'''
            <div class="clustering-section">
                {bigscape_stats_html if bigscape_stats_html else ''}
                {_build_gcf_support_section(gcf_support_rows)}
                {gcf_visualization_html if gcf_visualization_html else ''}
                <div class="info-box" style="margin-top: 20px;">
                    <p>BiG-SCAPE clusters biosynthetic gene clusters into gene cluster families based on sequence similarity.</p>
                    <p style="margin-top: 15px;"><strong>To view interactive results:</strong></p>
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
        'Transaminase', 'SMCOG1019 (Aminotran_1_2 / PF00155)',
        'L-Phosphonoalanine', 'PnaA, PalB'),
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


def build_pepm_section(pepm_b64, pepm_summary):
    """pepM identity vs gene-neighbourhood similarity, for the GCF Analysis tab.

    This is the evidence that the GCF assignments above it can be trusted: it
    reproduces Yu et al. (PNAS 2013;110(51):20759) Fig. 2B on the run's own data,
    correlating each BGC pair's pepM amino-acid identity against the domain
    content BiG-SCAPE actually clustered on.

    Returns '' when the analysis did not run, so the tab simply omits it.
    """
    if not pepm_summary:
        return ''
    stats = (pepm_summary.get('bigscape_similarity') or {}).get('regression') or {}
    jac = (pepm_summary.get('jaccard') or {}).get('regression') or {}
    n_seq = pepm_summary.get('sequences')
    n_pairs = pepm_summary.get('pairs')

    fig = ''
    if pepm_b64:
        fig = (f'<img src="data:image/svg+xml;base64,{pepm_b64}" '
               f'alt="pepM identity against BiG-SCAPE gene-cluster similarity, '
               f'{n_pairs:,} pairwise comparisons" '
               f'style="max-width:100%;height:auto;display:block;margin:0 auto;">')

    def row(label, r):
        if not r:
            return ''
        return (f'<tr><td style="padding:6px 10px;">{label}</td>'
                f'<td style="padding:6px 10px;text-align:right;">{r.get("r", 0):+.3f}</td>'
                f'<td style="padding:6px 10px;text-align:right;">{r.get("r2", 0):.3f}</td>'
                f'<td style="padding:6px 10px;text-align:right;">{r.get("slope", 0):+.2f}</td>'
                f'<td style="padding:6px 10px;text-align:right;">{r.get("n", 0):,}</td></tr>')

    return f'''
    <div class="section">
        <h3>pepM Identity vs Gene-Cluster Similarity</h3>
        <p style="color:#555;max-width:70ch;">
            Every pair of pepM (PEP mutase) sequences in this run compared against each
            other, plotted against how similar BiG-SCAPE found their gene neighbourhoods
            &mdash; a replication of Yu et&nbsp;al.
            (<a href="https://doi.org/10.1073/pnas.1315107110">PNAS 2013;110(51):20759</a>)
            Fig.&nbsp;2B on this dataset. Identity is computed from one alignment with
            pairwise deletion of missing sites, as in the paper, not from BLAST.
            <strong>{n_seq} pepM sequences, {n_pairs:,} pairwise comparisons.</strong>
        </p>
        {fig}
        <table style="border-collapse:collapse;margin-top:14px;font-size:0.9em;">
            <thead><tr style="background:#f6f7f8;">
                <th style="padding:6px 10px;text-align:left;">Neighbourhood measure</th>
                <th style="padding:6px 10px;">r</th><th style="padding:6px 10px;">r&sup2;</th>
                <th style="padding:6px 10px;">slope</th><th style="padding:6px 10px;">pairs</th>
            </tr></thead>
            <tbody>
                {row('BiG-SCAPE similarity (1 &minus; distance)', stats)}
                {row('Shared domain content (Jaccard)', jac)}
            </tbody>
        </table>
        <p style="color:#777;font-size:0.85em;margin-top:10px;max-width:70ch;">
            Fitted over the paper's 0.6&ndash;1.0 identity window. The correlation is
            expected to be strong in a taxonomically diverse set and weak in one dominated
            by a single closely-related family, where nearly every pepM pair is either
            near-identical or unrelated with little in between.
        </p>
    </div>
    '''


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


# ─── Tab bodies ────────────────────────────────────────────────────────────────
# These were inline in a single 300-line f-string inside generate_html_report,
# which made every tab edit a careful string match into a wall of markup. Each
# returns the inner HTML of one `.tab-content` div; the surrounding div and the
# tab nav stay in visualize_results.py, where the tab numbering lives.

_MISSING = ('<div style="color: #999; padding: 20px; background: #f8f9fa; '
            'border-radius: 8px; text-align: center; font-size: 0.9em;">{}</div>')


def build_gcf_analysis_tab(coupling_table_rows, bigscape_section_html, pepm_section_html,
                           priority_html=''):
    """Clustering statistics, the coupling-enzyme class table and the pepM figure.

    The trees themselves are in the GCF Trees tab; the class table stays here
    because it is a classification reference, and the tree figures carry their
    own colour legends.
    """
    no_clustering = ('<div class="info-box" style="background-color: #f8f9fa; '
                     'border-left: 4px solid #6c757d;"><p style="color: #666;">'
                     'No clustering analysis was performed. To enable clustering, run the '
                     'pipeline with <code>--clustering bigscape</code>.</p></div>')
    th = ('text-align: left; padding: 8px 12px; border-bottom: 2px solid #dee2e6;')
    return f'''
            {priority_html}

            <h3>GCF Biosynthetic Phylogeny</h3>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Classification of phosphonate BGCs by the coupling enzyme acting on phosphonopyruvate — the branching step
                immediately downstream of PEP mutase that determines the downstream biosynthetic pathway.
                The trees themselves are in the <strong>GCF Trees</strong> tab.</em>
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
            </div>

            <hr class="tab-section-divider">

            {bigscape_section_html}
            {pepm_section_html}
            {'' if bigscape_section_html else no_clustering}'''


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
            <h3>Gene Cluster Family Trees</h3>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Branch colours are coupling enzyme classes; the class table is in the
                <strong>GCF Analysis</strong> tab.</em>
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


def build_priority_section(ranking_path):
    """Gene cluster families ranked by how much they warrant laboratory follow-up.

    Distance and evidence are shown beside the priority they multiply to, deliberately.
    The weights behind them are reasoned, not fitted — there is no set of leads that
    panned out to fit against — and a single number would hide that. With the components
    visible a reader can disagree with the weighting and re-order by eye: the score
    orders the list, the components justify the order.

    Families whose coupling enzyme matches no characterised class are listed first and
    carry no score at all. An enzyme resembling nothing described is the strongest
    novelty signal in the data, and it is the one thing a distance cannot express.
    """
    import csv as _csv
    from pathlib import Path as _Path
    if not ranking_path or not _Path(ranking_path).exists():
        return ''
    rows = list(_csv.DictReader(_Path(ranking_path).open(), delimiter='\t'))
    if not rows:
        return ''

    unc = [r for r in rows if r['status'] == 'unclassified']
    ranked = [r for r in rows if r['status'] == 'ranked']

    def bar(frac, tone):
        pct = max(0.0, min(1.0, frac)) * 100
        return (f'<div style="display:flex;align-items:center;gap:.45rem">'
                f'<div style="flex:0 0 46px;height:6px;background:#e9ecef;border-radius:3px;'
                f'overflow:hidden"><div style="width:{pct:.0f}%;height:100%;'
                f'background:{tone}"></div></div>'
                f'<span style="font-variant-numeric:tabular-nums">{frac:.2f}</span></div>')

    unc_html = ''
    if unc:
        items = ''.join(
            f'<tr><td style="padding:6px 10px;white-space:nowrap;">'
            f'<a href="javascript:void(0)" onclick="showGCF(\'{r["gcf"]}\')" '
            f'style="color:#8a5a0c;font-weight:600;text-decoration:none;'
            f'border-bottom:1px dotted #8a5a0c;">GCF-{r["gcf"]}</a></td>'
            f'<td style="padding:6px 10px;text-align:right;">{r["members"]}</td>'
            f'<td style="padding:6px 10px;text-align:right;">{r["genomes"]}</td>'
            f'<td style="padding:6px 10px;text-align:right;">{float(r["intact"]):.0%}</td></tr>'
            for r in unc)
        unc_html = f'''
        <div style="background:#fdf6e3;border-left:3px solid #9a6b0f;padding:14px 18px;margin:0 0 22px;">
            <h4 style="margin:0 0 6px;">Unclassifiable coupling enzyme — {len(unc)} famil{'y' if len(unc)==1 else 'ies'}</h4>
            <p style="margin:0 0 10px;color:#555;font-size:.9em;max-width:66ch;">
                These matched no characterised coupling class, so they carry no distance and
                are not ranked. That is either a truncated cluster or chemistry with no
                described analogue — worth a look by hand before anything below.
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

    body = ''.join(
        f'<tr>'
        f'<td style="padding:7px 10px;color:#888;text-align:right;">{r["rank"]}</td>'
        f'<td style="padding:7px 10px;white-space:nowrap;">'
        f'<a href="javascript:void(0)" onclick="showGCF(\'{r["gcf"]}\')" '
        f'style="color:#2c5aa0;font-weight:600;text-decoration:none;'
        f'border-bottom:1px dotted #2c5aa0;" '
        f'title="Show this family in Gene Cluster Families">GCF-{r["gcf"]}</a></td>'
        f'<td style="padding:7px 10px;font-weight:600;text-align:right;'
        f'font-variant-numeric:tabular-nums;">{float(r["priority"]):.3f}</td>'
        f'<td style="padding:7px 10px;">{bar(float(r["distance"]), "#0e5c6b")}</td>'
        f'<td style="padding:7px 10px;">{bar(float(r["evidence"]), "#166b47")}</td>'
        f'<td style="padding:7px 10px;text-align:right;">{r["members"]}</td>'
        f'<td style="padding:7px 10px;text-align:right;">{r["genera"]}</td>'
        f'<td style="padding:7px 10px;text-align:right;">{float(r["intact"]):.0%}</td>'
        f'<td style="padding:7px 10px;">{r["coupling_class"]}</td>'
        f'</tr>'
        for r in ranked)

    return f'''
    <div class="section">
        <h3>Priority for Laboratory Follow-Up</h3>
        <p style="color:#555;max-width:70ch;">
            Families ordered by <strong>distance × evidence</strong>. Distance is how far the
            coupling enzyme sits from any characterised one; evidence is how confident we can
            be the family is real rather than an assembly artefact — independent genomes,
            independent genera, and the share of regions not truncated at a contig edge.
            They multiply because both are necessary.
        </p>
        <p style="color:#555;max-width:70ch;font-size:.92em;">
            <strong>The components are shown deliberately.</strong> Their weights are reasoned,
            not fitted to any set of leads that panned out, so the ordering is a considered
            opinion rather than a measurement. Read across the row, not just down the score:
            a family with high distance and low evidence is a different proposition from a
            middling one on both.
        </p>
        {unc_html}
        <div class="table-container">
        <table style="width:100%;border-collapse:collapse;font-size:.9em;">
            <thead><tr style="background:#e9ecef;">
                <th style="text-align:right;padding:6px 10px;">#</th>
                <th style="text-align:left;padding:6px 10px;">Family</th>
                <th style="text-align:right;padding:6px 10px;">Priority</th>
                <th style="text-align:left;padding:6px 10px;">Distance</th>
                <th style="text-align:left;padding:6px 10px;">Evidence</th>
                <th style="text-align:right;padding:6px 10px;">BGCs</th>
                <th style="text-align:right;padding:6px 10px;">Genera</th>
                <th style="text-align:right;padding:6px 10px;">Intact</th>
                <th style="text-align:left;padding:6px 10px;">Coupling class</th>
            </tr></thead>
            <tbody>{body}</tbody>
        </table>
        </div>
    </div>'''


def build_novelty_tab(priority_html, all_regions_html, n_regions=0):
    """BGC Novelty: what to look at, then everything else.

    The ranking is 17 rows; the full region list is 333. Presenting them as equals
    buries the actionable part under a table where every row says the same thing —
    which is what the old "Novel BGCs" tab did, at 30.7% of the whole report. The list
    is kept, because it is genuinely useful, but folded into a `<details>` so it is one
    click away rather than the first thing in the section.
    """
    count = f' ({n_regions:,} regions)' if n_regions else ''
    listing = f'''
        <details style="margin-top:26px;border:1px solid #dee2e6;border-radius:8px;padding:0;">
            <summary style="cursor:pointer;padding:13px 18px;font-weight:600;background:#f8f9fa;
                            border-radius:8px;">
                All detected regions{count}
            </summary>
            <div style="padding:4px 18px 18px;">
                <p style="color:#666;font-size:.9em;max-width:70ch;">
                    Every phosphonate region found, whether or not its family ranked above.
                    Use this to locate a specific contig or genome; use the ranking to decide
                    what to work on.
                </p>
                {all_regions_html}
            </div>
        </details>''' if all_regions_html else ''

    if not priority_html and not all_regions_html:
        return ('<div class="info-box"><p style="color:#666;">No BGC novelty analysis '
                'available. Run with <code>--clustering bigscape</code> to enable.</p></div>')

    return f'''
            <h2>BGC Novelty</h2>
            <p style="color:#666;max-width:70ch;">
                <em>Which gene cluster families are worth taking into the laboratory, and why.
                Ordered by divergence from characterised chemistry, discounted by how well
                evidenced each family is.</em>
            </p>
            {priority_html}
            {listing}'''
