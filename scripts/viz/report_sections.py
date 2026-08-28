"""HTML section builders for the BGC report.

Each function returns a self-contained block of markup that generate_html_report()
drops into the page, keeping the report generator readable. The _build_* helpers are
pure string building; build_coupling_table_rows additionally reads the BiG-SCAPE DB so
the coupling table reflects the current run rather than hardcoded family IDs.
"""

import sqlite3

from utils.constants import load_coupling_classes, COUPLING_COLORS


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
                                    f'padding: 2px 8px; border-radius: 4px; font-size: 0.85em;">GCF-{fid}</span>'
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


def _build_bigscape_overview_cards(gcf_data):
    """Return overview-grid HTML cards for BiG-SCAPE GCF summary stats."""
    if not (gcf_data and isinstance(gcf_data, dict) and gcf_data.get('gcfs')):
        return ''
    gcfs = gcf_data['gcfs']
    gcf_summary = gcf_data.get('summary', {})
    total_gcfs = gcf_summary.get('total', len(gcfs))
    singletons = gcf_summary.get('singletons', sum(1 for g in gcfs if g.get('is_singleton')))
    clusters = total_gcfs - singletons
    largest_gcf = max(gcfs, key=lambda x: x.get('member_count', 0)) if gcfs else None
    largest_gcf_card = ''
    if largest_gcf:
        largest_gcf_card = f'''
                <div style="background: rgba(122, 104, 85, 0.15); border: 1px solid rgba(122, 104, 85, 0.3); padding: 10px 14px; border-radius: 8px;">
                    <div style="font-size: 1.5em; font-weight: bold; color: #7a6855;">{largest_gcf.get("member_count", 0)}</div>
                    <div style="color: #555; font-size: 0.85em;">Largest GCF Size</div>
                    <div style="font-size: 0.78em; margin-top: 2px; color: #777;">{largest_gcf.get("product", "")[:30]}</div>
                </div>'''
    return f'''
                <div style="background: rgba(92, 107, 122, 0.15); border: 1px solid rgba(92, 107, 122, 0.3); padding: 10px 14px; border-radius: 8px;">
                    <div style="font-size: 1.5em; font-weight: bold; color: #5c6b7a;">{total_gcfs}</div>
                    <div style="color: #555; font-size: 0.85em;">Gene Cluster Families</div>
                </div>
                <div style="background: rgba(90, 122, 107, 0.15); border: 1px solid rgba(90, 122, 107, 0.3); padding: 10px 14px; border-radius: 8px;">
                    <div style="font-size: 1.5em; font-weight: bold; color: #5a7a6b;">{clusters}</div>
                    <div style="color: #555; font-size: 0.85em;">Multi-member GCFs</div>
                </div>
                <div style="background: rgba(122, 122, 122, 0.15); border: 1px solid rgba(122, 122, 122, 0.3); padding: 10px 14px; border-radius: 8px;">
                    <div style="font-size: 1.5em; font-weight: bold; color: #7a7a7a;">{singletons}</div>
                    <div style="color: #555; font-size: 0.85em;">Singletons</div>
                </div>
                {largest_gcf_card}'''


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
                variant unlike the one characterised example — both warrant a look. Classes
                with few references (see <em>refs</em>) give weaker evidence either way.
            </p>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>GCF</th><th>Dominant coupling class</th><th>Members</th>
                            <th>Median support</th><th>Range</th><th>Refs</th><th></th>
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
                <h3>BiG-SCAPE Gene Cluster Families</h3>
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
_COUPLING_META = {
    'Synthase':                        ('Synthase',                        'SMCOG1271 (HMGL-like)',          '→ phosphonomethylmalate → phosphinothricin-type', 'FrbC, HvrC'),
    'Reductase':                       ('Reductase',                       'Fe-ADH rule',                    '→ phosphonolactate (reductase route)',            'VlpB'),
    'Decarboxylase-Nucleotidyltransferase': ('Decarboxylase-Nucleotidyltransferase', 'SMCOG1055 + NTP_transf_3', '→ phosphonolipid (CDP-pathway)',             'DhpF, Fom2, Ppd'),
    'Decarboxylase':                   ('Decarboxylase',                   'SMCOG1055 (ThDP-dependent)',      '→ 2-phosphonoacetaldehyde → 2-AEP',              'DhpF, Fom2, Ppd'),
    'Transaminase':                    ('Transaminase',                    'SMCOG1019 (Aminotran_1_2/PF00155)', '→ L-phosphonoalanine',                         'PnaA, PalB'),
    'Unknown':                         ('Unknown',                         '—',                               '—',                                              '—'),
}
_COUPLING_ROW_ORDER = ['Synthase', 'Reductase', 'Decarboxylase-Nucleotidyltransferase', 'Decarboxylase', 'Transaminase', 'Unknown']

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
        for row in csv.DictReader(lines, delimiter='\t'):
            fid = bgc_to_gcf.get(row['bgc'])
            if fid is None:
                continue
            try:
                per_gcf[fid].append(float(row['assigned_pct_id']))
            except ValueError:
                continue
            n_refs[fid] = row.get('assigned_n_refs', '?')
        if not per_gcf:
            return None

        rows_html = []
        for i, fid in enumerate(sorted(per_gcf, key=lambda k: -len(per_gcf[k]))):
            vals = per_gcf[fid]
            cls = gcf_class.get(fid, 'Unknown')
            med = statistics.median(vals)
            colour = COUPLING_COLORS.get(cls, '#999999')
            # Support is advisory: a low value may mean the wrong class, or a novel
            # variant unlike the single characterised example. Flagged, never filtered.
            note = ('well evidenced' if med >= 90 else
                    'moderate' if med >= 50 else
                    'weak — review')
            bg = ' style="background:#fafafa;"' if i % 2 else ''
            td = 'padding: 7px 12px; border-bottom: 1px solid #eee;'
            rows_html.append(
                f'<tr{bg}>'
                f'<td style="{td}"><span style="background: {colour}; color: white; '
                f'padding: 2px 8px; border-radius: 4px; font-size: 0.85em;">GCF-{fid}</span></td>'
                f'<td style="{td}">{cls}</td>'
                f'<td style="{td} text-align: center;">{len(vals)}</td>'
                f'<td style="{td} text-align: center;">{med:.1f}%</td>'
                f'<td style="{td} text-align: center;">{min(vals):.1f}–{max(vals):.1f}%</td>'
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
            gcf_label = '/'.join(f'GCF-{g}' for g in gcf_ids) if gcf_ids else '—'
            display, marker, pathway, refs = _COUPLING_META[cls_id]
            bg = ' style="background:#fafafa;"' if i % 2 == 1 else ''
            is_last = (i == len(_COUPLING_ROW_ORDER) - 1)
            border = '' if is_last else 'border-bottom: 1px solid #eee; '
            td = f'padding: 7px 12px; {border}'
            rows_html.append(
                f'<tr{bg}>'
                f'<td style="{td}">{display}</td>'
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
