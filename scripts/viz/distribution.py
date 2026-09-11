"""GCF distribution across taxonomy: heatmap, genus-specific and widespread tables."""

import json
import re
from collections import defaultdict

import pandas as pd


def extract_assembly_id_from_genome_name(genome_name):
    """Extract GCA/GCF assembly ID from genome name."""
    # Look for GCA_XXXXXXXXX.X or GCF_XXXXXXXXX.X pattern
    match = re.search(r'(GC[AF]_\d+\.\d+)', genome_name)
    if match:
        return match.group(1)
    return None

def build_gcf_taxonomy_distribution(gcf_data, taxonomy_map):
    """
    Build GCF × Taxonomy distribution data.

    Returns:
        dict with:
        - 'heatmap_data': {gcf_id: {genus: count, ...}, ...}
        - 'gcf_taxonomy': {gcf_id: {'genera': {genus: count}, 'species': {species: count}}, ...}
        - 'genera': sorted list of all genera
        - 'species': sorted list of all species
        - 'top_gcfs': list of top GCFs by member count
    """
    if not gcf_data or not taxonomy_map:
        return None

    genome_gcf_mapping = gcf_data.get('genome_gcf_mapping', {})
    family_metadata = gcf_data.get('family_metadata', {})

    if not genome_gcf_mapping:
        return None

    # Build genome -> taxonomy lookup
    # Try genome_name directly first (works for GTDB-Tk keyed by genome name),
    # then fall back to assembly_id extraction (works for NCBI taxonomy map).
    genome_to_taxonomy = {}
    for genome_name in genome_gcf_mapping.keys():
        tax_info = None
        if genome_name in taxonomy_map:
            tax_info = taxonomy_map[genome_name]
        else:
            assembly_id = extract_assembly_id_from_genome_name(genome_name)
            if assembly_id and assembly_id in taxonomy_map:
                tax_info = taxonomy_map[assembly_id]
        if tax_info:
            lineage = tax_info.get('lineage', {})
            genus = lineage.get('genus', {}).get('name', 'Unknown')
            species = lineage.get('species', {}).get('name', 'Unknown')
            genome_to_taxonomy[genome_name] = {
                'genus': genus,
                'species': species
            }

    # Build GCF -> taxonomy distribution
    gcf_taxonomy = {}  # {gcf_id: {'genera': {genus: count}, 'species': {species: count}}}
    all_genera = set()
    all_species = set()

    for genome_name, gcf_ids in genome_gcf_mapping.items():
        tax = genome_to_taxonomy.get(genome_name)
        if not tax:
            continue

        genus = tax['genus']
        species = tax['species']
        all_genera.add(genus)
        all_species.add(species)

        for gcf_id in gcf_ids:
            gcf_id_str = str(gcf_id)
            if gcf_id_str not in gcf_taxonomy:
                gcf_taxonomy[gcf_id_str] = {'genera': {}, 'species': {}}

            gcf_taxonomy[gcf_id_str]['genera'][genus] = gcf_taxonomy[gcf_id_str]['genera'].get(genus, 0) + 1
            gcf_taxonomy[gcf_id_str]['species'][species] = gcf_taxonomy[gcf_id_str]['species'].get(species, 0) + 1

    # Build heatmap data (GCF × Genus)
    heatmap_data = {}
    for gcf_id, tax_dist in gcf_taxonomy.items():
        heatmap_data[gcf_id] = tax_dist['genera']

    # Get top GCFs by member count
    top_gcfs = sorted(
        [(gcf_id, family_metadata.get(gcf_id, {}).get('member_count', 0), family_metadata.get(gcf_id, {}).get('product', 'Unknown'))
         for gcf_id in gcf_taxonomy.keys()],
        key=lambda x: -x[1]
    )[:50]  # Top 50 GCFs

    return {
        'heatmap_data': heatmap_data,
        'gcf_taxonomy': gcf_taxonomy,
        'genera': sorted(all_genera),
        'species': sorted(all_species),
        'top_gcfs': top_gcfs,
        'family_metadata': family_metadata
    }

def generate_bgc_distribution_html(gcf_data, taxonomy_map, gtdbtk_summary_path=None, gcf_heatmap_b64=None):
    """
    Generate HTML for BGC Distribution tab (replaces Tree View).
    Shows GCF × Taxonomy heatmap and distribution insights.
    """
    # Try to use GTDB-Tk taxonomy if available (more accurate)
    taxonomy_source = "NCBI"
    if gtdbtk_summary_path:
        try:
            gtdbtk_df = pd.read_csv(gtdbtk_summary_path, sep='\t')
            # Build taxonomy map from GTDB-Tk classification
            gtdb_taxonomy = {}
            for _, row in gtdbtk_df.iterrows():
                user_genome = row.get('user_genome', '')
                classification = row.get('classification', '')
                if user_genome and classification:
                    # Parse GTDB classification string: d__Bacteria;p__Proteobacteria;c__...;g__Genus;s__Species
                    parts = classification.split(';')
                    lineage = {}
                    for part in parts:
                        if part.startswith('g__'):
                            lineage['genus'] = part[3:] or 'Unknown'
                        elif part.startswith('s__'):
                            lineage['species'] = part[3:] or 'Unknown'

                    # GTDB-Tk prefixes user genomes with "usr_" to avoid clashes with
                    # reference sequences. Strip it to match genome_gcf_mapping keys.
                    genome_name = user_genome[4:] if user_genome.startswith('usr_') else user_genome
                    tax_entry = {
                        'lineage': {
                            'genus': {'name': lineage.get('genus', 'Unknown')},
                            'species': {'name': lineage.get('species', 'Unknown')}
                        }
                    }
                    gtdb_taxonomy[genome_name] = tax_entry
                    # Also index by assembly_id as fallback
                    assembly_id = extract_assembly_id_from_genome_name(genome_name)
                    if assembly_id:
                        gtdb_taxonomy[assembly_id] = tax_entry

            if gtdb_taxonomy:
                taxonomy_map = gtdb_taxonomy
                taxonomy_source = "GTDB-Tk"
                print(f"Using GTDB-Tk taxonomy for {len(gtdb_taxonomy)} genomes")
        except Exception as e:
            print(f"Warning: Could not load GTDB-Tk taxonomy: {e}")

    dist_data = build_gcf_taxonomy_distribution(gcf_data, taxonomy_map)

    if not dist_data:
        return '''
        <div class="info-box warning">
            <p><strong>BGC Distribution analysis requires:</strong></p>
            <ul>
                <li>BiG-SCAPE clustering results (GCF assignments)</li>
                <li>Taxonomy data (from NCBI or GTDB-Tk)</li>
            </ul>
            <p>Run with <code>--clustering bigscape</code> to enable GCF analysis.</p>
        </div>
        '''

    genera = dist_data['genera']
    top_gcfs = dist_data['top_gcfs']
    gcf_taxonomy = dist_data['gcf_taxonomy']
    family_metadata = dist_data['family_metadata']

    # Limit to top genera by total BGC count
    genus_totals = {}
    for gcf_id, tax_dist in gcf_taxonomy.items():
        for genus, count in tax_dist['genera'].items():
            genus_totals[genus] = genus_totals.get(genus, 0) + count

    top_genera = sorted(genus_totals.items(), key=lambda x: -x[1])[:20]
    top_genera_names = [g[0] for g in top_genera]

    # Build heatmap data for JavaScript
    heatmap_rows = []
    for gcf_id, member_count, product in top_gcfs[:30]:  # Top 30 GCFs for heatmap
        row_data = []
        for genus in top_genera_names:
            count = gcf_taxonomy.get(gcf_id, {}).get('genera', {}).get(genus, 0)
            row_data.append(count)
        heatmap_rows.append({
            'gcf_id': gcf_id,
            'product': product[:30],
            'member_count': member_count,
            'values': row_data
        })

    heatmap_json = json.dumps({
        'rows': heatmap_rows,
        'columns': top_genera_names
    })

    # Calculate distribution statistics
    gcf_specificity = []  # GCFs that are genus-specific
    gcf_widespread = []   # GCFs found in many genera

    for gcf_id, tax_dist in gcf_taxonomy.items():
        genera_with_gcf = [g for g, c in tax_dist['genera'].items() if c > 0]
        num_genera = len(genera_with_gcf)
        total_bgcs = sum(tax_dist['genera'].values())

        meta = family_metadata.get(gcf_id, {})
        product = meta.get('product', 'Unknown')

        if num_genera == 1 and total_bgcs >= 3:
            gcf_specificity.append({
                'gcf_id': gcf_id,
                'genus': genera_with_gcf[0],
                'count': total_bgcs,
                'product': product
            })
        elif num_genera >= 5:
            gcf_widespread.append({
                'gcf_id': gcf_id,
                'num_genera': num_genera,
                'count': total_bgcs,
                'product': product
            })

    # Sort by count
    gcf_specificity.sort(key=lambda x: -x['count'])
    gcf_widespread.sort(key=lambda x: -x['num_genera'])

    # Build specificity insights HTML
    specificity_rows = ''
    for item in gcf_specificity[:10]:
        specificity_rows += f'''
            <tr>
                <td style="padding: 6px 10px;">GCF-{item['gcf_id']}</td>
                <td style="padding: 6px 10px;">{item['product'][:35]}</td>
                <td style="padding: 6px 10px; font-weight: bold;">{item['genus']}</td>
                <td style="padding: 6px 10px; text-align: right;">{item['count']}</td>
            </tr>
        '''

    widespread_rows = ''
    for item in gcf_widespread[:10]:
        widespread_rows += f'''
            <tr>
                <td style="padding: 6px 10px;">GCF-{item['gcf_id']}</td>
                <td style="padding: 6px 10px;">{item['product'][:35]}</td>
                <td style="padding: 6px 10px; text-align: right;">{item['num_genera']}</td>
                <td style="padding: 6px 10px; text-align: right;">{item['count']}</td>
            </tr>
        '''

    return f'''
    <p style="color: #666; margin-bottom: 20px;">
        <em>Analysis of BGC distribution across taxonomic groups. Taxonomy source: <strong>{taxonomy_source}</strong>.
        Shows which Gene Cluster Families (GCFs) are taxon-specific vs widespread.</em>
    </p>

    <h3>GCF × Genus Heatmap</h3>
    <p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">
        {'Presence/absence of each Gene Cluster Family across genera. Rows: GCFs ordered by Jaccard-distance hierarchical clustering; columns: genera ordered by GTDB-Tk phylogeny.' if gcf_heatmap_b64 else 'Top 30 GCFs (rows) vs top 20 genera (columns). Color intensity = number of BGCs.'}
    </p>
    {'<div style="margin-bottom: 30px;"><img src="data:image/svg+xml;base64,' + gcf_heatmap_b64 + '" alt="GCF × Genus Heatmap" style="max-width: 100%; height: auto; display: block;"></div>' if gcf_heatmap_b64 else f"""
    <div id="heatmap-container" style="width: 100%; overflow-x: auto; margin-bottom: 30px;">
        <canvas id="heatmap-canvas" style="max-width: 100%;"></canvas>
    </div>
    <script>
        (function() {{
            const data = {heatmap_json};
            const canvas = document.getElementById('heatmap-canvas');
            const ctx = canvas.getContext('2d');
            const cellWidth = 45; const cellHeight = 22;
            const labelWidth = 180; const labelHeight = 120;
            const rows = data.rows; const cols = data.columns;
            canvas.width = labelWidth + cols.length * cellWidth + 80;
            canvas.height = labelHeight + rows.length * cellHeight + 20;
            let maxVal = 1;
            rows.forEach(row => {{ row.values.forEach(v => {{ if (v > maxVal) maxVal = v; }}); }});
            function getColor(value) {{
                if (value === 0) return '#f8f9fa';
                const intensity = Math.min(1, value / maxVal);
                return `rgb(${{Math.round(255 - intensity*212)}}, ${{Math.round(255 - intensity*165)}}, ${{Math.round(255 - intensity*95)}})`;
            }}
            ctx.fillStyle = 'white'; ctx.fillRect(0, 0, canvas.width, canvas.height);
            ctx.save(); ctx.font = '10px sans-serif'; ctx.fillStyle = '#333';
            cols.forEach((col, i) => {{
                ctx.save(); ctx.translate(labelWidth + i*cellWidth + cellWidth/2, labelHeight-5);
                ctx.rotate(-Math.PI/3); ctx.textAlign = 'left';
                ctx.fillText(col.length > 15 ? col.slice(0,15)+'...' : col, 0, 0); ctx.restore();
            }}); ctx.restore();
            rows.forEach((row, i) => {{
                const y = labelHeight + i * cellHeight;
                ctx.font = '10px sans-serif'; ctx.fillStyle = '#333'; ctx.textAlign = 'right';
                ctx.fillText(`GCF-${{row.gcf_id}} ${{row.product}}`, labelWidth-5, y+cellHeight/2+3);
                row.values.forEach((value, j) => {{
                    const x = labelWidth + j * cellWidth;
                    ctx.fillStyle = getColor(value); ctx.fillRect(x, y, cellWidth-1, cellHeight-1);
                    if (value > 0) {{
                        ctx.fillStyle = value > maxVal*0.5 ? 'white' : '#333';
                        ctx.font = '9px sans-serif'; ctx.textAlign = 'center';
                        ctx.fillText(value.toString(), x+cellWidth/2, y+cellHeight/2+3);
                    }}
                }});
                ctx.fillStyle='#666'; ctx.font='9px sans-serif'; ctx.textAlign='left';
                ctx.fillText(`(${{row.member_count}})`, labelWidth+cols.length*cellWidth+5, y+cellHeight/2+3);
            }});
            const legendY = labelHeight + rows.length*cellHeight + 10;
            ctx.font='10px sans-serif'; ctx.fillStyle='#666'; ctx.textAlign='left';
            ctx.fillText('Color: BGC count (darker = more)', labelWidth, legendY);
        }})();
    </script>""" }

    <div style="display: flex; gap: 30px; flex-wrap: wrap; margin-top: 20px;">
        <div style="flex: 1; min-width: 400px;">
            <h3>Genus-Specific GCFs</h3>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">
                GCFs found in only one genus (potential taxon-specific metabolites)
            </p>
            <table style="width: 100%; border-collapse: collapse; font-size: 0.9em;">
                <thead>
                    <tr style="background: #f0f4f8;">
                        <th style="padding: 8px 10px; text-align: left;">GCF</th>
                        <th style="padding: 8px 10px; text-align: left;">Product</th>
                        <th style="padding: 8px 10px; text-align: left;">Genus</th>
                        <th style="padding: 8px 10px; text-align: right;">BGCs</th>
                    </tr>
                </thead>
                <tbody>
                    {specificity_rows if specificity_rows else '<tr><td colspan="4" style="padding: 10px; color: #666;">No genus-specific GCFs found</td></tr>'}
                </tbody>
            </table>
        </div>

        <div style="flex: 1; min-width: 400px;">
            <h3>Widespread GCFs</h3>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">
                GCFs found across 5+ genera (conserved or horizontally transferred)
            </p>
            <table style="width: 100%; border-collapse: collapse; font-size: 0.9em;">
                <thead>
                    <tr style="background: #f0f4f8;">
                        <th style="padding: 8px 10px; text-align: left;">GCF</th>
                        <th style="padding: 8px 10px; text-align: left;">Product</th>
                        <th style="padding: 8px 10px; text-align: right;">Genera</th>
                        <th style="padding: 8px 10px; text-align: right;">BGCs</th>
                    </tr>
                </thead>
                <tbody>
                    {widespread_rows if widespread_rows else '<tr><td colspan="4" style="padding: 10px; color: #666;">No widespread GCFs found</td></tr>'}
                </tbody>
            </table>
        </div>
    </div>

    '''


# Coupling class display metadata (class_id → (display_name, marker, pathway, reference_genes))
