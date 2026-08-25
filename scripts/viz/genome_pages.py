"""Per-genome HTML pages linked from the main BGC report."""

import json
import re
from pathlib import Path

import pandas as pd


def create_genome_metadata_pages(counts_file, assembly_info, name_map, outdir, taxon, taxonomy_map_data=None, tabulation_file=None):
    '''Create individual HTML pages for each genome with metadata and KCB data'''
    genome_dir = outdir / 'genomes'
    genome_dir.mkdir(exist_ok=True)

    # Clean taxon name for URL - match Nextflow sanitizeTaxon function
    taxon_clean = re.sub(r'[^a-zA-Z0-9_]', '_', taxon)
    taxon_clean = re.sub(r'_+', '_', taxon_clean).strip('_')

    # Read counts to get genome list
    counts_df = pd.read_csv(counts_file, sep='\t', skiprows=lambda i: i == 0)

    # Read assembly info
    assembly_df = pd.read_csv(assembly_info, sep='\t')

    # Load tabulation data for KCB information
    tab_df = None
    if tabulation_file:
        try:
            tab_df = pd.read_csv(tabulation_file, sep='\t')
            # Fill NaN values with empty strings for consistent filtering
            for col in ['KCB_hit', 'KCB_acc', 'KCB_sim']:
                if col in tab_df.columns:
                    tab_df[col] = tab_df[col].fillna('')
        except Exception as e:
            print(f"Warning: Could not load tabulation file: {e}")
            tab_df = None

    # Load name map and create reverse mapping (genome_name -> assembly_id)
    with open(name_map, 'r') as f:
        name_map_data = json.load(f)
    reverse_map = {v: k for k, v in name_map_data.items()}

    for _, row in counts_df.iterrows():
        genome_name = row['record'].replace('.gbff', '')

        # Look up assembly ID using reverse name map
        assembly_id = reverse_map.get(genome_name, None)

        # Find matching assembly metadata
        metadata_rows = assembly_df[assembly_df.iloc[:, 0] == assembly_id] if assembly_id else pd.DataFrame()

        if not metadata_rows.empty:
            metadata = metadata_rows.iloc[0]

            def get_value(idx):
                if len(metadata) > idx:
                    val = metadata.iloc[idx]
                    if pd.isna(val) or val == '' or str(val) == 'nan':
                        return 'N/A'
                    return str(val)
                return 'N/A'

            organism_name = get_value(1)
            strain = get_value(2)
            isolation_source = get_value(3)
            isolate = get_value(4)
            notes = get_value(5)
            checkm_completeness = get_value(6)
            checkm_contamination = get_value(7)
            organism_tax_id = get_value(8)

            # Display taxonomy lineage if available
            taxonomy_html = ''
            if taxonomy_map_data and assembly_id and assembly_id in taxonomy_map_data:
                genome_taxonomy = taxonomy_map_data[assembly_id]
                lineage = genome_taxonomy.get('lineage', {})

                lineage_parts = []
                for rank in ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                    if rank in lineage:
                        tax_node = lineage[rank]
                        name = tax_node.get('name', 'N/A')
                        if name not in ['N/A', '', 'nan', 'None']:
                            lineage_parts.append(f"<strong>{rank.capitalize()}:</strong> {name}")

                if lineage_parts:
                    lineage_html = '<br>'.join(lineage_parts)
                    taxonomy_html = f'''
        <tr><td colspan="2"><strong>Taxonomic Lineage</strong></td></tr>
        <tr><td colspan="2">{lineage_html}</td></tr>'''

            assembly_table = f'''
    <table>
        <tr><th>Field</th><th>Value</th></tr>
        <tr><td>Assembly Accession</td><td>{assembly_id}</td></tr>
        <tr><td>Organism Tax ID</td><td>{organism_tax_id}</td></tr>
        <tr><td>Organism Name</td><td>{organism_name}</td></tr>
        <tr><td>Strain</td><td>{strain}</td></tr>
        <tr><td>Isolate</td><td>{isolate}</td></tr>
        <tr><td>Isolation Source</td><td>{isolation_source}</td></tr>
        <tr><td>Assembly Notes</td><td>{notes}</td></tr>
        <tr><td>CheckM Completeness</td><td>{checkm_completeness}</td></tr>
        <tr><td>CheckM Contamination</td><td>{checkm_contamination}</td></tr>
        {taxonomy_html}
    </table>'''
        else:
            assembly_table = '<p><em>Assembly metadata not found for this genome</em></p>'

        # Build KCB (KnownClusterBlast) section if tabulation data available
        kcb_section = ''
        if tab_df is not None and 'KCB_hit' in tab_df.columns:
            # Filter tabulation data for this genome
            genome_tab = tab_df[tab_df['file'] == genome_name]
            total_regions = len(genome_tab)
            kcb_hits = genome_tab[genome_tab['KCB_hit'] != '']
            regions_with_hits = len(kcb_hits)
            regions_without_hits = total_regions - regions_with_hits
            unique_clusters = kcb_hits['KCB_hit'].nunique() if not kcb_hits.empty else 0

            kcb_section = f'''
    <h2>Known Cluster Matches (KnownClusterBlast)</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
        <tr><td>Total BGC Regions</td><td>{total_regions}</td></tr>
        <tr><td>Regions with Known Cluster Matches</td><td>{regions_with_hits}</td></tr>
        <tr><td>Regions without Known Cluster Matches (Potentially Novel)</td><td>{regions_without_hits}</td></tr>
        <tr><td>Unique Known Clusters Matched</td><td>{unique_clusters}</td></tr>
    </table>
    <p style="color: #666; margin-top: 15px;"><em>View detailed cluster matches in the antiSMASH results for this genome.</em></p>
    '''

        html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <title>Genome Metadata - {genome_name}</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; background: #f8f9fa; color: #333; }}
        h1 {{ color: #333; }}
        h2 {{ color: #333; margin-top: 30px; border-bottom: 2px solid #5b8ac5; padding-bottom: 10px; }}
        table {{ border-collapse: collapse; width: 100%; background: white; margin-top: 20px; border: 1px solid #ddd; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background: #2c5aa0; color: white; font-weight: bold; }}
        tr:nth-child(even) {{ background-color: #f8f9fa; }}
        tr:hover {{ background-color: #e9ecef; }}
        .back-link {{ margin-top: 20px; display: inline-block; padding: 10px 20px; background: #2c5aa0; color: white; font-weight: bold; text-decoration: none; border-radius: 4px; }}
        .back-link:hover {{ background: #5b8ac5; }}
        .antismash-link {{ display: inline-block; padding: 10px 20px; background: #28a745; color: white; font-weight: bold; text-decoration: none; border-radius: 4px; margin-right: 10px; }}
        .antismash-link:hover {{ background: #34ce57; color: white; text-decoration: none; }}
        .button-container {{ margin-top: 30px; }}
        h3 {{ color: #333; margin-top: 25px; }}
        table a {{ color: #2c5aa0; text-decoration: none; }}
        table a:hover {{ text-decoration: underline; }}
    </style>
</head>
<body>
    <h1>Genome Metadata: {genome_name}</h1>

    <h2>BGC Statistics</h2>
    <table>
        <tr><th>Field</th><th>Value</th></tr>
        <tr><td>Total BGC Count</td><td>{row.get('total_count', 'N/A')}</td></tr>
    </table>

    <h2>Assembly Metadata</h2>
    {assembly_table}

    {kcb_section}

    <div class="button-container">
        <a href="../../../antismash_results/{taxon_clean}/{genome_name}/index.html" class="antismash-link" target="_blank">View antiSMASH Results</a>
        <a href="../bgc_report.html" class="back-link">← Back to Main Report</a>
    </div>
</body>
</html>
'''

        output_file = genome_dir / f'{genome_name}.html'
        with open(output_file, 'w') as f:
            f.write(html_content)

    return len(counts_df)
