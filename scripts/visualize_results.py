#!/usr/bin/env python3
"""
BGC Analysis Visualization Script

Generates interactive HTML reports and visualizations for biosynthetic gene cluster analysis.
Includes:
- Tab-based HTML report with sections (Overview, Taxonomy, BGC Distribution, Genomes, Clustering)
- Summary statistics dashboard with tabulation stats
- Taxonomic distribution tree (interactive expandable)
- GCF × Taxonomy heatmap and distribution analysis (replaces tree visualization)
- Donut chart for BGC type distribution
- BGC histogram showing distribution across genomes
- Searchable genome table
- Individual genome metadata pages
- Clustering statistics (BiG-SCAPE)
"""

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from utils.parsers import parse_trace_file
from utils.report_lint import check_report
from utils.trace import aggregate_trace_by_process, generate_resource_usage_html
from viz.clustering import generate_bigscape_stats_html, generate_gcf_visualization_html
from viz.tables import (calculate_summary_statistics, create_bgc_distribution_table,
                        generate_genome_table_html)
from viz.taxonomy import generate_taxonomy_tree_html
from viz.tree_viz import prepare_phylo_tree_for_js
from viz.report_assets import REPORT_CSS, REPORT_JS
from viz.distribution import generate_bgc_distribution_html
from viz.genome_pages import create_genome_metadata_pages
from viz.rarefaction import generate_rarefaction_curve
from viz.report_sections import (_build_bigscape_section_html,
                                 _build_kcb_content, _build_rarefaction_section,
                                 gcf_coupling_classes, build_gcf_support_rows,
                                 build_overview_stats,
                                 _build_versions_html, build_coupling_table_rows)


def generate_html_report(outdir, taxon, table_header, table_rows, stats, tree_html='',
                         bigscape_stats_html='', gcf_visualization_html='',
                         phylo_tree_generated=False, genome_table=None,
                         resource_usage_html='', phylo_tree_data=None, gcf_data=None, taxonomy_map=None,
                         versions_data=None, rarefaction_stats=None,
                         gtdbtk_summary_path=None, gcf_tree_b64=None,
                         gcf_tree_mime='image/png',
                         all_bgcs_tree_b64=None, all_bgcs_tree_mime='image/png',
                         gcf_heatmap_b64=None,
                         coupling_table_rows=None, gcf_classes=None,
                         gcf_support_rows=None):
    '''Generate tab-based HTML report combining all visualizations'''

    # Clean taxon name for URLs - match Nextflow sanitizeTaxon function
    taxon_clean = re.sub(r'[^a-zA-Z0-9_]', '_', taxon)
    taxon_clean = re.sub(r'_+', '_', taxon_clean).strip('_')

    kcb_stats = stats.get('kcb_stats', {})
    kcb = _build_kcb_content(kcb_stats, taxon_clean, gcf_data, gcf_classes)
    kcb_mapping_section    = kcb['kcb_mapping_section']
    novel_bgcs_tab_content = kcb['novel_bgcs_tab_content']
    kcb_hits_tab_content   = kcb['kcb_hits_tab_content']

    overview_stats = build_overview_stats(stats, kcb_stats, gcf_data, rarefaction_stats)

    # Provenance under the title. A report that gets shared or archived should say when it
    # was made and over what — the software versions in the Pipeline tab do not answer
    # "how many genomes was this?" or "when?".
    _generated = datetime.now().strftime('%Y-%m-%d %H:%M')
    _n_genomes = stats.get('total_genomes', 0) or 0
    _versions = versions_data or {}
    _as_ver = ''
    for _k, _v in (_versions.items() if isinstance(_versions, dict) else []):
        if 'antismash' in str(_k).lower():
            _as_ver = f' · antiSMASH {_v}'
            break
    provenance = (f'{_n_genomes:,} genomes · generated {_generated}{_as_ver}'
                  if _n_genomes else f'generated {_generated}{_as_ver}')

    bigscape_section_html = _build_bigscape_section_html(bigscape_stats_html, gcf_visualization_html,
                                                        taxon_clean, gcf_support_rows)

    versions_html = _build_versions_html(versions_data)

    rarefaction_section = _build_rarefaction_section(rarefaction_stats)

    # Build BGC Distribution section HTML (replaces Tree View)
    distribution_section = generate_bgc_distribution_html(gcf_data, taxonomy_map, gtdbtk_summary_path, gcf_heatmap_b64=gcf_heatmap_b64)

    # Use distribution_section for the "BGC Distribution" tab (replaces old tree view)
    tree_section = distribution_section  # Keep variable name for template compatibility

    # Build genome table section (avoid nested f-strings)
    # Only the first slice is in the HTML; the rest ships as JSON and renders on demand
    # — see viz/tables.generate_genome_table_html for why.
    _gt = genome_table or {}
    genome_table_rows = (_gt.get('initial_rows')
                         or '<tr><td colspan="6">No genome data available</td></tr>')
    genome_data_json = _gt.get('data_json', '[]')
    genome_total = _gt.get('total', 0)
    genome_shown = _gt.get('shown', 0)

    # Coupling table rows: use dynamic data when available, otherwise placeholder
    if coupling_table_rows is None:
        coupling_table_rows = '<tr><td colspan="4" style="padding: 7px 12px; color: #999;">Coupling enzyme data not available. Run the pipeline with <code>--clustering bigscape</code>.</td></tr>'

    html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <!-- Without this, mobile browsers lay the page out against a ~980px virtual
         viewport and the max-width media queries below never fire. -->
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>BGC Analysis Report - {taxon}</title>
    <style>{REPORT_CSS}</style>
</head>
<body>
    <h1>BGC Analysis Report</h1>
    <p class="subtitle">Taxon: <strong>{taxon}</strong></p>
    <p style="text-align: center; color: #888; font-size: 0.85em; margin: -6px 0 4px;">
        {provenance}
    </p>

    <div class="tabs">
        <input type="radio" id="tab1" name="tabs" checked>
        <label for="tab1">Overview</label>

        <input type="radio" id="tab2" name="tabs">
        <label for="tab2">Phylogeny</label>

        <input type="radio" id="tab3" name="tabs">
        <label for="tab3">Genomes</label>

        <input type="radio" id="tab4" name="tabs">
        <label for="tab4">GCF Analysis</label>

        <input type="radio" id="tab5" name="tabs">
        <label for="tab5">Novel BGCs</label>

        <input type="radio" id="tab6" name="tabs">
        <label for="tab6">KCB Hits</label>

        <input type="radio" id="tab7" name="tabs">
        <label for="tab7">Pipeline</label>

        <!-- Tab 1: Overview -->
        <div class="tab-content" id="content1">
            <p style="color: #666; font-size: 0.9em; margin: 4px 0 2px;">
                <em>Detection is restricted to the antiSMASH <strong>phosphonate</strong> rule
                (<code>--hmmdetection-limit-to-rule-names phosphonate</code>), so every region below is a
                phosphonate BGC and no other BGC class was searched for. "No MIBiG match" should be read
                against that: MIBiG holds few characterised phosphonate pathways, so a miss is expected
                and is weaker evidence of novelty than it would be for a well-represented class.</em>
            </p>
            {overview_stats}
            {kcb_mapping_section}
            {rarefaction_section}

        </div>

        <!-- Tab 2: Phylogeny (Taxonomy tree + GTDB-Tk BGC distribution) -->
        <div class="tab-content" id="content2">
            <h3>Taxonomic Distribution of BGCs</h3>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Expandable NCBI taxonomy tree showing BGC statistics at each taxonomic level. Click on nodes to expand/collapse.
                Species nodes expand to show individual genomes with their BGC counts.</em>
            </p>
            {tree_html if tree_html else '<div class="info-box warning"><p>Taxonomy tree data not available.</p></div>'}

            <hr class="tab-section-divider">

            <h3>GCF Distribution Across Taxa</h3>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Which Gene Cluster Families are confined to one genus and which are widespread,
                using GTDB-Tk taxonomy where available. The GTDB-Tk tree itself is not drawn here —
                the Newick files are published under <code>gtdbtk_results/</code> for iTOL, FigTree
                or Dendroscope.</em>
            </p>
            {tree_section}
        </div>

        <!-- Tab 3: Genomes -->
        <div class="tab-content" id="content3">
            <h2>All Genomes</h2>
            <p style="color: #666; margin-bottom: 15px;">
                <em>Searchable table of all analyzed genomes. Click genome names for detailed metadata pages.</em>
            </p>
            <div class="search-box">
                <input type="text" id="genomeSearch" placeholder="Search by genome, strain, assembly or taxonomy" onkeyup="filterGenomes()">
            </div>
            <p id="genomeTableStatus" style="color: #777; font-size: 0.85em; margin: 2px 0 10px;">
                Showing {genome_shown} of {genome_total} genomes.
                <button type="button" onclick="showAllGenomes()"
                        style="background: none; border: none; color: #2c5aa0; cursor: pointer;
                               padding: 0; font: inherit; text-decoration: underline;">Show all</button>
            </p>
            <div class="table-container">
                <table id="genomeTable">
                    <thead>
                        <tr>
                            <th>Genome Name</th>
                            <th>Assembly ID</th>
                            <th>Organism</th>
                            <th>Taxonomy</th>
                            <th>Total BGCs</th>
                            <th>Top BGC Types</th>
                        </tr>
                    </thead>
                    <tbody id="genomeTableBody">
                    {genome_table_rows}
                    </tbody>
                </table>
            </div>
            <script id="genomeData" type="application/json">{genome_data_json}</script>
        </div>

        <!-- Tab 4: GCF Analysis (Clustering + Coupling Enzyme) -->
        <div class="tab-content" id="content4">

            <h3>GCF Biosynthetic Phylogeny</h3>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Classification of phosphonate BGCs by the coupling enzyme acting on phosphonopyruvate — the branching step
                immediately downstream of PEP mutase that determines the downstream biosynthetic pathway.</em>
            </p>

            <!-- GCF medoid tree (top half) -->
            <div class="plot" style="margin-top: 20px;">
                <h4 style="margin: 0 0 8px 0; color: #333;">GCF-Level Tree</h4>
                <p style="color: #666; font-size: 0.85em; margin: 0 0 12px 0;">
                    One representative (medoid) per Gene Cluster Family. Circle size ∝ GCF membership.
                </p>
                {'<img src="data:' + gcf_tree_mime + ';base64,' + gcf_tree_b64 + '" alt="GCF Biosynthetic Phylogeny" style="max-width: 100%; height: auto; display: block;">' if gcf_tree_b64 else '<div style="color: #999; padding: 20px; background: #f8f9fa; border-radius: 8px; text-align: center; font-size: 0.9em;">GCF biosynthetic tree not generated.<br>Run with <code>--clustering bigscape</code> to enable.</div>'}
            </div>

            <hr class="tab-section-divider">

            <!-- All-BGCs circular tree (bottom half) -->
            <div class="plot">
                <h4 style="margin: 0 0 8px 0; color: #333;">All-BGCs Circular Tree</h4>
                <p style="color: #666; font-size: 0.85em; margin: 0 0 12px 0;">
                    Every BGC as a leaf. NJ tree from the full BiG-SCAPE pairwise distance matrix, colored by coupling enzyme class.
                </p>
                {'<img src="data:' + all_bgcs_tree_mime + ';base64,' + all_bgcs_tree_b64 + '" alt="All-BGCs Biosynthetic Tree" style="max-width: 100%; height: auto; display: block;">' if all_bgcs_tree_b64 else '<div style="color: #999; padding: 20px; background: #f8f9fa; border-radius: 8px; text-align: center; font-size: 0.9em;">All-BGCs tree not generated.<br>Run with <code>--clustering bigscape</code> to enable.</div>'}
            </div>

            <div style="margin-top: 24px; background: #f8f9fa; padding: 20px; border-radius: 10px;">
                <h4 style="margin-top: 0;">Coupling Enzyme Classes</h4>
                <table style="width: 100%; border-collapse: collapse; font-size: 0.9em;">
                    <thead>
                        <tr style="background: #e9ecef;">
                            <th style="text-align: left; padding: 8px 12px; border-bottom: 2px solid #dee2e6;">Class</th>
                            <th style="text-align: left; padding: 8px 12px; border-bottom: 2px solid #dee2e6;">Marker</th>
                            <th style="text-align: left; padding: 8px 12px; border-bottom: 2px solid #dee2e6;">Product</th>
                            <th style="text-align: left; padding: 8px 12px; border-bottom: 2px solid #dee2e6;">Reference genes</th>
                            <th style="text-align: left; padding: 8px 12px; border-bottom: 2px solid #dee2e6;">GCFs</th>
                        </tr>
                    </thead>
                    <tbody>
                        {coupling_table_rows}
                    </tbody>
                </table>
            </div>

            <hr class="tab-section-divider">

            {bigscape_section_html}
            {f'<div class="info-box" style="background-color: #f8f9fa; border-left: 4px solid #6c757d;"><p style="color: #666;">No clustering analysis was performed. To enable clustering, run the pipeline with <code>--clustering bigscape</code>.</p></div>' if not bigscape_section_html else ''}
        </div>

        <!-- Tab 5: Novel BGCs -->
        <div class="tab-content" id="content5">
            {novel_bgcs_tab_content}
        </div>

        <!-- Tab 6: KCB Hits -->
        <div class="tab-content" id="content6">
            {kcb_hits_tab_content}
        </div>

        <!-- Tab 7: Pipeline Info -->
        <div class="tab-content" id="content7">
            <h2>Pipeline Information</h2>
            <h3 style="margin-top: 20px;">Resource Usage</h3>
            <p style="color: #666; margin-bottom: 20px; font-size: 0.9em;">
                <em>Resource consumption metrics from Nextflow trace data, showing CPU, memory, and runtime for each pipeline process.</em>
            </p>
            {resource_usage_html if resource_usage_html else '<div class="info-box warning"><p>No resource usage data available. Trace data will appear here after running the pipeline.</p></div>'}
            {versions_html}
        </div>

    </div>

    <script>{REPORT_JS}</script>
</body>
</html>
'''

    # Lint before writing. The report's JS is assembled from several Python string
    # constants, so this is the first point at which all of it exists together — and
    # a handler wired to a missing function fails silently in the browser.
    problems = check_report(html_content)
    if problems:
        for problem in problems:
            print(f"ERROR: {problem}", file=sys.stderr)
        raise SystemExit("refusing to write a report with broken JavaScript handlers")

    with open(f'{outdir}/bgc_report.html', 'w', encoding='utf-8') as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Visualize BGC analysis results')
    parser.add_argument('--counts', type=Path, help='Path to region_counts.tsv')
    parser.add_argument('--tabulation', type=Path, help='Path to region_tabulation.tsv')
    parser.add_argument('--assembly_info', type=Path, help='Path to assembly_info_table.txt')
    parser.add_argument('--name_map', type=Path, help='Path to name_map.json')
    parser.add_argument('--taxonomy_map', type=Path, help='Path to taxonomy_map.json')
    parser.add_argument('--taxonomy_tree', type=Path, help='Path to taxonomy_tree.json')
    parser.add_argument('--phylo_tree', type=Path, help='Path to Newick phylogenetic tree file from GTDB-Tk')
    parser.add_argument('--gtdbtk_summary', type=Path, help='Path to GTDB-Tk summary TSV file')
    parser.add_argument('--bigscape_stats', type=Path, help='Path to bigscape_statistics.json')
    parser.add_argument('--bigscape_db', type=Path, help='Path to BiG-SCAPE SQLite database for rarefaction curve')
    parser.add_argument('--gcf_data', type=Path, help='Path to gcf_representatives.json')
    parser.add_argument('--trace', type=Path, help='Path to Nextflow pipeline_trace.tsv file')
    parser.add_argument('--outdir', type=Path, required=True, help='Output directory for plots')
    parser.add_argument('--taxon', type=str, default='Unknown', help='Taxon name for report')
    parser.add_argument('--mibig_included', action='store_true', help='Whether MIBiG references were included in BiG-SCAPE analysis')
    parser.add_argument('--versions', type=Path, help='Path to software_versions.json')
    parser.add_argument('--skip_tree', action='store_true', help='Skip phylogenetic tree visualization (useful for very large datasets)')
    parser.add_argument('--outgroup', type=str, help='Outgroup taxon pattern for tree pruning (e.g., "g__Escherichia")')
    parser.add_argument('--gcf_tree', type=Path, help='Path to GCF biosynthetic NJ tree PNG from GCF_BIOSYNTHETIC_TREE')
    parser.add_argument('--gcf_tree_svg', type=Path, help='Path to GCF biosynthetic NJ tree SVG (preferred over PNG for quality)')
    parser.add_argument('--all_bgcs_tree', type=Path, help='Path to all-BGCs circular NJ tree PNG from GCF_BIOSYNTHETIC_TREE')
    parser.add_argument('--all_bgcs_tree_svg', type=Path, help='Path to all-BGCs circular NJ tree SVG (preferred over PNG: vector, and ~45%% smaller once base64-encoded)')
    parser.add_argument('--gcf_heatmap_svg', type=Path, help='Path to GCF × species heatmap SVG from GCF_BIOSYNTHETIC_TREE')
    parser.add_argument('--coupling_annotation', type=Path, help='Path to phosphonate_itol_coupling.txt from GCF_BIOSYNTHETIC_TREE')
    parser.add_argument('--coupling_support', type=Path,
                        help='phosphonate_coupling_support.tsv from GCF_BIOSYNTHETIC_TREE')
    parser.add_argument('--seed', type=int, default=0, help='RNG seed for the rarefaction resampling; fixed by default so reports are reproducible')

    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load taxonomy data if available
    taxonomy_map_data = None
    if args.taxonomy_map and args.taxonomy_map.exists():
        print(f"Loading taxonomy map...")
        with open(args.taxonomy_map, 'r') as f:
            taxonomy_map_data = json.load(f)

    taxonomy_tree_data = None
    if args.taxonomy_tree and args.taxonomy_tree.exists():
        print(f"Loading taxonomy tree...")
        with open(args.taxonomy_tree, 'r') as f:
            taxonomy_tree_data = json.load(f)

    table_header = ''
    table_rows = ''
    stats = {}
    tree_html = ''
    genome_table = None

    if args.counts:
        print(f"Generating count visualizations...")

        # Calculate summary statistics (including tabulation stats if available)
        print(f"Calculating summary statistics...")
        tabulation_path = str(args.tabulation) if args.tabulation and args.tabulation.exists() else None
        stats = calculate_summary_statistics(args.counts, tabulation_path)

        # Create genome metadata pages and genome table HTML
        if args.assembly_info and args.name_map:
            print(f"Creating genome metadata pages...")
            num_pages = create_genome_metadata_pages(args.counts, args.assembly_info, args.name_map, args.outdir, args.taxon, taxonomy_map_data, tabulation_path)
            print(f"Created {num_pages} genome metadata pages")

            # Generate searchable genome table for Genomes tab
            print(f"Generating genome table HTML...")
            genome_table = generate_genome_table_html(args.counts, args.assembly_info,
                                                      args.name_map, taxonomy_map_data)

        # Create interactive table
        print(f"Creating interactive BGC distribution table...")
        table_header, table_rows = create_bgc_distribution_table(args.counts, args.outdir)

    # Generate taxonomy tree visualization
    if taxonomy_tree_data:
        print(f"Generating taxonomy tree visualization...")
        tree_html = generate_taxonomy_tree_html(taxonomy_tree_data)

    # Generate phylogenetic tree data for JavaScript visualization
    phylo_tree_generated = False
    phylo_tree_data = None
    gtdbtk_summary_path = str(args.gtdbtk_summary) if args.gtdbtk_summary and args.gtdbtk_summary.exists() else None
    if args.skip_tree:
        print("Skipping phylogenetic tree visualization (--skip_tree enabled)")
        phylo_tree_data = {'skipped': True, 'reason': 'user_disabled'}
    elif args.phylo_tree and args.phylo_tree.exists():
        print(f"Preparing phylogenetic tree for visualization...")
        gtdbtk_summary_path = str(args.gtdbtk_summary) if args.gtdbtk_summary and args.gtdbtk_summary.exists() else None
        counts_path = str(args.counts) if args.counts and args.counts.exists() else None

        # Prepare pruned tree data for JavaScript visualization
        phylo_tree_data = prepare_phylo_tree_for_js(
            str(args.phylo_tree),
            gtdbtk_summary_path,
            counts_path,
            str(args.outdir),
            outgroup=args.outgroup
        )

        if phylo_tree_data:
            phylo_tree_generated = True
            print(f"Prepared interactive tree with {phylo_tree_data.get('leaf_count', 0)} genomes")

    # Generate BiG-SCAPE statistics visualization
    bigscape_stats_html = ''
    rarefaction_stats = None
    if args.bigscape_stats and args.bigscape_stats.exists():
        print(f"Generating BiG-SCAPE statistics visualization...")
        bigscape_stats_html = generate_bigscape_stats_html(str(args.bigscape_stats), args.mibig_included)

        # Generate rarefaction curve from BiG-SCAPE database
        bigscape_db = None
        if args.bigscape_db and args.bigscape_db.exists():
            bigscape_db = args.bigscape_db
        else:
            # Fallback: try to find .db file in the same directory as stats
            bigscape_dir = args.bigscape_stats.parent
            db_files = list(bigscape_dir.glob('*.db'))
            if db_files:
                bigscape_db = db_files[0]

        if bigscape_db:
            print(f"Generating rarefaction curve from {bigscape_db.name}...")
            rarefaction_stats = generate_rarefaction_curve(str(bigscape_db), args.outdir, args.taxon,
                                                           seed=args.seed, counts_file=args.counts)
            if rarefaction_stats:
                _c = rarefaction_stats["chao2"]
                print(f"Rarefaction curve generated: {rarefaction_stats['total_gcfs']} GCFs observed, "
                      f"{_c['s_est']:.0f} estimated (Chao2), {_c['coverage']:.0f}% coverage "
                      f"over {rarefaction_stats['n_genomes']} genomes "
                      f"[{rarefaction_stats['denominator']}]")
        else:
            print("Warning: BiG-SCAPE database not found, skipping rarefaction curve generation")

    # Generate GCF visualization HTML and load GCF data for overview
    gcf_visualization_html = ''
    gcf_data_dict = None
    if args.gcf_data and args.gcf_data.exists():
        print(f"Generating GCF visualization...")
        gcf_visualization_html = generate_gcf_visualization_html(str(args.gcf_data), args.taxon)
        # Also load as dict for overview sections
        try:
            with open(args.gcf_data, 'r') as f:
                gcf_data_dict = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load GCF data for overview: {e}")

    # Load taxonomy map for overview sections
    taxonomy_map_dict = None
    if args.taxonomy_map and args.taxonomy_map.exists():
        try:
            with open(args.taxonomy_map, 'r') as f:
                taxonomy_map_dict = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load taxonomy map: {e}")

    # Generate resource usage HTML from trace file
    resource_usage_html = ''
    if args.trace and args.trace.exists():
        print(f"Processing pipeline trace data...")
        tasks = parse_trace_file(str(args.trace))
        if tasks:
            processes = aggregate_trace_by_process(tasks)
            resource_usage_html = generate_resource_usage_html(tasks, processes)
            print(f"Processed {len(tasks)} tasks from trace file")

    # Load software versions if available
    versions_data = None
    if args.versions and args.versions.exists():
        print(f"Loading software versions...")
        try:
            with open(args.versions, 'r') as f:
                versions_data = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load versions file: {e}")

    # Load GCF biosynthetic tree image as base64 for embedding in report
    # Prefer SVG (vector, publication quality) over PNG when available
    import base64
    gcf_tree_b64 = None
    gcf_tree_mime = 'image/png'
    if args.gcf_tree_svg and args.gcf_tree_svg.exists():
        with open(args.gcf_tree_svg, 'rb') as f:
            gcf_tree_b64 = base64.b64encode(f.read()).decode('ascii')
        gcf_tree_mime = 'image/svg+xml'
    elif args.gcf_tree and args.gcf_tree.exists():
        with open(args.gcf_tree, 'rb') as f:
            gcf_tree_b64 = base64.b64encode(f.read()).decode('ascii')

    # Load GCF × species heatmap SVG as base64
    gcf_heatmap_b64 = None
    if args.gcf_heatmap_svg and args.gcf_heatmap_svg.exists():
        with open(args.gcf_heatmap_svg, 'rb') as f:
            gcf_heatmap_b64 = base64.b64encode(f.read()).decode('ascii')

    # Load all-BGCs circular tree as base64. SVG first: it is vector (this figure has
    # ~320 leaves and is unreadable without zoom) and much smaller once base64-encoded
    # — 1.8 MB PNG vs 972 KB SVG on the Pantoea genus run.
    all_bgcs_tree_b64 = None
    all_bgcs_tree_mime = 'image/png'
    _all_bgcs_svg = args.all_bgcs_tree_svg
    if not (_all_bgcs_svg and _all_bgcs_svg.exists()) and args.all_bgcs_tree:
        # fall back to an SVG sitting beside the PNG (kept for standalone invocation;
        # under Nextflow only declared inputs are staged, so the sibling is usually absent)
        sibling = args.all_bgcs_tree.with_suffix('.svg')
        _all_bgcs_svg = sibling if sibling.exists() else None
    if _all_bgcs_svg and _all_bgcs_svg.exists():
        with open(_all_bgcs_svg, 'rb') as f:
            all_bgcs_tree_b64 = base64.b64encode(f.read()).decode('ascii')
        all_bgcs_tree_mime = 'image/svg+xml'
    elif args.all_bgcs_tree and args.all_bgcs_tree.exists():
        with open(args.all_bgcs_tree, 'rb') as f:
            all_bgcs_tree_b64 = base64.b64encode(f.read()).decode('ascii')

    # Build coupling enzyme table rows from live BiG-SCAPE data
    coupling_table_rows = None
    gcf_classes = None
    gcf_support_rows = None
    bigscape_db_for_coupling = None
    if args.bigscape_db and args.bigscape_db.exists():
        bigscape_db_for_coupling = args.bigscape_db
    elif args.bigscape_stats and args.bigscape_stats.exists():
        import glob as _glob
        db_candidates = _glob.glob(str(args.bigscape_stats.parent / '*.db'))
        if db_candidates:
            bigscape_db_for_coupling = Path(db_candidates[0])
    if (args.coupling_annotation and args.coupling_annotation.exists()
            and bigscape_db_for_coupling):
        print("Building coupling enzyme table from live data...")
        coupling_table_rows = build_coupling_table_rows(
            args.coupling_annotation, bigscape_db_for_coupling)
        # Same map drives the GCF badge colours in the Novel BGCs table
        gcf_classes = gcf_coupling_classes(
            args.coupling_annotation, bigscape_db_for_coupling)
        if args.coupling_support and args.coupling_support.exists():
            gcf_support_rows = build_gcf_support_rows(
                args.coupling_support, args.coupling_annotation, bigscape_db_for_coupling)

    if args.counts or args.tabulation:
        print(f"Generating HTML report...")
        generate_html_report(args.outdir, args.taxon, table_header, table_rows, stats, tree_html,
                            bigscape_stats_html, gcf_visualization_html,
                            phylo_tree_generated,
                            genome_table, resource_usage_html, phylo_tree_data,
                            gcf_data=gcf_data_dict, taxonomy_map=taxonomy_map_dict,
                            versions_data=versions_data,
                            rarefaction_stats=rarefaction_stats,
                            gtdbtk_summary_path=gtdbtk_summary_path,
                            gcf_tree_b64=gcf_tree_b64,
                            gcf_tree_mime=gcf_tree_mime,
                            all_bgcs_tree_b64=all_bgcs_tree_b64,
                            all_bgcs_tree_mime=all_bgcs_tree_mime,
                            gcf_heatmap_b64=gcf_heatmap_b64,
                            coupling_table_rows=coupling_table_rows,
                            gcf_classes=gcf_classes,
                            gcf_support_rows=gcf_support_rows)
        print(f"Visualizations complete! Open {args.outdir}/bgc_report.html in a browser.")

if __name__ == '__main__':
    main()
