process VISUALIZE_RESULTS {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: 'copy'

    input:
    val taxon
    path counts_file
    path tabulation_file
    path assembly_info
    path name_map
    path taxonomy_map
    path taxonomy_tree
    path bigscape_stats
    path bigscape_db
    path gcf_data
    path phylo_tree
    path gtdbtk_summary
    path trace_file
    path versions_file
    path gcf_tree_png
    path gcf_tree_svg
    path all_bgcs_tree_png
    path all_bgcs_tree_svg
    path gcf_heatmap_svg
    path coupling_annotation
    path coupling_support

    output:
    path "*.png", emit: plots, optional: true
    path "*.html", emit: reports, optional: true
    path "*.nwk", emit: newick_files, optional: true
    path "genomes/*.html", emit: genome_pages, optional: true

    script:

    def counts_arg              = Utils.optArg('--counts',              counts_file)
    def tab_arg                 = Utils.optArg('--tabulation',          tabulation_file)
    def assembly_arg            = Utils.optArg('--assembly_info',       assembly_info)
    def name_map_arg            = Utils.optArg('--name_map',            name_map)
    def taxonomy_map_arg        = Utils.optArg('--taxonomy_map',        taxonomy_map)
    def taxonomy_tree_arg       = Utils.optArg('--taxonomy_tree',       taxonomy_tree)
    def bigscape_stats_arg      = Utils.optArg('--bigscape_stats',      bigscape_stats)
    def bigscape_db_arg         = Utils.optArg('--bigscape_db',         bigscape_db)
    def gcf_data_arg            = Utils.optArg('--gcf_data',            gcf_data)
    def phylo_tree_arg          = Utils.optArg('--phylo_tree',          phylo_tree)
    def gtdbtk_summary_arg      = Utils.optArg('--gtdbtk_summary',      gtdbtk_summary)
    def trace_arg               = Utils.optArg('--trace',               trace_file)
    def versions_arg            = Utils.optArg('--versions',            versions_file)
    def gcf_tree_arg            = Utils.optArg('--gcf_tree',            gcf_tree_png)
    def gcf_tree_svg_arg        = Utils.optArg('--gcf_tree_svg',        gcf_tree_svg)
    def all_bgcs_tree_arg       = Utils.optArg('--all_bgcs_tree',       all_bgcs_tree_png)
    def all_bgcs_tree_svg_arg   = Utils.optArg('--all_bgcs_tree_svg',   all_bgcs_tree_svg)
    def gcf_heatmap_svg_arg     = Utils.optArg('--gcf_heatmap_svg',     gcf_heatmap_svg)
    def coupling_annotation_arg = Utils.optArg('--coupling_annotation', coupling_annotation)
    def coupling_support_arg    = Utils.optArg('--coupling_support',    coupling_support)

    def mibig_arg     = params.bigscape_mibig_version ? "--mibig_included" : ""
    def skip_tree_arg = params.skip_tree ? "--skip_tree" : ""
    def outgroup_arg  = params.gtdbtk_outgroup ? "--outgroup '${params.gtdbtk_outgroup}'" : ""
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['visualize_results.py', 'utils', 'viz'])}
    python ${projectDir}/scripts/visualize_results.py ${counts_arg} ${tab_arg} ${assembly_arg} ${name_map_arg} ${taxonomy_map_arg} ${taxonomy_tree_arg} ${bigscape_stats_arg} ${bigscape_db_arg} ${gcf_data_arg} ${phylo_tree_arg} ${gtdbtk_summary_arg} ${trace_arg} ${versions_arg} ${mibig_arg} ${skip_tree_arg} ${outgroup_arg} ${gcf_tree_arg} ${gcf_tree_svg_arg} ${all_bgcs_tree_arg} ${all_bgcs_tree_svg_arg} ${gcf_heatmap_svg_arg} ${coupling_annotation_arg} ${coupling_support_arg} --outdir . --taxon "${taxon}"
    """
}
