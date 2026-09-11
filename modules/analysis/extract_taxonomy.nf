process EXTRACT_TAXONOMY {
    tag "Extracting taxonomy from NCBI metadata"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path assembly_report
    path taxonomy_report
    path taxdump_dir

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "taxonomy_map.json", emit: taxonomy_map

    script:
    """
    export TAXONKIT_DB='${taxdump_dir}'
    python ${projectDir}/scripts/taxonomy/extract_taxonomy.py ${assembly_report} ${taxdump_dir} taxonomy_map.json
    """
}
