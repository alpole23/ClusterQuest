process EXTRACT_TAXONOMY {
    tag "Extracting taxonomy from NCBI metadata"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path assembly_report
    // NOT consumed. extract_taxonomy.py takes assembly_report, taxdump_dir and
    // output_file -- it never opens this. It was declared as a plain input, and
    // NCBI_FETCH_METADATA emits it `optional: true`, so when NCBI returns no
    // taxonomy report the channel is EMPTY and this process never fires. That
    // silently deletes taxonomy_map, which VISUALIZE_RESULTS needs, so the run
    // produces no report at all and Nextflow still prints "SUCCESS, failed=0" --
    // a skipped process is not a failed one.
    //
    // That is exactly what happened to the whole-Bacteria pilot: 4,193 genomes,
    // 1,345 BGCs, 1,020 families and no bgc_report.html. NCBI emits the file for
    // a named clade (Enterobacterales got one -- 104 bytes, a single line naming
    // the queried taxon) and not for "Bacteria".
    //
    // Staged through a placeholder so the channel always emits. Kept as an input
    // rather than deleted because the staged file documents what the run queried.
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
