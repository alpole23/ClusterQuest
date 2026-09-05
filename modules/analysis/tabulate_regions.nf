process TABULATE_REGIONS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path "antismash_results/*"

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "region_tabulation.tsv", emit: tabulation

    script:
    """
    python ${projectDir}/scripts/analysis/tabulate_regions.py antismash_results region_tabulation.tsv
    """
}
