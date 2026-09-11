process EXTRACT_GCF_REPRESENTATIVES {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path bigscape_dir
    path "antismash_input/*", stageAs: 'antismash_input/*'
    path tabulation_file

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "gcf_representatives.json", emit: gcf_data

    script:
    def taxon_clean = Utils.sanitizeTaxon(params.taxon)
    def tabulation_arg = Utils.optArg('--tabulation', tabulation_file)
    """
    python ${projectDir}/scripts/clustering/extract_gcf_representatives.py ${bigscape_dir} antismash_input gcf_representatives.json --taxon "${taxon_clean}" ${tabulation_arg}
    """
}
