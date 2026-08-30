/**
 * Extract statistics from BiG-SCAPE clustering output.
 */
process EXTRACT_CLUSTERING_STATS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path input_dir

    output:
    path "bigscape_statistics.json", emit: stats_json

    script:
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['clustering/extract_bigscape_stats.py'])}
    python ${projectDir}/scripts/clustering/extract_bigscape_stats.py ${input_dir} bigscape_statistics.json
    """
}
