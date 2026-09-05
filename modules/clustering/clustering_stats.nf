/**
 * Clustering statistics from the BiG-SCAPE database.
 *
 * Replaces the TSV-parsing EXTRACT_CLUSTERING_STATS, which needed BiG-SCAPE's
 * output *directory* — partitioned runs have one per partition and merge only
 * the databases. Verified to emit byte-identical JSON on an unpartitioned run,
 * so both paths use this.
 */
process CLUSTERING_STATS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode

    input:
    val taxon
    path bigscape_db

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "bigscape_statistics.json", emit: stats_json

    script:
    """
    python ${projectDir}/scripts/clustering/stats_from_db.py \\
        ${bigscape_db} bigscape_statistics.json
    """
}
