/**
 * Merge per-partition BiG-SCAPE databases into the single {taxon}.db that every
 * downstream consumer expects, remapping the per-partition autoincrement ids.
 *
 * The merged `distance` table holds only within-partition comparisons. That is
 * the point of partitioning, not a defect: the omitted pairs are the ones pepM
 * identity established cannot cluster together. Readers must treat a missing
 * pair as "not compared" rather than "distance zero".
 */
process MERGE_BIGSCAPE {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode

    input:
    val taxon
    path partition_dbs
    path partitions

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "${Utils.sanitizeTaxon(params.taxon)}.db", emit: bigscape_db
    // EXTRACT_GCF_REPRESENTATIVES globs *.db out of a directory, so give it one
    // rather than changing a consumer that the unpartitioned path also uses.
    path "merged_dir",                              emit: bigscape_dir
    path "merge_report.txt",                        emit: report

    script:
    """
    python ${projectDir}/scripts/clustering/merge_bigscape_dbs.py \\
        --inputs ${partition_dbs} \\
        --out ${Utils.sanitizeTaxon(params.taxon)}.db \\
        --partitions ${partitions} | tee merge_report.txt

    mkdir -p merged_dir
    cp ${Utils.sanitizeTaxon(params.taxon)}.db merged_dir/
    """
}
