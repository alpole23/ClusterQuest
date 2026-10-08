/**
 * Reference distances for a partitioned run.
 *
 * BIGSCAPE_PARTITION already did the comparing, one reference pass per partition, so
 * nothing here runs BiG-SCAPE. This only reads those databases and joins each measured
 * reference x BGC pair back to the family the MERGED run assigned, which is why it
 * needs the merged database too.
 *
 * The join is on (genome, region), not on record id: merge_bigscape_dbs.py offsets
 * record ids per partition, so a partition's record 7 is not the merged database's
 * record 7. reference_distances.py carries that reasoning in full.
 *
 * Output is identical in shape to BIGSCAPE_REFERENCES, so the report reads one channel
 * on both paths and needs to know nothing about partitioning.
 */
process PARTITION_REFERENCE_DISTANCES {
    tag "$taxon"
    label 'process_low'            // reads SQLite; no alignment, no BiG-SCAPE
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'reference_*'

    input:
    val taxon
    path refpass_dbs
    path merged_db
    path reference_dir

    // Digest of the reference GenBanks; see BIGSCAPE_PARTITION for why this is a val.
    val reference_version

    // Digest of the Python this process runs. See CLAUDE.md.
    val scripts_version

    output:
    path "reference_distances.tsv", emit: distances
    path "reference_summary.json", emit: summary

    script:
    def cutoffs = params.bigscape_cutoffs ?: "0.30"
    def cutoff = cutoffs.tokenize(',')[0]
    """
    python ${projectDir}/scripts/clustering/reference_distances.py \\
        --pass-db ${refpass_dbs} \\
        --main-db ${merged_db} \\
        --reference-dir ${reference_dir} \\
        --cutoff ${cutoff} \\
        --distances reference_distances.tsv \\
        --summary reference_summary.json
    """
}
