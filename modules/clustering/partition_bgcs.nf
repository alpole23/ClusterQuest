/**
 * Assign BGCs to BiG-SCAPE partitions by pepM sequence identity.
 *
 * BiG-SCAPE memory is quadratic in BGC count above ~4,000 (~1.9 TB at 121,000).
 * Splitting on pepM identity first rebuilds the identical GCF network — ARI
 * 1.0000 on three independent sets — with the largest job at ~84 GB.
 *
 * Runs before BIGSCAPE, so it reads pepM out of the region GenBanks via
 * hmmsearch rather than from a clustering database that does not exist yet.
 */
process PARTITION_BGCS {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'bgc_partitions.tsv'

    input:
    val taxon
    path "antismash_input/*", stageAs: 'antismash_input/*'
    path pfam_db

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "bgc_partitions.tsv", emit: partitions

    script:
    // The cap exists to keep one partition inside its own allocation, so
    // derive it from that allocation rather than from a guessed constant.
    def mem_gb = task.memory ? task.memory.toGiga() : 32
    """
    python ${projectDir}/scripts/clustering/partition_bgcs.py \\
        --antismash_dir antismash_input \\
        --pfam ${pfam_db}/Pfam-A.hmm \\
        --out bgc_partitions.tsv \\
        --threshold ${params.bigscape_partition_identity} \\
        --max_partition_size ${params.bigscape_partition_max_size} \\
        --partition_threshold ${params.bigscape_partition_threshold} \\
        --memory_gb ${mem_gb} \\
        --memory_margin ${params.bigscape_partition_memory_margin} \\
        --cpus ${task.cpus} \\
        --hmmfetch \$(which hmmfetch) \\
        --hmmsearch \$(which hmmsearch) \\
        --hmmalign \$(which hmmalign)
    """
}
