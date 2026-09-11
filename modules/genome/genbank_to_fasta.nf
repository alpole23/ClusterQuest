/**
 * Convert a batch of GenBank files to FASTA.
 *
 * Batched (params.task_batch_size) because the per-genome work is ~1.5 s.
 * Conversion failures are per-genome: a bad genome is skipped, not fatal for the batch.
 */
process GENBANK_TO_FASTA {
    tag "${genomes instanceof List ? genomes.size() : 1} genomes"
    label 'process_low'
    label 'tolerant'

    input:
    path genomes

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "*.fna", emit: fasta

    script:
    """
    python ${projectDir}/scripts/genome/genbank_to_fasta.py ${genomes}
    """
}
