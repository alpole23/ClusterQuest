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

    output:
    path "*.fna", emit: fasta

    script:
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['genome/genbank_to_fasta.py'])}
    python ${projectDir}/scripts/genome/genbank_to_fasta.py ${genomes}
    """
}
