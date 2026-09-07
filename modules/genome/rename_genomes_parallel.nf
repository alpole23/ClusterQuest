/**
 * Rename a batch of genome files using the name map.
 *
 * Genomes are processed in batches (params.task_batch_size) because the per-genome
 * work is ~0.2 s — far below the scheduler overhead of one job per genome.
 * Files are staged as genome1.gbff, genome2.gbff, ... (their NCBI names all collide
 * on "genomic.gbff"); staging order matches the assembly_ids list.
 */
process RENAME_GENOMES {
    tag "${assembly_ids.size()} genomes"
    label 'process_low'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(params.taxon)}/renamed_genomes", mode: params.publish_mode

    input:
    val taxon
    tuple val(assembly_ids), path(genome_files, stageAs: 'genome?.gbff')
    path name_map

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "*.gbff", emit: renamed_genome

    script:
    def staged = [genome_files].flatten().collect { it.name }
    def manifest = [assembly_ids, staged].transpose().collect { id, f -> "${id}\t${f}" }.join('\n')
    """
    cat > manifest.tsv <<'MANIFEST_EOF'
${manifest}
MANIFEST_EOF

    python ${projectDir}/scripts/genome/rename_genome.py manifest.tsv ${name_map}
    """
}
