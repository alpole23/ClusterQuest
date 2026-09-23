/*
 * Batching regression test for the short per-genome processes.
 *
 * Checks that RENAME_GENOMES, GENBANK_TO_FASTA and COPY_ANTISMASH_RESULT still:
 *   - pair each genome with the right assembly ID across a batch (staging order)
 *   - emit one file per genome after .flatten()
 *   - tolerate a single unparseable genome without losing the rest of the batch
 *   - copy reused result directories faithfully (hidden + nested files)
 *
 * Run via tests/run_tests.sh — it builds the fixtures and asserts on the output.
 */

nextflow.enable.dsl=2

include { RENAME_GENOMES }        from '../modules/genome/rename_genomes_parallel'
include { GENBANK_TO_FASTA }      from '../modules/genome/genbank_to_fasta'
include { COPY_ANTISMASH_RESULT } from '../modules/analysis/check_antismash_reuse'
include { batchSize }             from '../subworkflows/helpers'

workflow {
    // Mirrors the wiring in subworkflows/download_genomes.nf
    genome_batches = Channel.fromPath("${params.fixtures}/data/*/genomic.gbff")
        .map { gbff -> tuple(gbff.parent.name, gbff) }
        .collate(batchSize())
        .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }

    // NO_PEPM_DB placeholder: the screen inside RENAME_GENOMES needs diamond and a
    // reference database, neither of which the fixture harness has. This test is
    // about batch pairing and per-genome failure isolation, so it exercises the
    // unscreened path; the screen itself is covered by the held-out clade data.
    RENAME_GENOMES(params.taxon, genome_batches, file("${params.fixtures}/name_map.json"),
                   file('NO_PEPM_DB'),
                   Utils.scriptsHash(projectDir,
                       ['genome/rename_genome.py', 'analysis/pepm_prescreen.py']))
    renamed = RENAME_GENOMES.out.renamed_genome.flatten()
    renamed.view { "RENAMED: ${it.name}" }

    // Mirrors subworkflows/phylogeny.nf
    GENBANK_TO_FASTA(renamed.collate(batchSize()),
                     Utils.scriptsHash(projectDir, ['genome/genbank_to_fasta.py']))
    GENBANK_TO_FASTA.out.fasta.flatten().view { "FASTA: ${it.name}" }

    // Mirrors the reuse branch of subworkflows/antismash_analysis.nf
    reuse_batches = Channel.fromPath("${params.fixtures}/reuse/*", type: 'dir')
        .map { d -> tuple(d.name, d) }
        .collate(batchSize())
    COPY_ANTISMASH_RESULT(params.taxon, reuse_batches)
    COPY_ANTISMASH_RESULT.out.result_dir.flatten().view { "COPIED: ${it.name}" }
}
