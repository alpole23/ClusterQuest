/*
 * Batching regression test for the short per-genome processes.
 *
 * Checks that GENBANK_TO_FASTA and COPY_ANTISMASH_RESULT still:
 *   - emit one file per genome after .flatten()
 *   - tolerate a single unparseable genome without losing the rest of the batch
 *   - copy reused result directories faithfully (hidden + nested files)
 *
 * Run via tests/run_tests.sh — it builds the fixtures and asserts on the output.
 */

nextflow.enable.dsl=2

include { GENBANK_TO_FASTA }      from '../modules/genome/genbank_to_fasta'
include { COPY_ANTISMASH_RESULT } from '../modules/analysis/check_antismash_reuse'
include { batchSize }             from '../subworkflows/helpers'

workflow {
    // Renamed genomes come from the fixtures directly. They used to come from
    // RENAME_GENOMES, which was folded into FETCH_RENAME_SCREEN when the download
    // was batched -- that process fetches from NCBI, so it cannot run offline here.
    // The assembly-ID pairing this once guarded is no longer order-dependent:
    // FETCH_RENAME_SCREEN reads each accession from its own download directory
    // name rather than from staging order.
    renamed = Channel.fromPath("${params.fixtures}/renamed/*.gbff")
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
