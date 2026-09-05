include { GENBANK_TO_FASTA } from '../modules/genome/genbank_to_fasta'
include { DOWNLOAD_GTDBTK_DB } from '../modules/databases/download_gtdbtk_db'
include { GTDBTK_CLASSIFY } from '../modules/phylogeny/gtdbtk'
include { CHECK_GTDBTK_REUSE; FILTER_GTDBTK_RESULTS } from '../modules/phylogeny/check_gtdbtk_reuse'
include { batchSize; placeholder } from './helpers'

/*
 * Subworkflow: GTDB-Tk phylogenetic classification (with optional result reuse)
 */
workflow PHYLOGENY {
    take:
        taxon
        renamed_genomes
        counts

    main:
        phylo_tree_ch = placeholder('NO_PHYLO_TREE')
        gtdbtk_summary_ch = placeholder('NO_GTDBTK_SUMMARY')

        if (params.run_gtdbtk) {
            // Determine which genomes to process
            if (params.gtdbtk_bgc_genomes_only) {
                // Filter to genomes with BGCs
                genomes_with_bgcs_ch = counts
                    .splitCsv(header: true, sep: '\t', skip: 1)
                    .filter { row -> (row.total_count ?: '0').toInteger() > 0 }
                    .map { row -> tuple(row.record, true) }

                renamed_genomes_tuples = renamed_genomes
                    .map { genome -> tuple(genome.name, genome) }

                genomes_for_gtdbtk = renamed_genomes_tuples
                    .join(genomes_with_bgcs_ch)
                    .map { name, genome, flag -> genome }
            } else {
                genomes_for_gtdbtk = renamed_genomes
            }

            // Convert GenBank to FASTA (batched — ~1.5 s per genome)
            GENBANK_TO_FASTA(genomes_for_gtdbtk.collate(batchSize()),
                             Utils.scriptsHash(projectDir,
                                 ['genome/genbank_to_fasta.py']))
            fasta_ch = GENBANK_TO_FASTA.out.fasta.flatten()
            fasta_files = fasta_ch.collect()

            if (params.reuse_gtdbtk_from) {
                // === GTDB-Tk REUSE MODE ===
                genome_list_ch = fasta_ch
                    .map { it.toString() }
                    .collectFile(name: 'genome_list.txt', newLine: true)

                CHECK_GTDBTK_REUSE(taxon, params.reuse_gtdbtk_from, genome_list_ch)

                check_result = CHECK_GTDBTK_REUSE.out.check_result
                    .branch {
                        reuse: it[0] == "REUSE"
                        run: it[0] == "RUN"
                    }

                // REUSE path
                FILTER_GTDBTK_RESULTS(
                    taxon,
                    genome_list_ch,
                    check_result.reuse.map { it[1] },
                    check_result.reuse.map { it[2] },
                    Utils.scriptsHash(projectDir,
                        ['phylogeny/filter_gtdbtk_results.py'])
                )

                // RUN path
                DOWNLOAD_GTDBTK_DB()
                fasta_for_fresh_run = check_result.run
                    .combine(fasta_files)
                    .map { status, summary, tree, files -> files }
                GTDBTK_CLASSIFY(taxon, fasta_for_fresh_run, DOWNLOAD_GTDBTK_DB.out.db_dir)

                // Combine outputs
                phylo_tree_ch = FILTER_GTDBTK_RESULTS.out.bacterial_tree
                    .mix(GTDBTK_CLASSIFY.out.bacterial_tree)
                    .ifEmpty(file('NO_PHYLO_TREE'))
                gtdbtk_summary_ch = FILTER_GTDBTK_RESULTS.out.bacterial_summary
                    .mix(GTDBTK_CLASSIFY.out.bacterial_summary)
                    .ifEmpty(file('NO_GTDBTK_SUMMARY'))
            } else {
                // === GTDB-Tk NORMAL MODE ===
                DOWNLOAD_GTDBTK_DB()
                GTDBTK_CLASSIFY(taxon, fasta_files, DOWNLOAD_GTDBTK_DB.out.db_dir)

                phylo_tree_ch = GTDBTK_CLASSIFY.out.bacterial_tree.ifEmpty(file('NO_PHYLO_TREE'))
                gtdbtk_summary_ch = GTDBTK_CLASSIFY.out.bacterial_summary.ifEmpty(file('NO_GTDBTK_SUMMARY'))
            }
        }

    emit:
        tree    = phylo_tree_ch
        summary = gtdbtk_summary_ch
}
