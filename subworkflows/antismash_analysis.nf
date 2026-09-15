include { DOWNLOAD_ANTISMASH_DBS } from '../modules/databases/download_antismash_dbs'
include { GET_ANTISMASH_VERSION; ANTISMASH } from '../modules/analysis/antismash'
include { CHECK_ANTISMASH_REUSE; COPY_ANTISMASH_RESULT } from '../modules/analysis/check_antismash_reuse'
include { PEPM_MAKEDB; PEPM_PRESCREEN } from '../modules/analysis/pepm_prescreen'
include { BUILD_PROTEIN_POOL; RECOVER_ORFS } from '../modules/genome/recover_orfs'
include { batchSize; antismashBatchSize; pepmBatchSize; placeholder; sortedBatches; sortedTupleBatches } from './helpers'

/*
 * Subworkflow: Run antiSMASH on genomes (with optional result reuse)
 */
workflow ANTISMASH_ANALYSIS {
    take:
        taxon
        renamed_genomes

    main:
        DOWNLOAD_ANTISMASH_DBS()

        // Get antiSMASH version for tracking
        GET_ANTISMASH_VERSION()
        antismash_version = GET_ANTISMASH_VERSION.out.version

        // Generate hash of current antiSMASH parameters for tracking
        antismash_params_hash = Utils.antismashParamsHash(params)

        // pepM pre-screen. A genome with no PEP mutase cannot carry a phosphonate
        // BGC, and finding out costs ~0.9 CPU-s against antiSMASH's 41.4. Measured
        // on Erwiniaceae: keeps 11.0% of genomes, loses none of the 298 BGC-positive.
        // On by default since 2026-09-09; set --pepm_prescreen false to send every
        // genome to antiSMASH.
        prescreen_report_ch = Channel.empty()
        if (params.pepm_prescreen) {
            PEPM_MAKEDB(file("${projectDir}/assets/reference_sequences/reference_pepM.faa"))
            PEPM_PRESCREEN(
                taxon,
                sortedBatches(renamed_genomes, pepmBatchSize()),
                PEPM_MAKEDB.out.db,
                Utils.scriptsHash(projectDir, ['analysis/pepm_prescreen.py'])
            )
            prescreen_report_ch = PEPM_PRESCREEN.out.report

            // Same join pattern PHYLOGENY uses to restrict GTDB-Tk to BGC-positive
            // genomes: key both sides on the bare genome name and inner-join.
            pepm_pass_ch = PEPM_PRESCREEN.out.report
                .splitCsv(header: true, sep: '\t')
                .filter { row -> row.pass == 'yes' }
                .map { row -> tuple(row.genome, true) }

            renamed_genomes = renamed_genomes
                .map { g -> tuple(g.baseName, g) }
                .join(pepm_pass_ch)
                .map { name, genome, flag -> genome }
        }

        // Recover genes the submitted annotations left out, as GFF3 for antiSMASH.
        // After the pre-screen deliberately: the pool is then drawn from genomes that
        // carry a phosphonate pathway, which is both cheaper and a better reference
        // set for the genes being recovered.
        //
        // Pairing is by basename. Genomes are <name>.gbff and their GFF3 <name>.gff3,
        // joined before batching so a task's genomes and GFF3s cannot drift apart --
        // collating two channels independently would pair them only by luck of order.
        recovered_gff_ch = placeholder('NO_RECOVERED_GFF')
        if (params.recover_orfs) {
            BUILD_PROTEIN_POOL(
                taxon,
                renamed_genomes.collect(),
                Utils.scriptsHash(projectDir, ['genome/build_protein_pool.py'])
            )
            RECOVER_ORFS(
                taxon,
                sortedBatches(renamed_genomes, pepmBatchSize()),
                BUILD_PROTEIN_POOL.out.pool,
                Utils.scriptsHash(projectDir, ['genome/recover_orfs.py'])
            )
            recovered_gff_ch = RECOVER_ORFS.out.gff.flatten()
        }

        // Batch genomes with their GFF3s together. join() on basename keeps each
        // genome with its own recovered calls; without it a batch could carry one
        // genome's GFF3 next to another's.
        paired_ch = params.recover_orfs
            ? renamed_genomes.map { g -> tuple(g.baseName, g) }
                // NOT simpleName: it strips every extension, so a genome named
                // ..._GCA_963520565.1 would key as ..._GCA_963520565 and fail to
                // join. 652 of 2,771 Erwiniaceae genome names contain a dot; they
                // would have been dropped from antiSMASH silently.
                .join(recovered_gff_ch.map { f -> tuple(f.name.replaceAll(/\.gff3$/, ''), f) })
            : renamed_genomes.map { g -> tuple(g.baseName, g, file('NO_RECOVERED_GFF')) }

        if (params.reuse_antismash_from) {
            // === REUSE MODE ===
            CHECK_ANTISMASH_REUSE(
                taxon,
                params.reuse_antismash_from,
                renamed_genomes,
                antismash_version,
                antismash_params_hash
            )

            // Split: genomes to run vs genomes to reuse
            genomes_to_run = CHECK_ANTISMASH_REUSE.out.check_result
                .filter { genome, status, path -> status == "RUN" }
                .map { genome, status, path -> genome }

            genomes_to_reuse = CHECK_ANTISMASH_REUSE.out.check_result
                .filter { genome, status, path -> status == "REUSE" }
                .map { genome, status, path -> tuple(genome.baseName, file(path)) }

            // Run antiSMASH on genomes that need it
            run_batches = genomes_to_run.map { g -> tuple(g.baseName, g) }
                .join(paired_ch.map { n, g, f -> tuple(n, f) })
            run_batches = sortedTupleBatches(run_batches, antismashBatchSize())
            ANTISMASH(taxon,
                      run_batches.map { rows -> rows.collect { it[1] } },
                      run_batches.map { rows -> rows.collect { it[2] } },
                      DOWNLOAD_ANTISMASH_DBS.out.db_dir, antismash_version, antismash_params_hash)

            // Copy reused results in batches (each copy is ~1 s — one job per genome
            // is almost entirely scheduler overhead)
            COPY_ANTISMASH_RESULT(taxon, sortedTupleBatches(genomes_to_reuse, batchSize()))

            // Combine all results
            antismash_results = ANTISMASH.out.result_dir.flatten()
                .mix(COPY_ANTISMASH_RESULT.out.result_dir.flatten())
                .collect()
        } else {
            // === NORMAL MODE ===
            batches = sortedTupleBatches(paired_ch, antismashBatchSize())
            ANTISMASH(taxon,
                      batches.map { rows -> rows.collect { it[1] } },
                      batches.map { rows -> rows.collect { it[2] } },
                      DOWNLOAD_ANTISMASH_DBS.out.db_dir, antismash_version, antismash_params_hash)
            antismash_results = ANTISMASH.out.result_dir.flatten().collect()
        }

    emit:
        results          = antismash_results
        prescreen_report = prescreen_report_ch
        recovered_gff    = recovered_gff_ch
}
