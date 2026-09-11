include { DOWNLOAD_ANTISMASH_DBS } from '../modules/databases/download_antismash_dbs'
include { GET_ANTISMASH_VERSION; ANTISMASH } from '../modules/analysis/antismash'
include { CHECK_ANTISMASH_REUSE; COPY_ANTISMASH_RESULT } from '../modules/analysis/check_antismash_reuse'
include { PEPM_MAKEDB; PEPM_PRESCREEN } from '../modules/analysis/pepm_prescreen'
include { batchSize; antismashBatchSize; pepmBatchSize } from './helpers'

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
        // Off by default because it *removes* genomes from the analysis.
        prescreen_report_ch = Channel.empty()
        if (params.pepm_prescreen) {
            PEPM_MAKEDB(file("${projectDir}/assets/reference_sequences/reference_pepM.faa"))
            PEPM_PRESCREEN(
                taxon,
                renamed_genomes.collate(pepmBatchSize()),
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
            ANTISMASH(taxon, genomes_to_run.collate(antismashBatchSize()),
                      DOWNLOAD_ANTISMASH_DBS.out.db_dir, antismash_version, antismash_params_hash)

            // Copy reused results in batches (each copy is ~1 s — one job per genome
            // is almost entirely scheduler overhead)
            COPY_ANTISMASH_RESULT(taxon, genomes_to_reuse.collate(batchSize()))

            // Combine all results
            antismash_results = ANTISMASH.out.result_dir.flatten()
                .mix(COPY_ANTISMASH_RESULT.out.result_dir.flatten())
                .collect()
        } else {
            // === NORMAL MODE ===
            ANTISMASH(taxon, renamed_genomes.collate(antismashBatchSize()),
                      DOWNLOAD_ANTISMASH_DBS.out.db_dir, antismash_version, antismash_params_hash)
            antismash_results = ANTISMASH.out.result_dir.flatten().collect()
        }

    emit:
        results          = antismash_results
        prescreen_report = prescreen_report_ch
}
