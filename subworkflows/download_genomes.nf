include { NCBI_FETCH_METADATA } from '../modules/genome/ncbi_fetch_metadata'
include { FETCH_RENAME_SCREEN } from '../modules/genome/fetch_rename_screen'
include { CREATE_NAME_MAP } from '../modules/genome/create_name_map'
include { DOWNLOAD_TAXONKIT_DB } from '../modules/databases/download_taxonkit_db'
include { PEPM_MAKEDB } from '../modules/analysis/pepm_prescreen'
include { EXTRACT_TAXONOMY } from '../modules/analysis/extract_taxonomy'
include { accessionBatches; placeholder } from './helpers'

/*
 * Subworkflow: Download and prepare genomes from NCBI
 *
 * Two stages, split so peak disk stops tracking taxon size:
 *
 *   NCBI_FETCH_METADATA   resolve the taxon to accessions + metadata, no payload
 *   FETCH_RENAME_SCREEN   per batch: fetch, rename, screen, keep the survivors
 *
 * The old single NCBI_DATASETS_DOWNLOAD held every genome in the taxon in one
 * work directory before anything downstream could run -- 1.3 TB for a
 * 150,000-genome taxon, and unavoidable at that shape, because the screen cannot
 * filter what has not finished downloading. See each module's header.
 */
workflow DOWNLOAD_GENOMES {
    take:
        taxon

    main:
        NCBI_FETCH_METADATA(taxon)
        CREATE_NAME_MAP(taxon, NCBI_FETCH_METADATA.out.assembly_info,
                        Utils.scriptsHash(projectDir, ['genome/create_name_map.py']))
        DOWNLOAD_TAXONKIT_DB()

        pepm_db_ch = placeholder('NO_PEPM_DB')
        if (params.pepm_prescreen) {
            PEPM_MAKEDB(file("${projectDir}/assets/reference_sequences/reference_pepM.faa"))
            pepm_db_ch = PEPM_MAKEDB.out.db
        }

        // One file of accessions per batch. Batches are assigned by hash of the
        // accession, not by position, so adding a genome to the taxon dirties one
        // batch instead of shifting every subsequent one and invalidating the
        // whole download on -resume.
        accession_batches = accessionBatches(NCBI_FETCH_METADATA.out.accessions,
                                             params.download_batch_size)

        FETCH_RENAME_SCREEN(taxon, accession_batches, CREATE_NAME_MAP.out.name_map,
                            pepm_db_ch,
                            Utils.scriptsHash(projectDir,
                                ['genome/rename_genome.py', 'analysis/pepm_prescreen.py']))

        EXTRACT_TAXONOMY(
            taxon,
            NCBI_FETCH_METADATA.out.assembly_data_report,
            NCBI_FETCH_METADATA.out.taxonomy_report,
            DOWNLOAD_TAXONKIT_DB.out.taxdump_dir,
            Utils.scriptsHash(projectDir, ['taxonomy/extract_taxonomy.py'])
        )

    emit:
        // flatten: one list per batch, consumers want one genome each
        renamed_genomes      = FETCH_RENAME_SCREEN.out.renamed_genome.flatten()
        // Non-empty only when the screen ran here; ANTISMASH_ANALYSIS reads the
        // matching `prescreened` flag to know it must not screen a second time.
        prescreen_report     = FETCH_RENAME_SCREEN.out.report.collect().ifEmpty([])
        assembly_info        = NCBI_FETCH_METADATA.out.assembly_info
        name_map             = CREATE_NAME_MAP.out.name_map
        taxonomy_map         = EXTRACT_TAXONOMY.out.taxonomy_map
        assembly_data_report = NCBI_FETCH_METADATA.out.assembly_data_report
        taxonomy_report      = NCBI_FETCH_METADATA.out.taxonomy_report
}
