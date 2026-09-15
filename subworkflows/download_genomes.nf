include { NCBI_DATASETS_DOWNLOAD } from '../modules/genome/ncbi_datasets_download'
include { CREATE_NAME_MAP } from '../modules/genome/create_name_map'
include { RENAME_GENOMES } from '../modules/genome/rename_genomes_parallel'
include { DOWNLOAD_TAXONKIT_DB } from '../modules/databases/download_taxonkit_db'
include { EXTRACT_TAXONOMY } from '../modules/analysis/extract_taxonomy'
include { batchSize; sortedTupleBatches } from './helpers'

/*
 * Subworkflow: Download and prepare genomes from NCBI
 */
workflow DOWNLOAD_GENOMES {
    take:
        taxon

    main:
        NCBI_DATASETS_DOWNLOAD(taxon)
        CREATE_NAME_MAP(taxon, NCBI_DATASETS_DOWNLOAD.out.assembly_info,
                        Utils.scriptsHash(projectDir, ['genome/create_name_map.py']))
        DOWNLOAD_TAXONKIT_DB()

        // Prepare genome pairs (assembly_id, genome_file), batched — renaming is a
        // ~0.2 s copy, so one task per genome is pure scheduler overhead
        genome_pairs = NCBI_DATASETS_DOWNLOAD.out.genomes
            .flatten()
            .map { gbff -> tuple(gbff.parent.name, gbff) }
        genome_batches = sortedTupleBatches(genome_pairs, batchSize())
            .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }

        RENAME_GENOMES(taxon, genome_batches, CREATE_NAME_MAP.out.name_map,
                       Utils.scriptsHash(projectDir, ['genome/rename_genome.py']))

        EXTRACT_TAXONOMY(
            taxon,
            NCBI_DATASETS_DOWNLOAD.out.assembly_data_report,
            NCBI_DATASETS_DOWNLOAD.out.taxonomy_report,
            DOWNLOAD_TAXONKIT_DB.out.taxdump_dir,
            Utils.scriptsHash(projectDir, ['taxonomy/extract_taxonomy.py'])
        )

    emit:
        // flatten: RENAME_GENOMES emits one list per batch, consumers want one genome each
        renamed_genomes      = RENAME_GENOMES.out.renamed_genome.flatten()
        assembly_info        = NCBI_DATASETS_DOWNLOAD.out.assembly_info
        name_map             = CREATE_NAME_MAP.out.name_map
        taxonomy_map         = EXTRACT_TAXONOMY.out.taxonomy_map
        assembly_data_report = NCBI_DATASETS_DOWNLOAD.out.assembly_data_report
        taxonomy_report      = NCBI_DATASETS_DOWNLOAD.out.taxonomy_report
}
