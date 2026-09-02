include { DOWNLOAD_PFAM } from '../modules/databases/download_pfam'
include { BIGSCAPE } from '../modules/clustering/bigscape'
include { EXTRACT_CLUSTERING_STATS } from '../modules/clustering/extract_clustering_stats'
include { EXTRACT_GCF_REPRESENTATIVES } from '../modules/clustering/extract_gcf_representatives'
include { clusteringEnabled; placeholder } from './helpers'

/*
 * Subworkflow: BiG-SCAPE GCF clustering
 */
workflow CLUSTERING {
    take:
        taxon
        antismash_results
        taxonomy_map
        name_map
        tabulation

    main:
        bigscape_stats_ch = placeholder('NO_BIGSCAPE_STATS')
        bigscape_db_ch = placeholder('NO_BIGSCAPE_DB')
        bigscape_dir_ch = placeholder('NO_BIGSCAPE_DIR')
        gcf_data_ch = placeholder('NO_GCF_DATA')
        pfam_db_ch = placeholder('NO_PFAM_DB')

        if (clusteringEnabled("bigscape")) {
            DOWNLOAD_PFAM()
            pfam_db_ch = DOWNLOAD_PFAM.out.pfam_db
            BIGSCAPE(taxon, antismash_results, pfam_db_ch)
            EXTRACT_CLUSTERING_STATS(taxon, BIGSCAPE.out.bigscape_dir)
            bigscape_stats_ch = EXTRACT_CLUSTERING_STATS.out.stats_json
            bigscape_db_ch = BIGSCAPE.out.bigscape_db
            bigscape_dir_ch = BIGSCAPE.out.bigscape_dir

            // Extract GCF representatives (needs tabulation for KCB hit lookup)
            if (tabulation.name != 'NO_TABULATION') {
                EXTRACT_GCF_REPRESENTATIVES(taxon, bigscape_dir_ch, antismash_results, tabulation)
                gcf_data_ch = EXTRACT_GCF_REPRESENTATIVES.out.gcf_data
            }
        }

    emit:
        bigscape_stats = bigscape_stats_ch
        bigscape_db    = bigscape_db_ch
        bigscape_dir   = bigscape_dir_ch
        gcf_data       = gcf_data_ch
        pfam_db        = pfam_db_ch
}
