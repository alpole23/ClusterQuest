include { DOWNLOAD_PFAM } from '../modules/databases/download_pfam'
include { BIGSCAPE } from '../modules/clustering/bigscape'
include { CLUSTERING_STATS } from '../modules/clustering/clustering_stats'
include { PARTITION_BGCS } from '../modules/clustering/partition_bgcs'
include { BIGSCAPE_PARTITION } from '../modules/clustering/bigscape_partition'
include { MERGE_BIGSCAPE } from '../modules/clustering/merge_bigscape'
include { BIGSCAPE_CENTERS } from '../modules/clustering/bigscape_centers'
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
        centers_db_ch = placeholder('NO_CENTERS_DB')
        partition_dbs_ch = Channel.empty()

        if (clusteringEnabled("bigscape")) {
            DOWNLOAD_PFAM()
            pfam_db_ch = DOWNLOAD_PFAM.out.pfam_db

            if (params.bigscape_partition) {
                // Split by pepM identity so no single BiG-SCAPE job goes
                // quadratic. Verified to rebuild the identical GCF network
                // (ARI 1.0000); see CLAUDE.md. The partitioner falls back to one
                // partition below params.bigscape_partition_min, so enabling this
                // on a small taxon costs only the pepM alignment.
                PARTITION_BGCS(taxon, antismash_results, pfam_db_ch,
                               Utils.scriptsHash(projectDir,
                                   ['clustering/partition_bgcs.py', 'utils']))

                partition_ch = PARTITION_BGCS.out.partitions
                    .splitCsv(header: true, sep: '\t')
                    .map { row -> tuple(row.partition, file(row.gbk)) }
                    .groupTuple()

                BIGSCAPE_PARTITION(taxon, partition_ch, pfam_db_ch,
                                   PARTITION_BGCS.out.partitions)
                MERGE_BIGSCAPE(taxon,
                               BIGSCAPE_PARTITION.out.db.collect(),
                               PARTITION_BGCS.out.partitions,
                               Utils.scriptsHash(projectDir,
                                   ['clustering/merge_bigscape_dbs.py']))
                // Keyed by the partition id in the filename, so the
                // per-partition trees can be labelled and published apart.
                partition_dbs_ch = BIGSCAPE_PARTITION.out.db
                    .map { db -> tuple((db.name =~ /part_(\d+)\.db/)[0][1], db) }
                bigscape_db_ch  = MERGE_BIGSCAPE.out.bigscape_db
                bigscape_dir_ch = MERGE_BIGSCAPE.out.bigscape_dir
            } else {
                BIGSCAPE(taxon, antismash_results, pfam_db_ch)
                bigscape_db_ch  = BIGSCAPE.out.bigscape_db
                bigscape_dir_ch = BIGSCAPE.out.bigscape_dir
            }

            // Only partitioned runs have an incomplete distance table, so only
            // they need centre distances measured separately.
            if (params.bigscape_partition) {
                BIGSCAPE_CENTERS(taxon, bigscape_db_ch, pfam_db_ch,
                                 Utils.scriptsHash(projectDir,
                                     ['clustering/extract_family_centers.py']))
                centers_db_ch = BIGSCAPE_CENTERS.out.centers_db.ifEmpty(file('NO_CENTERS_DB'))
            }

            // Reads the database, so it is identical on both paths.
            CLUSTERING_STATS(taxon, bigscape_db_ch,
                             Utils.scriptsHash(projectDir,
                                 ['clustering/stats_from_db.py']))
            bigscape_stats_ch = CLUSTERING_STATS.out.stats_json

            // Extract GCF representatives (needs tabulation for KCB hit lookup)
            if (tabulation.name != 'NO_TABULATION') {
                EXTRACT_GCF_REPRESENTATIVES(taxon, bigscape_dir_ch, antismash_results,
                                            tabulation,
                                            Utils.scriptsHash(projectDir,
                                                ['clustering/extract_gcf_representatives.py',
                                                 'utils']))
                gcf_data_ch = EXTRACT_GCF_REPRESENTATIVES.out.gcf_data
            }
        }

    emit:
        bigscape_stats = bigscape_stats_ch
        bigscape_db    = bigscape_db_ch
        bigscape_dir   = bigscape_dir_ch
        gcf_data       = gcf_data_ch
        pfam_db        = pfam_db_ch
        centers_db     = centers_db_ch
        partition_dbs  = partition_dbs_ch
}
