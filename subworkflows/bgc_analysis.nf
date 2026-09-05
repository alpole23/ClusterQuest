include { COUNT_REGIONS } from '../modules/analysis/count_regions'
include { TABULATE_REGIONS } from '../modules/analysis/tabulate_regions'
include { AGGREGATE_TAXONOMY } from '../modules/analysis/aggregate_taxonomy'
include { VISUALIZE_RESULTS } from '../modules/visualization/visualize_results'
include { GCF_BIOSYNTHETIC_TREE } from '../modules/visualization/gcf_biosynthetic_tree'
include { PEPM_ALL_BY_ALL } from '../modules/analysis/pepm_all_by_all'
include { PARTITION_TREES } from '../modules/clustering/partition_trees'
include { COLLECT_VERSIONS } from '../modules/utilities/collect_versions'

include { ANTISMASH_ANALYSIS } from './antismash_analysis'
include { CLUSTERING } from './clustering'
include { PHYLOGENY } from './phylogeny'
include { clusteringEnabled; placeholder } from './helpers'

/*
 * Workflow: Complete BGC detection and analysis pipeline
 */
workflow BGC_ANALYSIS {
    take:
        taxon
        renamed_genomes
        assembly_info
        name_map
        taxonomy_map

    main:
        // --- BGC Detection ---
        ANTISMASH_ANALYSIS(taxon, renamed_genomes)
        antismash_results = ANTISMASH_ANALYSIS.out.results

        // --- Region Analysis ---
        counts_ch = placeholder('NO_COUNTS')
        tabulation_ch = placeholder('NO_TABULATION')
        taxonomy_tree_ch = placeholder('NO_TAXONOMY_TREE')

        if (params.run_analysis) {
            COUNT_REGIONS(taxon, antismash_results)
            counts_ch = COUNT_REGIONS.out.counts

            AGGREGATE_TAXONOMY(taxon, taxonomy_map, COUNT_REGIONS.out.counts, name_map)
            taxonomy_tree_ch = AGGREGATE_TAXONOMY.out.taxonomy_tree

            TABULATE_REGIONS(taxon, antismash_results)
            tabulation_ch = TABULATE_REGIONS.out.tabulation
        }

        // --- Clustering ---
        CLUSTERING(taxon, antismash_results, taxonomy_map, name_map, tabulation_ch)

        // --- Phylogenetic Analysis ---
        PHYLOGENY(taxon, renamed_genomes, counts_ch)

        // --- Visualization ---
        if (params.run_analysis) {
            COLLECT_VERSIONS(antismash_results, CLUSTERING.out.bigscape_db, PHYLOGENY.out.summary)
            versions_ch = COLLECT_VERSIONS.out.versions

            trace_file_ch = file("${params.pipeline_info_dir}/pipeline_trace.tsv").exists()
                ? Channel.value(file("${params.pipeline_info_dir}/pipeline_trace.tsv"))
                : Channel.value(file('NO_TRACE_FILE'))

            // --- GCF Biosynthetic Tree (runs before visualization so its output can be embedded) ---
            gcf_tree_png_ch         = placeholder('NO_GCF_TREE')
            gcf_tree_svg_ch         = placeholder('NO_GCF_TREE_SVG')
            all_bgcs_tree_ch        = placeholder('NO_ALL_BGCS_TREE')
            all_bgcs_tree_svg_ch    = placeholder('NO_ALL_BGCS_TREE_SVG')
            gcf_heatmap_svg_ch      = placeholder('NO_GCF_HEATMAP_SVG')
            coupling_annotation_ch  = placeholder('NO_COUPLING_ANNOTATION')
            coupling_support_ch     = placeholder('NO_COUPLING_SUPPORT')
            pepm_svg_ch             = placeholder('NO_PEPM_SVG')
            pepm_json_ch            = placeholder('NO_PEPM_JSON')
            if (clusteringEnabled("bigscape")) {
                GCF_BIOSYNTHETIC_TREE(
                    taxon,
                    CLUSTERING.out.bigscape_db,
                    antismash_results,
                    PHYLOGENY.out.tree,
                    PHYLOGENY.out.summary,
                    CLUSTERING.out.centers_db
                )
                // pepM all-by-all: reproduces Yu et al. 2013 Fig. 2B on this run's
                // data and reports whether pepM identity could partition
                // BiG-SCAPE. Independent of the tree above, so Nextflow runs
                // them concurrently.
                PEPM_ALL_BY_ALL(
                    taxon,
                    CLUSTERING.out.bigscape_db,
                    CLUSTERING.out.pfam_db
                )
                // The figure goes in GCF Analysis, the partitioning table in
                // Pipeline Info. Both optional: PEPM_ALL_BY_ALL emits nothing
                // when there are too few pepMs to compare.
                pepm_svg_ch  = PEPM_ALL_BY_ALL.out.svgs
                    .flatten().filter { it.name.contains('bigscape_similarity') }
                    .ifEmpty(file('NO_PEPM_SVG'))
                pepm_json_ch = PEPM_ALL_BY_ALL.out.summary.ifEmpty(file('NO_PEPM_JSON'))

                // Drill-down beside the global centre tree. Each partition's own
                // database has complete within-partition distances, so these
                // trees substitute nothing — unlike a global all-BGCs tree on a
                // partitioned run.
                PARTITION_TREES(
                    taxon,
                    CLUSTERING.out.partition_dbs,
                    GCF_BIOSYNTHETIC_TREE.out.coupling_annotation
                        .ifEmpty(file('NO_COUPLING_ANNOTATION'))
                )

                gcf_tree_png_ch        = GCF_BIOSYNTHETIC_TREE.out.gcf_tree_png.ifEmpty(file('NO_GCF_TREE'))
                gcf_tree_svg_ch        = GCF_BIOSYNTHETIC_TREE.out.gcf_tree_svg.ifEmpty(file('NO_GCF_TREE_SVG'))
                all_bgcs_tree_ch       = GCF_BIOSYNTHETIC_TREE.out.all_bgcs_tree_png.ifEmpty(file('NO_ALL_BGCS_TREE'))
                all_bgcs_tree_svg_ch   = GCF_BIOSYNTHETIC_TREE.out.all_bgcs_tree_svg.ifEmpty(file('NO_ALL_BGCS_TREE_SVG'))
                gcf_heatmap_svg_ch     = GCF_BIOSYNTHETIC_TREE.out.heatmap_svg.ifEmpty(file('NO_GCF_HEATMAP_SVG'))
                coupling_annotation_ch = GCF_BIOSYNTHETIC_TREE.out.coupling_annotation.ifEmpty(file('NO_COUPLING_ANNOTATION'))
                coupling_support_ch    = GCF_BIOSYNTHETIC_TREE.out.coupling_support.ifEmpty(file('NO_COUPLING_SUPPORT'))
            }

            VISUALIZE_RESULTS(
                taxon,
                counts_ch,
                tabulation_ch,
                assembly_info,
                name_map,
                taxonomy_map,
                taxonomy_tree_ch,
                CLUSTERING.out.bigscape_stats,
                CLUSTERING.out.bigscape_db,
                CLUSTERING.out.gcf_data,
                PHYLOGENY.out.tree,
                PHYLOGENY.out.summary,
                trace_file_ch,
                versions_ch,
                gcf_tree_png_ch,
                gcf_tree_svg_ch,
                all_bgcs_tree_ch,
                all_bgcs_tree_svg_ch,
                gcf_heatmap_svg_ch,
                coupling_annotation_ch,
                coupling_support_ch,
                pepm_svg_ch,
                pepm_json_ch
            )
        }
}
