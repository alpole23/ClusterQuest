include { COUNT_REGIONS } from '../modules/analysis/count_regions'
include { TABULATE_REGIONS } from '../modules/analysis/tabulate_regions'
include { AGGREGATE_TAXONOMY } from '../modules/analysis/aggregate_taxonomy'
include { COUPLING_ENZYME_TREE } from '../modules/analysis/coupling_enzyme_tree'
include { VISUALIZE_RESULTS } from '../modules/visualization/visualize_results'
include { GCF_BIOSYNTHETIC_TREE } from '../modules/visualization/gcf_biosynthetic_tree'
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

            trace_file_ch = file("${params.outdir}/pipeline_info/pipeline_trace.tsv").exists()
                ? Channel.value(file("${params.outdir}/pipeline_info/pipeline_trace.tsv"))
                : Channel.value(file('NO_TRACE_FILE'))

            // --- GCF Biosynthetic Tree (runs before visualization so its output can be embedded) ---
            gcf_tree_png_ch         = placeholder('NO_GCF_TREE')
            gcf_tree_svg_ch         = placeholder('NO_GCF_TREE_SVG')
            all_bgcs_tree_ch        = placeholder('NO_ALL_BGCS_TREE')
            gcf_heatmap_svg_ch      = placeholder('NO_GCF_HEATMAP_SVG')
            coupling_annotation_ch  = placeholder('NO_COUPLING_ANNOTATION')
            if (clusteringEnabled("bigscape")) {
                GCF_BIOSYNTHETIC_TREE(
                    taxon,
                    CLUSTERING.out.bigscape_db,
                    antismash_results,
                    PHYLOGENY.out.tree,
                    PHYLOGENY.out.summary
                )
                gcf_tree_png_ch        = GCF_BIOSYNTHETIC_TREE.out.gcf_tree_png.ifEmpty(file('NO_GCF_TREE'))
                gcf_tree_svg_ch        = GCF_BIOSYNTHETIC_TREE.out.gcf_tree_svg.ifEmpty(file('NO_GCF_TREE_SVG'))
                all_bgcs_tree_ch       = GCF_BIOSYNTHETIC_TREE.out.all_bgcs_tree_png.ifEmpty(file('NO_ALL_BGCS_TREE'))
                gcf_heatmap_svg_ch     = GCF_BIOSYNTHETIC_TREE.out.heatmap_svg.ifEmpty(file('NO_GCF_HEATMAP_SVG'))
                coupling_annotation_ch = GCF_BIOSYNTHETIC_TREE.out.coupling_annotation.ifEmpty(file('NO_COUPLING_ANNOTATION'))

                // --- Coupling Enzyme Trees (pepM + per-class, anchored on reference sequences) ---
                // Runs off the metadata + coupling annotation emitted by GCF_BIOSYNTHETIC_TREE;
                // if either is missing the channels stay empty and the process is skipped.
                if (params.run_coupling_tree) {
                    COUPLING_ENZYME_TREE(
                        taxon,
                        antismash_results,
                        GCF_BIOSYNTHETIC_TREE.out.metadata,
                        GCF_BIOSYNTHETIC_TREE.out.coupling_annotation,
                        file("${projectDir}/assets/reference_sequences/reference_pepM.faa"),
                        file("${projectDir}/assets/reference_sequences/reference_coupling_enzymes.faa")
                    )
                }
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
                gcf_heatmap_svg_ch,
                coupling_annotation_ch
            )
        }
}
