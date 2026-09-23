include { COUNT_REGIONS } from '../modules/analysis/count_regions'
include { TABULATE_REGIONS } from '../modules/analysis/tabulate_regions'
include { AGGREGATE_TAXONOMY } from '../modules/analysis/aggregate_taxonomy'
include { VISUALIZE_RESULTS } from '../modules/visualization/visualize_results'
include { GCF_BIOSYNTHETIC_TREE } from '../modules/visualization/gcf_biosynthetic_tree'
include { NOVELTY_SCORE } from '../modules/analysis/novelty_score'
include { GCF_ANNOTATION_TRANSFER } from '../modules/analysis/gcf_annotation_transfer'
include { BRANCH_POINT_PREDICTION } from '../modules/analysis/branch_point_prediction'
include { BIOSYNTHETIC_PROFILE } from '../modules/analysis/biosynthetic_profile'
include { PEPM_ALL_BY_ALL } from '../modules/analysis/pepm_all_by_all'
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
            COUNT_REGIONS(taxon, antismash_results,
                          Utils.scriptsHash(projectDir, ['analysis/count_regions.py']))
            counts_ch = COUNT_REGIONS.out.counts

            // Needs the taxonomy map, which a bgc_analysis run pointed at
            // another run's genomes does not have. Skipping loses the report's
            // taxonomy tree and nothing else; failing lost the whole run.
            // main.nf warns when it substitutes the placeholder -- do not warn
            // again from inside a channel closure, where `log` is out of scope
            // and -preview cannot catch the NoSuchVariable because it never
            // executes operators.
            AGGREGATE_TAXONOMY(taxon,
                               taxonomy_map.filter { m -> Utils.isValidInput(m) },
                               COUNT_REGIONS.out.counts, name_map,
                               Utils.scriptsHash(projectDir,
                                   ['taxonomy/aggregate_taxonomy.py']))
            taxonomy_tree_ch = AGGREGATE_TAXONOMY.out.taxonomy_tree
                .ifEmpty(file('NO_TAXONOMY_TREE'))

            TABULATE_REGIONS(taxon, antismash_results,
                             Utils.scriptsHash(projectDir,
                                 ['analysis/tabulate_regions.py', 'utils']))
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
            gcf_heatmap_svg_ch      = placeholder('NO_GCF_HEATMAP_SVG')
            coupling_annotation_ch  = placeholder('NO_COUPLING_ANNOTATION')
            coupling_support_ch     = placeholder('NO_COUPLING_SUPPORT')
            novelty_ch              = placeholder('NO_NOVELTY')
            branch_point_ch            = placeholder('NO_BRANCH_POINT')
            bioprofile_ch              = placeholder('NO_BIOPROFILE')
            consensus_ch            = placeholder('NO_CONSENSUS')
            transfer_summary_ch     = placeholder('NO_TRANSFER_SUMMARY')
            pepm_svg_ch             = placeholder('NO_PEPM_SVG')
            pepm_json_ch            = placeholder('NO_PEPM_JSON')
            if (clusteringEnabled("bigscape")) {
                GCF_BIOSYNTHETIC_TREE(
                    taxon,
                    CLUSTERING.out.bigscape_db,
                    antismash_results,
                    PHYLOGENY.out.summary,
                    PHYLOGENY.out.gtdbtk_db,
                    CLUSTERING.out.centers_db,
                    Utils.scriptsHash(projectDir,
                        ['bgc_coupling_annotation.py', 'bgc_gcf_heatmap.py',
                         'bgc_gcf_tree.py', 'bgc_pfam_tree.py', 'utils'])
                )
                // Annotation transfer needs family membership, so it runs after
                // CLUSTERING. Independent of the tree and the all-by-all, so
                // Nextflow runs all three concurrently.
                if (params.annotation_transfer) {
                    GCF_ANNOTATION_TRANSFER(
                        taxon,
                        CLUSTERING.out.bigscape_db,
                        antismash_results,
                        Utils.scriptsHash(projectDir,
                            ['analysis/gcf_annotation_transfer.py', 'utils'])
                    )
                    consensus_ch = GCF_ANNOTATION_TRANSFER.out.consensus
                        .ifEmpty(file('NO_CONSENSUS'))
                    transfer_summary_ch = GCF_ANNOTATION_TRANSFER.out.summary
                        .ifEmpty(file('NO_TRANSFER_SUMMARY'))
                }

                // Secondary branch point: 2-AEP vs 2-HEP, decided by the enzyme
                // acting on phosphonoacetaldehyde. The coupling class names the
                // fate of phosphonopyruvate; this names the fate of its product.
                // Needs family membership, so it runs after CLUSTERING.
                BRANCH_POINT_PREDICTION(
                    taxon,
                    CLUSTERING.out.bigscape_db,
                    antismash_results,
                    file("${projectDir}/assets/reference_sequences/reference_branch_point_enzymes.faa"),
                    Utils.scriptsHash(projectDir,
                        ['analysis/branch_point_prediction.py', 'utils'])
                )
                branch_point_ch = BRANCH_POINT_PREDICTION.out.prediction

                // Reads the finished clustering database, so it runs beside the
                // branch-point call rather than after it.
                BIOSYNTHETIC_PROFILE(
                    taxon,
                    CLUSTERING.out.bigscape_db,
                    Utils.scriptsHash(projectDir,
                        ['analysis/biosynthetic_profile.py', 'utils'])
                )
                bioprofile_ch = BIOSYNTHETIC_PROFILE.out.profile
                    .ifEmpty(file('NO_BIOPROFILE'))

                // pepM all-by-all: reproduces Yu et al. 2013 Fig. 2B on this run's
                // data and reports whether pepM identity could partition
                // BiG-SCAPE. Independent of the tree above, so Nextflow runs
                // them concurrently.
                PEPM_ALL_BY_ALL(
                    taxon,
                    CLUSTERING.out.bigscape_db,
                    CLUSTERING.out.pfam_db,
                    Utils.scriptsHash(projectDir, ['analysis/pepm_all_by_all.py', 'utils'])
                )
                // The figure goes in GCF Analysis, the partitioning table in
                // Pipeline Info. Both optional: PEPM_ALL_BY_ALL emits nothing
                // when there are too few pepMs to compare.
                pepm_svg_ch  = PEPM_ALL_BY_ALL.out.svgs
                    .flatten().filter { it.name.contains('bigscape_similarity') }
                    .ifEmpty(file('NO_PEPM_SVG'))
                pepm_json_ch = PEPM_ALL_BY_ALL.out.summary.ifEmpty(file('NO_PEPM_JSON'))


                gcf_tree_png_ch        = GCF_BIOSYNTHETIC_TREE.out.gcf_tree_png.ifEmpty(file('NO_GCF_TREE'))
                gcf_tree_svg_ch        = GCF_BIOSYNTHETIC_TREE.out.gcf_tree_svg.ifEmpty(file('NO_GCF_TREE_SVG'))
                gcf_heatmap_svg_ch     = GCF_BIOSYNTHETIC_TREE.out.heatmap_svg.ifEmpty(file('NO_GCF_HEATMAP_SVG'))
                coupling_annotation_ch = GCF_BIOSYNTHETIC_TREE.out.coupling_annotation.ifEmpty(file('NO_COUPLING_ANNOTATION'))
                coupling_support_ch    = GCF_BIOSYNTHETIC_TREE.out.coupling_support.ifEmpty(file('NO_COUPLING_SUPPORT'))

                // Needs the coupling support, so it runs after the tree rather than
                // beside the clustering that produced the families.
                NOVELTY_SCORE(
                    taxon,
                    CLUSTERING.out.gcf_data,
                    tabulation_ch,
                    GCF_BIOSYNTHETIC_TREE.out.coupling_support,
                    CLUSTERING.out.bigscape_db,
                    Utils.scriptsHash(projectDir, ['analysis/novelty_score.py'])
                )
                novelty_ch = NOVELTY_SCORE.out.ranking.ifEmpty(file('NO_NOVELTY'))
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
                PHYLOGENY.out.summary,
                trace_file_ch,
                versions_ch,
                gcf_tree_png_ch,
                gcf_tree_svg_ch,
                gcf_heatmap_svg_ch,
                coupling_annotation_ch,
                coupling_support_ch,
                novelty_ch,
                branch_point_ch,
                bioprofile_ch,
                consensus_ch,
                transfer_summary_ch,
                pepm_svg_ch,
                pepm_json_ch,
                Utils.scriptsHash(projectDir, ['visualize_results.py', 'utils', 'viz'])
            )
        }
}
