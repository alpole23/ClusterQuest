process AGGREGATE_TAXONOMY {
    tag "Aggregating BGC statistics by taxonomy"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: 'copy'

    input:
    val taxon
    path taxonomy_map
    path region_counts
    path name_map

    output:
    path "taxonomy_tree.json", emit: taxonomy_tree

    script:
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['taxonomy/aggregate_taxonomy.py'])}
    python ${projectDir}/scripts/taxonomy/aggregate_taxonomy.py ${taxonomy_map} ${region_counts} ${name_map} taxonomy_tree.json --taxon "${taxon}"
    """
}
