process EXTRACT_TAXONOMY {
    tag "Extracting taxonomy from NCBI metadata"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: 'copy'

    input:
    val taxon
    path assembly_report
    path taxonomy_report
    path taxdump_dir

    output:
    path "taxonomy_map.json", emit: taxonomy_map

    script:
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['taxonomy/extract_taxonomy.py'])}
    export TAXONKIT_DB='${taxdump_dir}'
    python ${projectDir}/scripts/taxonomy/extract_taxonomy.py ${assembly_report} ${taxdump_dir} taxonomy_map.json
    """
}
