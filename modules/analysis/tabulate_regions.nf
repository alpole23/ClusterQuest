process TABULATE_REGIONS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: 'copy'

    input:
    val taxon
    path "antismash_results/*"

    output:
    path "region_tabulation.tsv", emit: tabulation

    script:
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['analysis/tabulate_regions.py', 'utils'])}
    python ${projectDir}/scripts/analysis/tabulate_regions.py antismash_results region_tabulation.tsv
    """
}
