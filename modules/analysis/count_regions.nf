process COUNT_REGIONS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path "antismash_results/*"

    output:
    path "region_counts.tsv", emit: counts

    script:
    def by_contig = params.count_per_contig ? "--by_contig" : ""
    def split_hybrids = params.split_hybrids ? "--split_hybrids" : ""
    """
    # Cache key. The scripts below are interpolated paths, not declared inputs,
    # so Nextflow would not otherwise notice when they change. Listed explicitly
    # rather than hashing all of scripts/ — see Utils.scriptsHash.
    # scripts-version: ${Utils.scriptsHash(projectDir, ['analysis/count_regions.py'])}
    python ${projectDir}/scripts/analysis/count_regions.py antismash_results region_counts.tsv ${by_contig} ${split_hybrids}
    """
}
