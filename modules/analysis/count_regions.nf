process COUNT_REGIONS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path "antismash_results/*"

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "region_counts.tsv", emit: counts

    script:
    def by_contig = params.count_per_contig ? "--by_contig" : ""
    def split_hybrids = params.split_hybrids ? "--split_hybrids" : ""
    """
    python ${projectDir}/scripts/analysis/count_regions.py antismash_results region_counts.tsv ${by_contig} ${split_hybrids}
    """
}
