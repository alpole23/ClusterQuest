/**
 * Rank gene cluster families by how much they warrant laboratory follow-up.
 *
 * The report's existing "Potentially Novel" designation is KnownClusterBlast absence,
 * which on Erwiniaceae is true of **0 of 333 regions** — MIBiG holds only a handful of
 * phosphonate clusters, all from actinomycetes, so an Enterobacterales BGC cannot match
 * one. That flag is structurally constant for this chemistry and orders nothing.
 *
 * This ranks on the axis that does vary: divergence from the characterised coupling
 * enzymes the pipeline already aligns against, discounted by how well evidenced the
 * family is. Both components are published beside the score, because the weights are
 * reasoned rather than fitted — there is no set of leads-that-panned-out to fit against,
 * and a lone number would launder that judgement into something that looks measured.
 */
process NOVELTY_SCORE {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'novelty_ranking.tsv'

    input:
    val taxon
    path gcf_representatives
    path tabulation
    path coupling_support

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "novelty_ranking.tsv", emit: ranking, optional: true

    script:
    """
    python ${projectDir}/scripts/analysis/novelty_score.py \\
        --gcf_representatives ${gcf_representatives} \\
        --tabulation ${tabulation} \\
        --coupling_support ${coupling_support} \\
        --out novelty_ranking.tsv
    """
}
