/**
 * Rank gene cluster families by how much they warrant laboratory follow-up.
 *
 * The report's existing "Potentially Novel" designation is KnownClusterBlast absence,
 * which on Erwiniaceae is true of **0 of 333 regions** — MIBiG holds only a handful of
 * phosphonate clusters, all from actinomycetes, so an Enterobacterales BGC cannot match
 * one. That flag is structurally constant for this chemistry and orders nothing.
 *
 * Divergence from the characterised references was the first attempt and is itself
 * biased: six of the seven references are Streptomyces, so identities are bimodal with
 * nothing between 45% and 94% — a readout of whether a same-taxon reference exists.
 * This ranks on ISOLATION instead: how far a family sits from every other family in the
 * run, over the all-pairs matrix BiG-SCAPE already computes. Reference identity is kept
 * as reported context and used only to zero a family whose chemistry is characterised.
 *
 * Both components are published beside the score, because the weights are reasoned
 * rather than fitted — there is no set of leads-that-panned-out to fit against, and a
 * lone number would launder that judgement into something that looks measured.
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
    path bigscape_db

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
        --bigscape_db ${bigscape_db} \\
        --cutoff ${(params.bigscape_cutoffs.toString().split(',')[0]).trim()} \\
        --out novelty_ranking.tsv
    """
}
