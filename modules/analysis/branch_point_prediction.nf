/**
 * Predict the phosphonate intermediate each GCF makes: 2-AEP, 2-HEP, or neither.
 *
 * 2-AEP and 2-HEP are the two common intermediates in phosphonate biosynthesis and
 * the choice between them is made by a single enzyme acting after Ppd: an
 * aepZ-family (class V) transaminase gives 2-AEP, a reductase gives 2-HEP. Both
 * routes run through phosphonoacetaldehyde, so Ppd is required for either.
 *
 * Evidence is tiered — sequence homology to a characterised reference first, domain
 * family second — and the tier is reported as part of the call. Scoring is by profile HMM,
 * built from the references at run time. Pairwise identity was tried first and
 * called 18 of 19 Erwiniaceae families "none"; on the same references a profile
 * recovered all 26 Enterobacterial AEP clusters where identity found 2. The thin
 * reference set was real but secondary — the method was the larger
 * identity floor.
 *
 * Carrier class is separate: the same intermediate goes onto a glycan, a lipid, or
 * nothing depending on machinery the cluster may not even encode.
 */
process BRANCH_POINT_PREDICTION {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'branch_point_prediction.tsv'

    input:
    val taxon
    path bigscape_db
    path antismash_results
    path branch_point_refs

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "branch_point_prediction.tsv", emit: prediction, optional: true

    script:
    """
    python ${projectDir}/scripts/analysis/branch_point_prediction.py \\
        --antismash ${antismash_results} \\
        --db ${bigscape_db} \\
        --references ${branch_point_refs} \\
        --cutoff ${(params.bigscape_cutoffs.toString().split(',')[0]).trim()} \\
        --threads ${task.cpus} \\
        --hmmbuild \$(which hmmbuild) \\
        --hmmalign \$(which hmmalign) \\
        --hmmsearch \$(which hmmsearch) \\
        --out branch_point_prediction.tsv
    """
}
