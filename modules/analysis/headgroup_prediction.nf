/**
 * Predict the phosphonate intermediate each GCF makes: 2-AEP, 2-HEP, or neither.
 *
 * 2-AEP and 2-HEP are the two common intermediates in phosphonate biosynthesis and
 * the choice between them is made by a single enzyme acting after Ppd: an
 * aepZ-family (class V) transaminase gives 2-AEP, a reductase gives 2-HEP. Both
 * routes run through phosphonoacetaldehyde, so Ppd is required for either.
 *
 * Evidence is tiered — sequence homology to a characterised reference first, domain
 * family second — and the tier is reported as part of the call. Sequence alone
 * called 18 of 19 Erwiniaceae families "none" because aepZ is the only characterised
 * 2-AEP transaminase in existence and distant orthologues fall under any usable
 * identity floor.
 *
 * Carrier class is separate: the same headgroup goes onto a glycan, a lipid, or
 * nothing depending on machinery the cluster may not even encode.
 */
process HEADGROUP_PREDICTION {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'headgroup_prediction.tsv'

    input:
    val taxon
    path bigscape_db
    path antismash_results
    path headgroup_refs

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "headgroup_prediction.tsv", emit: prediction, optional: true

    script:
    """
    python ${projectDir}/scripts/analysis/headgroup_prediction.py \\
        --antismash ${antismash_results} \\
        --db ${bigscape_db} \\
        --references ${headgroup_refs} \\
        --cutoff ${(params.bigscape_cutoffs.toString().split(',')[0]).trim()} \\
        --threads ${task.cpus} \\
        --diamond \$(which diamond) \\
        --out headgroup_prediction.tsv
    """
}
