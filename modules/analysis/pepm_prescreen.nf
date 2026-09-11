/**
 * Build the diamond database of pepM references once per run.
 *
 * The references ship with the repo (assets/reference_sequences/reference_pepM.faa):
 * seven characterised PEP mutases, four of them MIBiG entries. Expanding the set
 * with a further eight pepMs mined from MIBiG by HMM was tested and **made the
 * screen worse** — no sensitivity gain, 67 more false positives — so seven it is.
 */
process PEPM_MAKEDB {
    label 'process_low'

    input:
    path references

    output:
    path "pepm.dmnd", emit: db

    script:
    """
    diamond makedb --in ${references} -d pepm --quiet
    """
}

/**
 * Drop genomes that cannot contain a phosphonate BGC, before antiSMASH sees them.
 *
 * antiSMASH costs 41.4 CPU-s per genome; this costs ~0.9. On the measured
 * Erwiniaceae set it keeps 11.0% of genomes and loses none of the 298 that carry a
 * BGC, which is what turns an order-scale run from months into weeks.
 *
 * Batched for the same reason antiSMASH is: at one task per genome the scheduler's
 * submission rate, not the work, sets the wall time.
 */
process PEPM_PRESCREEN {
    tag "${genomes instanceof List ? genomes.size() + ' genomes' : genomes.baseName}"
    label 'process_medium'
    publishDir "${params.outdir}/prescreen_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'prescreen_*.tsv'

    input:
    val taxon
    path genomes
    path pepm_db

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "prescreen_${task.index}.tsv", emit: report

    script:
    """
    python ${projectDir}/scripts/analysis/pepm_prescreen.py \\
        --genomes ${genomes} \\
        --db ${pepm_db} \\
        --out prescreen_${task.index}.tsv \\
        --bitscore ${params.pepm_prescreen_bitscore} \\
        --min_density ${params.pepm_prescreen_min_density} \\
        --threads ${task.cpus} \\
        --diamond \$(which diamond)
    """
}
