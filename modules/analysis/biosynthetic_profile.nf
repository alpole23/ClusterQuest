/**
 * Per-family biosynthetic domain profile, and which family splits are chemistry.
 *
 * BiG-SCAPE compares whole regions, and a region is a rule core plus a symmetric
 * flank. At the 10 kb neighbourhood this pipeline uses that flank reaches ~9 kb
 * past the core in both directions and reliably catches chromosomal neighbours,
 * so two groups of genomes with different neighbours land in different families
 * even when their biosynthesis is identical.
 *
 * This reports each family's domain content filtered to the categories that are
 * about making a molecule -- core, tailoring, lipid, transport -- and flags
 * families whose filtered profiles are indistinguishable. Measured on the
 * Erwiniaceae verification run, the 215-member pantaphos family split 186/29 and
 * differed in ten domains across the whole region, zero inside the rule core,
 * and exactly one once filtered.
 *
 * Cheap by construction: it reads the finished clustering database and needs no
 * external tool, so it adds seconds rather than minutes.
 */
process BIOSYNTHETIC_PROFILE {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'biosynthetic_profile.tsv'

    input:
    val taxon
    path bigscape_db

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "biosynthetic_profile.tsv", emit: profile, optional: true

    script:
    """
    python ${projectDir}/scripts/analysis/biosynthetic_profile.py \\
        --db ${bigscape_db} \\
        --cutoff ${params.bigscape_cutoffs.toString().split(',')[0]} \\
        --out biosynthetic_profile.tsv
    """
}
