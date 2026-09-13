/**
 * Complete each BGC's gene annotation from its relatives in the same GCF.
 *
 * Only 41.5% of CDS in the Erwiniaceae run carry an informative product, and 170
 * of 333 regions carry none at all — those assemblies are GenBank-only, with no
 * functional annotation. Gene-content metrics built on that are reading NCBI
 * annotation pipelines, not biology.
 *
 * GCF members are homologous, and annotation quality across them is very uneven:
 * the typical family has a median member at 0% and a best member near 100%. So
 * orthologues are grouped within each family with diamond and a product found on
 * any member is propagated to the rest, with provenance on every transferred call.
 *
 * Runs after CLUSTERING because it needs family membership, and reads the region
 * GenBanks from antiSMASH because BiG-SCAPE's database stores no product names.
 */
process GCF_ANNOTATION_TRANSFER {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}/annotation_transfer",
        mode: params.publish_mode

    input:
    val taxon
    path bigscape_db
    path antismash_results

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "gcf_annotation_transfer.tsv",  emit: per_cds,   optional: true
    path "gcf_consensus_clusters.tsv",   emit: consensus, optional: true
    path "gcf_annotation_transfer.json", emit: summary,   optional: true

    script:
    def cutoff = (params.bigscape_cutoffs.toString().split(',')[0]).trim()
    """
    python ${projectDir}/scripts/analysis/gcf_annotation_transfer.py \\
        --db ${bigscape_db} \\
        --antismash ${antismash_results} \\
        --outdir . \\
        --cutoff ${cutoff} \\
        --min_identity ${params.transfer_min_identity} \\
        --min_coverage ${params.transfer_min_coverage} \\
        --threads ${task.cpus} \\
        --diamond \$(which diamond)
    """
}
