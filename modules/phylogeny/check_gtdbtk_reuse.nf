/**
 * Check if GTDB-Tk results can be reused from a previous taxon run.
 * For GTDB-Tk, we check if ALL genomes in the current run exist in the reuse results.
 * If yes, we filter the summary and prune the tree.
 * If no, we fall back to running GTDB-Tk normally.
 */
process CHECK_GTDBTK_REUSE {
    label 'process_local'
    executor 'local'
    cache false  // Always check filesystem

    input:
    val taxon
    val reuse_taxon
    path genome_list  // File with list of genome names (one per line)

    output:
    tuple env('STATUS'), env('REUSE_SUMMARY'), emit: check_result

    script:
    // Use Utils helper to build absolute reuse path
    def reuse_dir = Utils.buildReusePath(params, projectDir, "gtdbtk", reuse_taxon, "gtdbtk_output")
    def reuse_summary = "${reuse_dir}/gtdbtk.bac120.summary.tsv"
    """
    STATUS="RUN"
    REUSE_SUMMARY=""

    # Check if reuse summary exists
    if [ ! -f "${reuse_summary}" ]; then
        echo "RUN: Reuse summary not found at ${reuse_summary}"
        exit 0
    fi

    # Get list of genomes in reuse results (column 1, skip header)
    # GTDB-Tk adds 'usr_' prefix to user genomes
    tail -n +2 "${reuse_summary}" | cut -f1 | sed 's/^usr_//' | sort -u > reuse_genomes.txt

    # Get list of current genomes (from input file, extract basenames without extension)
    cat ${genome_list} | xargs -I {} basename {} .fna | sort -u > current_genomes.txt

    # Check if all current genomes exist in reuse results
    MISSING=\$(comm -23 current_genomes.txt reuse_genomes.txt | wc -l)

    if [ "\$MISSING" -eq 0 ]; then
        STATUS="REUSE"
        REUSE_SUMMARY="${reuse_summary}"
        CURRENT_COUNT=\$(wc -l < current_genomes.txt)
        REUSE_COUNT=\$(wc -l < reuse_genomes.txt)
        echo "REUSE: All \$CURRENT_COUNT genomes found in reuse results (\$REUSE_COUNT total in source)"
    else
        CURRENT_COUNT=\$(wc -l < current_genomes.txt)
        echo "RUN: \$MISSING of \$CURRENT_COUNT genomes not found in reuse results"
        echo "Missing genomes:"
        comm -23 current_genomes.txt reuse_genomes.txt | head -10
    fi
    """
}

/**
 * Filter GTDB-Tk results to only include genomes from current run.
 * Creates filtered summary TSV and pruned phylogenetic tree.
 */
process FILTER_GTDBTK_RESULTS {
    tag "${taxon}"
    label 'process_medium'
    publishDir "${params.outdir}/gtdbtk_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path genome_list      // File with list of genome names
    val reuse_summary     // Path to source summary TSV

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "gtdbtk_output", emit: output_dir
    path "gtdbtk_output/gtdbtk.bac120.summary.tsv", emit: bacterial_summary

    script:
    """
    python ${projectDir}/scripts/phylogeny/filter_gtdbtk_results.py \\
        ${genome_list} \\
        "${reuse_summary}" \\
        gtdbtk_output
    """
}
