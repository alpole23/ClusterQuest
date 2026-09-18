/**
 * Exact distances between GCF family centres.
 *
 * A partitioned run's merged `distance` table holds only within-partition
 * comparisons, so a global tree would substitute a constant for every
 * cross-partition pair — 92 of 171 centre pairs on Erwiniaceae, which leaves the
 * tree's backbone arbitrary. Re-running BiG-SCAPE over just the centres measures
 * all of them; 19 centres took 16 seconds.
 *
 * This scales where an all-BGCs comparison does not: centre count tracks
 * diversity rather than BGC count (19 for Erwiniaceae, 81 for Streptomyces).
 */
process BIGSCAPE_CENTERS {
    tag "$taxon"
    label 'process_high'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'family_centers.{db,tsv}'

    input:
    val taxon
    path bigscape_db
    path pfam_db

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "family_centers.db",  emit: centers_db, optional: true
    path "family_centers.tsv", emit: manifest,   optional: true

    script:
    """
    python ${projectDir}/scripts/clustering/extract_family_centers.py \\
        --db ${bigscape_db} \\
        --outdir center_input \\
        --manifest family_centers.tsv \\
        --cutoff ${params.bigscape_cutoffs.tokenize(',')[0]}

    # Three centres is the minimum for a meaningful tree; below that emit the
    # manifest alone and let the tree step fall back to the main database.
    n=\$(find center_input -name '*.gbk' | wc -l)
    if [ "\$n" -lt 3 ]; then
        echo "only \$n centres; skipping the centre comparison"
        exit 0
    fi

    export PFAM_PATH=\$(readlink -f ${pfam_db})/Pfam-A.hmm
    export PYTHONHASHSEED=0  # see BIGSCAPE: load order decides pair orientation
    bigscape cluster \\
        -i center_input -o center_out \\
        --pfam-path \$PFAM_PATH \\
        --alignment-mode ${params.bigscape_alignment_mode} \\
        --gcf-cutoffs ${params.bigscape_cutoffs} \\
        --include-singletons --cores ${task.cpus} \\
        ${params.bigscape_classify ? "--classify ${params.bigscape_classify}" : ''}

    cp \$(find center_out -name '*.db' | head -1) family_centers.db
    """
}
