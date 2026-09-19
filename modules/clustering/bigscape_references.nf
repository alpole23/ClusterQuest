/**
 * Distance from every BGC to the characterised phosphonate clusters.
 *
 * The references are deliberately NOT in the main clustering run. A reference that
 * falls within the GCF cutoff of a query joins its family, and everything downstream
 * reads the database without knowing which records are the dataset: measured on a
 * 14-genome subset, one reference raised `total_bgcs` from 18 to 19, grew a family
 * from 10 members to 11, and invented a genome named after the reference directory.
 *
 * So this runs BiG-SCAPE a second time over a COPY of the finished database, which
 * leaves the published clustering untouched and makes the leak impossible rather than
 * something ten downstream scripts must each remember to filter. BiG-SCAPE recognises
 * the work already in the database and computes only the pairs involving a reference.
 *
 * Measured on Erwiniaceae (334 BGCs, 8 cores): 9.4 s against 56.7 s for the main
 * clustering step, inside a pipeline run that takes hours. At 2,000 and 4,000
 * replicated BGCs it is 40 s and 134 s, 10-13% of the main run. `--db-only-output`
 * halves it by skipping HTML and tree regeneration nothing here reads, and was
 * verified to leave every reference distance identical.
 */
process BIGSCAPE_REFERENCES {
    tag "$taxon"
    label 'process_high'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: 'reference_*'
    cache 'lenient'  // directory inputs, as in BIGSCAPE

    input:
    val taxon
    path bigscape_db
    path antismash_results
    path pfam_db
    path reference_dir

    // Digest of the reference GenBanks. A val input, not an interpolated path:
    // Nextflow hashes a directory input by its own metadata rather than its contents,
    // so adding a reference would otherwise leave -resume reusing a run without it.
    val reference_version

    // Digest of the Python this process runs. See CLAUDE.md.
    val scripts_version

    output:
    path "reference_distances.tsv", emit: distances
    path "reference_summary.json", emit: summary

    script:
    def cutoffs = params.bigscape_cutoffs ?: "0.30"
    def cutoff = cutoffs.tokenize(',')[0]
    def alignment_mode = params.bigscape_alignment_mode ?: "auto"
    def classify = params.bigscape_classify ? "--classify ${params.bigscape_classify}" : ""
    def include_singletons = params.bigscape_include_singletons ? "--include-singletons" : ""
    def mix = params.bigscape_mix ? "--mix" : ""
    """
    export PFAM_PATH=\$(readlink -f ${pfam_db})/Pfam-A.hmm
    # BiG-SCAPE loads GBKs in the iteration order of a set keyed on a hash STRING, so
    # without this the load order — and with it the A/B orientation of each pair, which
    # its extended comparisons are not symmetric under — changes every run.
    export PYTHONHASHSEED=0

    mkdir -p antismash_input
    for dir in ${antismash_results}; do
        if [ -d "\$dir" ]; then
            ln -s \$(readlink -f "\$dir") antismash_input/
        fi
    done

    # Work on a copy. The published database stays dataset-only.
    cp ${bigscape_db} with_references.db

    bigscape cluster \\
        -i antismash_input \\
        -o reference_pass \\
        --db-path with_references.db \\
        --reference-dir ${reference_dir} \\
        --pfam-path \$PFAM_PATH \\
        --alignment-mode ${alignment_mode} \\
        --gcf-cutoffs ${cutoffs} \\
        ${include_singletons} \\
        ${mix} \\
        ${classify} \\
        --db-only-output \\
        --cores ${task.cpus}

    python ${projectDir}/scripts/clustering/reference_distances.py \\
        --pass-db with_references.db \\
        --main-db ${bigscape_db} \\
        --reference-dir ${reference_dir} \\
        --cutoff ${cutoff} \\
        --distances reference_distances.tsv \\
        --summary reference_summary.json
    """
}
