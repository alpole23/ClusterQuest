/**
 * One BGC tree per partition — the drill-down view.
 *
 * Each partition's own database holds *complete* within-partition distances, so
 * these trees have no substituted values at all. They are also small enough to
 * read, which the global all-BGCs tree stops being well before a million
 * genomes (it is already noted as unreadable without zoom at 320 leaves).
 *
 * Pair this with the family-centre tree: the centre tree gives the global
 * relationships between families, these give the structure inside one.
 */
process PARTITION_TREES {
    tag "${taxon} part ${part_id}"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}/partition_trees",
        mode: params.publish_mode

    input:
    val taxon
    tuple val(part_id), path(partition_db)
    path coupling_annotation

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "partition_${part_id}/*", emit: trees, optional: true

    script:
    def coupling_arg = Utils.optArg('--coupling_annotation', coupling_annotation)
    """
    mkdir -p partition_${part_id}

    # A partition of one or two BGCs has no tree to draw; skip rather than fail,
    # since singleton partitions are normal and expected.
    n=\$(python -c "
import sqlite3
c = sqlite3.connect('file:${partition_db}?mode=ro', uri=True)
print(c.execute(\\"SELECT COUNT(*) FROM bgc_record WHERE record_type='region'\\").fetchone()[0])
")
    if [ "\$n" -lt 3 ]; then
        echo "partition ${part_id} has \$n BGCs; no tree"
        exit 0
    fi

    python ${projectDir}/scripts/bgc_all_bgcs_tree.py \\
        --db ${partition_db} \\
        ${coupling_arg} \\
        --outdir partition_${part_id}
    """
}
