/**
 * Concatenate the per-shard GTDB-Tk summaries into one.
 *
 * Taxonomy assignments are one independent row per genome, so sharding is
 * lossless for them and the merge is a header-aware concatenation. This is the
 * output every downstream consumer reads; the per-shard trees are discarded
 * (see GTDBTK_CLASSIFY).
 */
process MERGE_GTDBTK {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/gtdbtk_results/${Utils.sanitizeTaxon(params.taxon)}/gtdbtk_output",
        mode: params.publish_mode

    input:
    val taxon
    path shard_summaries

    output:
    path "gtdbtk.bac120.summary.tsv", emit: bacterial_summary, optional: true

    script:
    """
    FIRST=1
    for f in ${shard_summaries}; do
        if [ "\$FIRST" -eq 1 ]; then
            cat "\$f" > gtdbtk.bac120.summary.tsv
            FIRST=0
        else
            tail -n +2 "\$f" >> gtdbtk.bac120.summary.tsv
        fi
    done

    if [ ! -s gtdbtk.bac120.summary.tsv ]; then
        echo "ERROR: no shard produced a summary"
        exit 1
    fi

    N=\$(tail -n +2 gtdbtk.bac120.summary.tsv | wc -l)
    echo "merged \$(echo ${shard_summaries} | wc -w) shard summaries -> \$N genomes"

    # A genome classified twice means the shards overlapped, which would silently
    # inflate every per-clade count downstream.
    DUPES=\$(tail -n +2 gtdbtk.bac120.summary.tsv | cut -f1 | sort | uniq -d | wc -l)
    if [ "\$DUPES" -gt 0 ]; then
        echo "ERROR: \$DUPES genome(s) appear in more than one shard"
        exit 1
    fi
    """
}
