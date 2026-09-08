/*
 * Small helpers shared by the subworkflows.
 *
 * Nextflow functions are file-scoped, so they are included like processes.
 */

/**
 * Create a placeholder channel for an optional input.
 */
def placeholder(name) {
    Channel.value(file(name))
}

/**
 * Check if a clustering method is enabled.
 */
def clusteringEnabled(method) {
    params.clustering == method
}

/**
 * Genomes per task for the short per-genome steps.
 * Coerced explicitly: params supplied on the command line arrive as strings, and
 * collate() only dispatches on a real Integer.
 */
def batchSize() {
    params.task_batch_size.toString().toInteger()
}

/**
 * Genomes per antiSMASH task.
 *
 * Separate from batchSize() because the trade-off is different: the other batched
 * steps take ~1 s per genome, antiSMASH takes ~41 s. At 1 genome per task a
 * million-genome run is bound by the scheduler's submission rate (20/min = 34.7
 * days) rather than by compute (2.4 days at 200 slots). 50 puts submission at
 * 0.7 days, comfortably under compute, so larger batches buy nothing and only
 * coarsen the retry granularity.
 */
def antismashBatchSize() {
    params.antismash_batch_size.toString().toInteger()
}

/**
 * Genomes per GTDB-Tk task.
 *
 * Runs below this stay a single task, so small analyses are unaffected.
 *
 * How GTDB-Tk scales with genome count is **not established**. The two benchmark
 * runs it was fitted on processed 285 and 298 genomes — 4.6% apart — because
 * gtdbtk_bgc_genomes_only restricts it to BGC-positive genomes. Cost rose 2.8%.
 * Those points cannot separate fixed cost from linear cost, and the two models
 * differ by ~94x at a million genomes (48 CPU-h against 4,500). CLAUDE.md
 * previously read the same numbers as "59% more genomes cost 3% more" by
 * comparing *total* genomes rather than the ones GTDB-Tk actually saw.
 *
 * Sharding is insurance against the linear case. If cost is really fixed-dominated
 * it multiplies a ~11.5 CPU-h overhead by the shard count — a few dollars — which
 * is the cheaper mistake to make.
 */
def gtdbtkShardSize() {
    params.gtdbtk_shard_size.toString().toInteger()
}
