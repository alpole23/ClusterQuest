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

/**
 * Genomes per pepM pre-screen task.
 *
 * Larger than the antiSMASH batch because the per-genome work is ~45x smaller
 * (0.9 CPU-s against 41.4), so the fixed cost of a task dominates sooner.
 */
def pepmBatchSize() {
    params.pepm_prescreen_batch_size.toString().toInteger()
}

/*
 * Deterministic batching.
 *
 * `collate()` slices a channel in arrival order, and arrival order is not stable
 * between runs: upstream tasks finish in whatever order the scheduler gives them. The
 * batches therefore differ run to run, every batched task gets a different hash, and
 * `-resume` cannot match any of them.
 *
 * Measured on Erwiniaceae: PEPM_PRESCREEN re-ran in full on three consecutive resumes,
 * and because its output gates everything downstream that also re-ran BUILD_PROTEIN_POOL,
 * RECOVER_ORFS and antiSMASH -- about 2.5 hours per resume, for work already done.
 * Nearly every prescreen work directory across those runs holds a different set of
 * genomes; only one pair of the ~24 shared a batch composition.
 *
 * Sorting first makes composition a function of the inputs alone. Sort on the BASENAME,
 * never the path: staged files live under `work/<hash>/`, so the path itself changes
 * every run and sorting by it would be just as unstable.
 */

/** Batches of files, ordered by basename so composition is reproducible. */
def sortedBatches(ch, n) {
    ch.toSortedList { a, b -> a.name <=> b.name }.flatMap { it }.collate(n)
}

/** Batches of tuples keyed on the first element (a genome name), same guarantee.
 *  flatMap rather than flatten: flatten() would tear the tuples apart as well. */
def sortedTupleBatches(ch, n) {
    ch.toSortedList { a, b -> a[0] <=> b[0] }.flatMap { it }.collate(n)
}
