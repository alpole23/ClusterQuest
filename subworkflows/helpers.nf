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
