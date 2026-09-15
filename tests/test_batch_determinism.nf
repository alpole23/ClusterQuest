/*
 * Batch composition must not depend on arrival order.
 *
 * collate() slices a channel as items arrive, and upstream tasks finish in whatever
 * order the scheduler gives them. The batches therefore differ between runs, every
 * batched task hashes differently, and -resume matches none of them. Measured on
 * Erwiniaceae: PEPM_PRESCREEN re-ran in full on three consecutive resumes, dragging
 * antiSMASH with it -- about 2.5 h of redundant work per resume.
 *
 * This feeds the same files in forward and reverse order and asserts the batches come
 * out identical. The FWD/REV pairs must match; the RAW pair must not, which is what
 * makes this a test of the fix rather than of collate().
 */
include { sortedBatches; sortedTupleBatches } from '../subworkflows/helpers'

workflow {
    names = (1..9).collect { "genome_${it}.gbff" }
    fwd = Channel.fromList(names.collect { file("/nonexistent/${it}") })
    rev = Channel.fromList(names.reverse().collect { file("/nonexistent/${it}") })

    sortedBatches(fwd, 3).map { b -> b.collect { it.name }.join(',') }.view { "FWD<${it}>" }
    sortedBatches(rev, 3).map { b -> b.collect { it.name }.join(',') }.view { "REV<${it}>" }

    // the unsorted path, to show the failure it is guarding against
    rev.collate(3).map { b -> b.collect { it.name }.join(',') }.view { "RAW<${it}>" }

    // tuple form, keyed on the first element
    tf = Channel.fromList(names.collect { tuple(it - '.gbff', file("/nonexistent/${it}")) })
    tr = Channel.fromList(names.reverse().collect { tuple(it - '.gbff', file("/nonexistent/${it}")) })
    sortedTupleBatches(tf, 3).map { b -> b.collect { it[0] }.join(',') }.view { "TFWD<${it}>" }
    sortedTupleBatches(tr, 3).map { b -> b.collect { it[0] }.join(',') }.view { "TREV<${it}>" }
}
