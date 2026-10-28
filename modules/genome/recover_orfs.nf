/**
 * Rebuild genes that a genome's submitted annotation left out.
 *
 * Partial annotation is common and silently destructive. Measured on Erwiniaceae:
 * 327 of 333 BGC regions contain an intergenic gap of 300+ bp, the median region is
 * 78.4% coding where a well-annotated operon is ~90%, and 914 kb inside BGC regions
 * carries no gene call. The worst case is P. ananatis LMG 5342 region 2 -- a
 * lab-confirmed phosphonolipid cluster -- where 9 of 15 genes were never called,
 * including the 2-AEP transaminase that decides the product's headgroup.
 *
 * antiSMASH will not fix this itself: it runs gene finding only on records with ZERO
 * CDS features, so one CDS anywhere on a chromosome disables prodigal entirely. But
 * --genefinding-gff3 merges extra calls into an already-annotated record, which is
 * what these processes feed it.
 *
 * Runs after the pepM pre-screen, so the pool is drawn from the genomes that carry a
 * phosphonate pathway rather than from everything downloaded. That is cheaper and
 * also better: screened genomes are each other's best reference for exactly the genes
 * this is trying to recover.
 */

/**
 * Pool every protein already annotated across the screened genomes.
 *
 * This is the reference set for the homology half of recovery. The premise is that
 * annotation quality is uneven between submissions -- LMG 5342's 2012 deposit is
 * missing aepZ, P. ananatis VY148's 2021 assembly has it -- so a gene absent from one
 * genome is usually present in a sibling.
 *
 * A gather step: every RECOVER_ORFS task waits on this one. It is cheap (concatenating
 * proteins out of a few hundred GenBanks) but it does serialise here.
 *
 * It also builds the DIAMOND database, which RECOVER_ORFS used to build once per GENOME
 * over this same pool -- 52% of that stage's wall time spent rebuilding an identical
 * 7.19 GB file, 200 times per batch. Here it is paid once per run.
 */
process BUILD_PROTEIN_POOL {
    tag "$taxon"
    label 'process_low'

    // Memory is a FUNCTION OF THE INPUT, not a constant, because the pool scales
    // with the run and this pipeline has to size itself on whatever machine it is
    // given. A literal tuned to one workstation either wastes a large node or
    // silently blocks scheduling on a small one.
    //
    // Scaled on INPUT BYTES, not on genome count, because the count does not say
    // how much annotation a genome carries and that is what costs memory. Across
    // the inputs this has been measured on, per-genome protein counts span 769
    // (Enterobacterales) to 4,726 (PGAP reference) -- 6.1x, which no per-genome
    // constant absorbs. It is picked up as a pure over-request on sparse input
    // and an OOM on dense input. File size tracks annotation directly, so it
    // absorbs that spread on its own.
    //
    // Calibrated end to end on the Bacteria GenBank run's OWN retained genomes,
    // running this process's script over a seeded random draw and reading peak
    // RSS as ru_maxrss of the single python process. Do NOT calibrate this off a
    // Nextflow `peak_rss`: that sums the process tree and spans the makedb phase
    // as well, and reading the slope off one gave 510 B/unique, 4.5x the truth.
    //
    //     genomes    input   gross CDS     unique    dedup   peak RSS
    //         200   1.37 GB     423,315    275,563   34.9%     126 MB
    //         800   6.15 GB   1,831,247  1,215,059   33.6%     234 MB
    //       3,200  25.25 GB   7,419,776  4,314,865   41.8%     565 MB
    //
    // which fits RSS = 99 MB + 113.5 BYTES per unique protein, over a 16x range.
    // Deduplication RISES with pool size (34.9% -> 41.8%) because the pool is
    // drawn from screened genomes that are each other's close relatives, so
    // unique proteins per input GB FALLS -- 201,140 -> 170,886 -- and the
    // marginal cost at 3,200 genomes is 113.5 B x 170,886 = 1.94% of input
    // bytes, still falling. 2.5% is declared, above the measured rate at every
    // size tested; the retry below covers the case where that is wrong.
    //
    // Beware the dereplicated case, which is the opposite extreme: the
    // whole-Bacteria REFERENCE pilot is one genome per species, so it dedupes
    // only 2.4% and yields 339,527 unique/GB -- 2.0x this run's rate. 2.5% still
    // covers it (58.37 GB in -> 4 GB declared against 2.4 GB needed), but a
    // broader dereplicated input is the direction that would break it.
    //
    // `diamond makedb` runs after the python exits, so it gets that memory back
    // and is not additive. Measured on the pilot's real 7.40 GB pool: peak
    // 1.52 GB, 20.5% of pool size, 108 s, producing a 7.72 GB database. The pool
    // is 373 bytes per unique protein and the database 389, so makedb stays well
    // under the python phase at every size this declaration covers.
    //
    //       387 genomes,    2.9 GB in ->  3 GB declared, 0.16 GB measured
    //     4,193 genomes,   58.4 GB in ->  4 GB declared, 2.4 GB needed
    //   204,000 genomes, 1,659.0 GB in -> 44 GB declared, 32 GB needed
    //                                   (~283 M unique proteins; makedb ~22 GB)
    //
    // NOTE the 204,000 figure is 6.088% retention (21,260 of 349,204 screened)
    // over 3,350,527 records at the 8.13 MB/genome this run's retained set
    // measures. That pool implies a ~110 GB DIAMOND database, which RECOVER_ORFS
    // below declares as 6 + ceil(db_GB) = ~116 GB -- re-check that directive
    // before a run of this size, not this one.
    //
    // Declared here rather than in conf/labels.config because a config `withName`
    // closure can see `task` but NOT the process inputs, and it would override
    // this if it also set memory.
    //
    // NOTE .size() is COUNT on a List and BYTES on a single file, hence the guard.
    memory {
        def inGB = (genomes instanceof List
                        ? (genomes.sum { it.size() } ?: 0)
                        : genomes.size()) / 1e9
        1.GB * (2 + Math.ceil(inGB * 0.025) * task.attempt)
    }
    maxRetries 1
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140] && task.attempt <= 1 ? 'retry' : 'terminate' }

    input:
    val taxon
    path genomes

    // Digest of the Python this process runs. A val input, not an
    // interpolation: Nextflow hashes the unevaluated script source plus the
    // input values, never the rendered text. See CLAUDE.md.
    val scripts_version

    output:
    path "protein_pool.faa",  emit: pool
    path "protein_pool.dmnd", emit: db

    script:
    """
    python ${projectDir}/scripts/genome/build_protein_pool.py \\
        --genomes ${genomes} \\
        --out protein_pool.faa

    # --threads explicitly, unlike the per-genome call this replaces: diamond defaults
    # to every core on the box, which under Nextflow meant 8 concurrent tasks each
    # claiming 16 cores. Measured on a 6.89 GB pool: 84 s / 268 CPU-s at 4 threads
    # against 59 s / 500 CPU-s at 16 -- the wall saving is not worth 1.9x the CPU.
    diamond makedb --in protein_pool.faa -d protein_pool --threads ${task.cpus} --quiet
    """
}

/**
 * Recover missing genes for a batch of genomes, as GFF3 for antiSMASH.
 *
 * Both prodigal and homology, because neither suffices alone. On LMG 5342 prodigal
 * found 8 of the 9 missing genes but truncated aepZ to 162 aa where the gene is 238 --
 * and aepZ is the gene the whole exercise is about. Homology recovers it at 100%
 * identity and names it.
 *
 * The GFF3s are published because they are the evidence for every downstream gene call
 * that was not in the source annotation.
 *
 * One task handles a batch of genomes with ONE blastx over all of their gaps -- see the
 * note in the script block, and run_batch() in recover_orfs.py for why that needs three
 * passes rather than two.
 */
process RECOVER_ORFS {
    tag "${genomes instanceof List ? genomes.size() + ' genomes' : genomes.baseName}"
    label 'process_medium'

    // Scales with the DATABASE, for the same portability reason as above.
    //
    // It is diamond that holds the memory. Measured against the whole-Bacteria
    // pool: a 7.19 GB database peaks at 7.26 GB on a full batch and 3.46 GB on a
    // single genome's gaps, so peak tracks database size closely at that point,
    // and 6 GB over it is comfortable headroom.
    //
    // CAVEAT, because one point is not a curve: that relationship is measured at
    // ONE database size. diamond chunks its database, so peak may well be BOUNDED
    // rather than linear above 7 GB -- in which case this over-requests on a large
    // pool. Over-requesting wastes an allocation; under-requesting is an OOM kill
    // of a long batch, so the error is taken in this direction, and the escalation
    // below absorbs the case where it is wrong the other way. Re-measure at a
    // larger pool before trusting the high end.
    memory { 1.GB * (6 + Math.ceil(protein_pool_db.size() / 1e9)) * task.attempt }
    maxRetries 1
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140] && task.attempt <= 1 ? 'retry' : 'terminate' }
    publishDir "${params.outdir}/recovered_orfs/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode, pattern: '*.gff3'

    input:
    val taxon
    path genomes
    path protein_pool
    path protein_pool_db
    val scripts_version

    output:
    path "*.gff3",            emit: gff
    path "recovery_stats.tsv", emit: stats, optional: true

    script:
    """
    # ONE call for the whole batch, not one per genome. diamond spends ~120 CPU-s
    # loading the 7.19 GB pool database before it searches anything, so the loop
    # this replaces paid that cost once per genome: measured on 20 genomes of the
    # whole-Bacteria pilot, 20 calls cost 2,536 CPU-s against 265 for one, and the
    # batched run reproduced all 20 genomes' hits exactly. That is ~90% of this
    # stage, which was 76.4% of the pilot's 428.5 CPU-h.
    #
    # Per-genome failures stay per-genome inside the script -- a genome that
    # cannot be parsed gets an empty GFF3 and the batch continues, exactly as the
    # shell loop did. A non-zero exit here therefore means the batch itself broke
    # (a missing binary, an unreadable database), which is not something a
    # fallback should paper over.
    python ${projectDir}/scripts/genome/recover_orfs.py \\
        --genomes ${genomes} \\
        --outdir . \\
        --pool ${protein_pool} \\
        --pool_db ${protein_pool_db} \\
        --prodigal \$(which prodigal) \\
        --diamond \$(which diamond) \\
        --threads ${task.cpus} \\
        --workdir _work \\
        > recovery_stats.tsv 2>&1
    rm -rf _work
    """
}
