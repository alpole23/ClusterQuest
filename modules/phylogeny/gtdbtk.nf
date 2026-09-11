/**
 * GTDB-Tk classification over one shard of genomes.
 *
 * Sharded so that wall time is not set by a single monolithic task; see
 * gtdbtkShardSize(). Runs below the shard size get exactly one shard, which is
 * the previous behaviour.
 *
 * **No tree is emitted.** classify_wf still builds one internally — pplacer
 * placement is how it classifies — but a per-shard tree spans a disjoint genome
 * set, and N such trees cannot be concatenated into one phylogeny. The taxonomy
 * assignments merge cleanly by row and are what every downstream consumer
 * actually uses.
 */
process GTDBTK_CLASSIFY {
    tag "${taxon} shard ${shard_id}"
    label 'process_high_memory'
    cache 'lenient'  // GTDB-Tk is memory-intensive - use lenient caching

    input:
    val taxon
    tuple val(shard_id), path(fasta_files)
    path gtdbtk_db

    output:
    // Summary files are at top-level in GTDB-Tk v2.x output. Renamed per shard so
    // MERGE_GTDBTK can stage them all in one directory without collision.
    path "gtdbtk.bac120.summary.${shard_id}.tsv", emit: bacterial_summary, optional: true
    path "gtdbtk.ar53.summary.${shard_id}.tsv",   emit: archaeal_summary,  optional: true

    script:
    def cpus = params.gtdbtk_cpus ?: task.cpus ?: 8  // Use dedicated param, fallback to task.cpus
    def pplacer_cpus = params.gtdbtk_pplacer_cpus ?: 1  // pplacer is memory-heavy, use 1 by default
    def skip_ani = "--skip_ani_screen"
    def min_perc_aa = params.gtdbtk_min_perc_aa ?: 10
    """
    echo "=============================================="
    echo "GTDB-Tk Classification"
    echo "=============================================="
    echo "Taxon: ${taxon}"
    echo "CPUs: ${cpus}"
    echo "Pplacer CPUs: ${pplacer_cpus}"
    echo ""

    # Set GTDB-Tk database path
    RELEASE_DIR=\$(ls -d ${gtdbtk_db}/release* 2>/dev/null | head -1)
    if [ -z "\$RELEASE_DIR" ]; then
        export GTDBTK_DATA_PATH="${gtdbtk_db}"
    else
        export GTDBTK_DATA_PATH="\$RELEASE_DIR"
    fi
    echo "Using GTDB-Tk database: \$GTDBTK_DATA_PATH"

    # Create input genome directory
    # Add unique prefix to avoid conflicts with GTDB reference genome names
    # This prevents "duplicate taxon labels" errors in tree generation
    mkdir -p genomes
    for f in ${fasta_files}; do
        BASENAME=\$(basename "\$f")
        cp "\$f" "genomes/usr_\${BASENAME}"
    done

    GENOME_COUNT=\$(ls genomes/*.fna 2>/dev/null | wc -l)
    echo "Processing \$GENOME_COUNT genomes"
    echo ""

    if [ "\$GENOME_COUNT" -eq 0 ]; then
        echo "ERROR: No FASTA files found"
        exit 1
    fi

    # Run GTDB-Tk classify workflow
    echo "Running GTDB-Tk classify_wf..."
    echo "This may take several hours for large datasets."
    echo ""

    # Create scratch directory to reduce peak memory usage
    mkdir -p scratch_tmp

    # Run GTDB-Tk - capture exit code but don't fail immediately
    # Tree generation can fail with duplicate taxon labels (known GTDB-Tk issue)
    # when user genomes match reference strain names exactly
    set +e
    gtdbtk classify_wf \\
        --genome_dir genomes \\
        --out_dir gtdbtk_output \\
        --extension fna \\
        --cpus ${cpus} \\
        --pplacer_cpus ${pplacer_cpus} \\
        --min_perc_aa ${min_perc_aa} \\
        --scratch_dir scratch_tmp \\
        ${skip_ani} \\
        --prefix gtdbtk
    GTDBTK_EXIT=\$?
    set -e

    echo ""
    echo "=============================================="
    echo "GTDB-Tk classification complete (exit code: \$GTDBTK_EXIT)"
    echo "=============================================="

    # Check if we have summary files (main output) even if tree generation failed
    if [ \$GTDBTK_EXIT -ne 0 ]; then
        if [ -f gtdbtk_output/gtdbtk.bac120.summary.tsv ] || [ -f gtdbtk_output/gtdbtk.ar53.summary.tsv ]; then
            echo "WARNING: GTDB-Tk had errors but summary files were generated."
            echo "This often occurs due to duplicate taxon labels in tree generation."
            echo "Classification results are still valid."
        else
            echo "ERROR: GTDB-Tk failed and no summary files were generated."
            exit 1
        fi
    fi

    # Check outputs
    if [ -f gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
        BAC_COUNT=\$(tail -n +2 gtdbtk_output/gtdbtk.bac120.summary.tsv | wc -l)
        echo "Bacterial genomes classified: \$BAC_COUNT"
    fi
    if [ -f gtdbtk_output/gtdbtk.ar53.summary.tsv ]; then
        AR_COUNT=\$(tail -n +2 gtdbtk_output/gtdbtk.ar53.summary.tsv | wc -l)
        echo "Archaeal genomes classified: \$AR_COUNT"
    fi

    # Shard-tagged copies for the merge. The tree classify_wf produced is left in
    # gtdbtk_output and not emitted: see the process docstring.
    for DOMAIN in bac120 ar53; do
        if [ -f "gtdbtk_output/gtdbtk.\${DOMAIN}.summary.tsv" ]; then
            cp "gtdbtk_output/gtdbtk.\${DOMAIN}.summary.tsv" \\
               "gtdbtk.\${DOMAIN}.summary.${shard_id}.tsv"
        fi
    done
    """
}
