/**
 * Get the installed antiSMASH version.
 * Used for tracking which version produced results.
 */
process GET_ANTISMASH_VERSION {
    label 'process_local'
    executor 'local'

    output:
    stdout emit: version

    script:
    """
    antismash --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 | tr -d '\\n'
    """
}

/**
 * antiSMASH over a batch of genomes.
 *
 * Batched because one task per genome makes a large run bound by the scheduler's
 * submission rate rather than by compute: measured at 41.4 s per genome, a million
 * genomes is 34.7 days of submission against 2.4 days of compute. See
 * antismashBatchSize().
 *
 * Each genome is run inside a shell loop that *continues* past a failure rather
 * than letting the task exit non-zero. That matters more once batched: previously
 * a failed genome cost one genome, now an aborting task would cost the whole
 * batch. The `tolerant` label and the retry on the conda-startup race still apply,
 * but they are the outer net, not the primary mechanism.
 *
 * Results are written under as_out/ so the output glob cannot match the staged
 * database directory, and `saveAs` strips that prefix so the published layout is
 * unchanged from the unbatched version.
 */
process ANTISMASH {
    tag "${genomes instanceof List ? genomes.size() + ' genomes' : genomes.baseName}"
    label 'process_medium'
    label 'tolerant'
    cache 'lenient'
    publishDir "${params.outdir}/antismash_results/${Utils.sanitizeTaxon(params.taxon)}",
        mode: params.publish_mode,
        saveAs: { fn -> fn.startsWith('as_out/') ? fn.substring(7) : fn }

    input:
    val taxon
    path genomes
    // Recovered gene calls, one <genome>.gff3 per genome, staged flat alongside the
    // genomes and paired by basename inside the loop. May be a placeholder when
    // recovery is disabled.
    path recovered_gff
    path antismash_db
    val antismash_version
    val antismash_params_hash

    // Digest of the Python this process runs. A val input, not an interpolation:
    // Nextflow hashes the unevaluated script source plus the input values. See CLAUDE.md.
    val scripts_version

    output:
    path "as_out/*", emit: result_dir, optional: true

    script:
    // Phosphonate-only detection — hardcoded
    def hmmdetection_flag = '--hmmdetection-limit-to-rule-names phosphonate'

    // Build minimal mode flag
    def minimal_flag = params.antismash_minimal ? '--minimal' : ''

    // When in minimal mode, enable HTML output for visualization compatibility
    def html_output_flag = params.antismash_minimal ? '--enable-html' : ''

    // Build antiSMASH analysis flags based on config (ClusterBlast options ignored in minimal mode)
    def cb_general_flag = params.antismash_minimal ? '' : (params.antismash_cb_general ? '--cb-general' : '')
    def cc_mibig_flag = params.antismash_minimal ? '' : (params.antismash_cc_mibig ? '--cc-mibig' : '')
    def cb_knownclusters_flag = params.antismash_minimal ? '' : '--cb-knownclusters'   // hardcoded: always compare vs MIBiG
    def smcog_trees_flag = params.antismash_minimal ? '' : (params.antismash_smcog_trees ? '--smcog-trees' : '')

    // Whole-genome summary GenBank, {genome}.gbk (~11 MB each, the single largest
    // file antiSMASH writes). Off by default: no pipeline step reads it, and
    // BiG-SCAPE only ingests .gbk filenames containing "cluster" or "region", so
    // the summary is filtered out of clustering regardless. Worth enabling on small
    // sets you want to open in a genome browser. Like --no-zip-output it stays out
    // of Utils.antismashParamsHash, so toggling it will NOT regenerate genomes
    // already present in a --reuse_antismash_from directory.
    def summary_gbk_flag = params.antismash_summary_gbk ? '--summary-gbk' : '--no-summary-gbk'

    // --no-zip-output (hardcoded below): antiSMASH otherwise writes {genome}.zip,
    // an archive of the very directory it sits in - ~6 MB per genome of pure
    // redundancy that nothing downstream reads. Deliberately absent from
    // Utils.antismashParamsHash: it changes packaging, not results, so results
    // generated before this flag stay reusable via --reuse_antismash_from.

    // Domain analysis flags - always enabled when not in minimal mode
    def clusterhmmer_flag = params.antismash_minimal ? '' : '--clusterhmmer'
    def tigrfam_flag = params.antismash_minimal ? '' : '--tigrfam'

    // antiSMASH sizes a region as the rule core plus a fixed neighbourhood, and the
    // strict phosphonate rule's 5 kb stops 4.2 kb short of the HiVir cluster. There is
    // no CLI option for it, so the installed rule file is edited in place first; the
    // value is part of antismashParamsHash, so results are never reused across it.
    def neighbourhood_patch = params.antismash_phosphonate_neighbourhood
        ? "python ${projectDir}/scripts/genome/patch_antismash_neighbourhood.py " +
          "--rule phosphonate --neighbourhood ${params.antismash_phosphonate_neighbourhood}"
        : "echo 'phosphonate neighbourhood: antiSMASH default'"

    """
    # Region size: see scripts/genome/patch_antismash_neighbourhood.py
    ${neighbourhood_patch}

    # Set antiSMASH database location
    export ANTISMASH_DB_PATH=\$(readlink -f ${antismash_db})

    # Create symlink for compatibility
    mkdir -p ~/.local/share
    ln -sf \$ANTISMASH_DB_PATH ~/.local/share/antismash

    mkdir -p as_out
    OK=0
    SKIPPED=0

    # One genome per iteration. Every failure path is a `continue`, never an exit:
    # an aborting task would forfeit the whole batch, not one genome.
    for GENOME in ${genomes}; do
        BASE=\$(basename "\$GENOME" .gbff)
        echo "=== \$BASE"

        if [ ! -s "\$GENOME" ]; then
            echo "  SKIP: input file is empty or missing"
            SKIPPED=\$((SKIPPED + 1))
            continue
        fi

        LOCUS_COUNT=\$(grep -c "^LOCUS" "\$GENOME" || echo "0")
        LOCUS_COUNT=\$(echo "\$LOCUS_COUNT" | tr -d '[:space:]')
        if [ "\$LOCUS_COUNT" -eq 0 ]; then
            echo "  SKIP: no LOCUS records"
            SKIPPED=\$((SKIPPED + 1))
            continue
        fi

        # antiSMASH needs genes; use its own caller only when the file has none.
        CDS_COUNT=\$(grep -c "^     CDS" "\$GENOME" 2>/dev/null || echo "0")
        CDS_COUNT=\$(echo "\$CDS_COUNT" | tr -d '[:space:]')
        if [ "\$CDS_COUNT" -gt 0 ]; then
            GENEFINDING="none"
        else
            GENEFINDING="prodigal"
        fi

        # Genes the submitter never called, recovered by RECOVER_ORFS. antiSMASH runs
        # its own gene finder ONLY on records with zero CDS features, so a partially
        # annotated genome would otherwise keep only the genes it shipped with. GFF3
        # features are merged before that check, which is why this works at all.
        # A header-only GFF3 (no recovered genes) is skipped: antiSMASH rejects a
        # GFF3 containing no CDS.
        GFF_FLAG=""
        if [ -s "\${BASE}.gff3" ] && grep -qv '^#' "\${BASE}.gff3"; then
            GFF_FLAG="--genefinding-gff3 \${BASE}.gff3"
            RECOVERED=\$(grep -vc '^#' "\${BASE}.gff3" || echo 0)
            echo "  +\$RECOVERED recovered genes"
        fi
        echo "  \$LOCUS_COUNT LOCUS, \$CDS_COUNT CDS, genefinding=\$GENEFINDING"

        if antismash \\
            --taxon bacteria \\
            --output-dir "as_out/\$BASE" \\
            --genefinding-tool \$GENEFINDING \\
            \$GFF_FLAG \\
            --databases \$ANTISMASH_DB_PATH \\
            --cpus ${task.cpus} \\
            --allow-long-headers \\
            --hmmdetection-strictness strict \\
            --no-zip-output \\
            ${summary_gbk_flag} \\
            ${minimal_flag} \\
            ${html_output_flag} \\
            ${hmmdetection_flag} \\
            ${cb_general_flag} \\
            ${cc_mibig_flag} \\
            ${cb_knownclusters_flag} \\
            ${smcog_trees_flag} \\
            ${clusterhmmer_flag} \\
            ${tigrfam_flag} \\
            "\$GENOME"; then
            # Consumed by CHECK_ANTISMASH_REUSE to decide whether a later run
            # can reuse this directory.
            cat > "as_out/\$BASE/.antismash_meta" << EOF
version=${antismash_version}
params_hash=${antismash_params_hash}
EOF
            OK=\$((OK + 1))
        else
            echo "  SKIP: antiSMASH failed"
            rm -rf "as_out/\$BASE"
            SKIPPED=\$((SKIPPED + 1))
        fi
    done

    echo "batch complete: \$OK succeeded, \$SKIPPED skipped"

    # Only a batch where nothing worked is worth failing — that signals a broken
    # environment (the conda startup race, a missing database) rather than bad
    # input, and the retry in conf/labels.config is the right response.
    if [ "\$OK" -eq 0 ]; then
        echo "ERROR: no genome in this batch produced output"
        exit 1
    fi
    exit 0
    """
}
