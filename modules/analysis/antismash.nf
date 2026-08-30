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

process ANTISMASH {
    tag "$genome.baseName"
    label 'process_medium'
    label 'tolerant'
    cache 'lenient'
    publishDir "${params.outdir}/antismash_results/${Utils.sanitizeTaxon(params.taxon)}", mode: params.publish_mode

    input:
    val taxon
    path genome
    path antismash_db
    val antismash_version
    val antismash_params_hash

    output:
    path "${genome.baseName}/", emit: result_dir, optional: true

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

    """
    # Set antiSMASH database location
    export ANTISMASH_DB_PATH=\$(readlink -f ${antismash_db})

    # Create symlink for compatibility
    mkdir -p ~/.local/share
    ln -sf \$ANTISMASH_DB_PATH ~/.local/share/antismash

    echo "Processing: ${genome.baseName}"

    # Check if input file exists and has content
    if [ ! -s ${genome} ]; then
        echo "ERROR: Input file is empty or missing"
        exit 1
    fi

    # Count LOCUS records
    LOCUS_COUNT=\$(grep -c "^LOCUS" ${genome} || echo "0")
    LOCUS_COUNT=\$(echo "\$LOCUS_COUNT" | tr -d '[:space:]')
    echo "Found \$LOCUS_COUNT LOCUS record(s)"

    if [ "\$LOCUS_COUNT" -eq 0 ]; then
        echo "ERROR: No LOCUS records found"
        exit 1
    fi

    # Check if GenBank file has CDS features
    CDS_COUNT=\$(grep -c "^     CDS" ${genome} 2>/dev/null || echo "0")
    CDS_COUNT=\$(echo "\$CDS_COUNT" | tr -d '[:space:]')
    echo "Found \$CDS_COUNT CDS features"

    # Determine gene finding tool
    if [ "\$CDS_COUNT" -gt 0 ]; then
        GENEFINDING="none"
        echo "Using existing gene annotations"
    else
        GENEFINDING="prodigal"
        echo "No CDS features found, using prodigal for gene finding"
    fi

    # Run antiSMASH
    antismash \\
        --taxon bacteria \\
        --output-dir ${genome.baseName} \\
        --genefinding-tool \$GENEFINDING \\
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
        ${genome}

    ANTISMASH_EXIT=\$?

    if [ \$ANTISMASH_EXIT -eq 0 ]; then
        echo "SUCCESS: antiSMASH completed for ${genome.baseName}"

        # Record version and params hash for result reuse tracking
        cat > "${genome.baseName}/.antismash_meta" << EOF
version=${antismash_version}
params_hash=${antismash_params_hash}
EOF
        exit 0
    else
        echo "WARNING: antiSMASH failed with exit code \$ANTISMASH_EXIT"
        exit 0
    fi
    """
}
