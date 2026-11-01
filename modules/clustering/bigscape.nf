process BIGSCAPE {
    tag "$taxon"
    label 'process_high'
    publishDir "${params.outdir}/bigscape_results", mode: params.publish_mode
    cache 'lenient'  // Use lenient caching for directory inputs

    input:
    val taxon
    val opts_hash
    path antismash_results
    path pfam_db
    
    output:
    path "${Utils.sanitizeTaxon(params.taxon)}/", emit: bigscape_dir
    path "${Utils.sanitizeTaxon(params.taxon)}/${Utils.sanitizeTaxon(params.taxon)}.db", emit: bigscape_db

    script:
    def safe_taxon = Utils.sanitizeTaxon(params.taxon)
    def cutoffs = params.bigscape_cutoffs ?: "0.30"
    def alignment_mode = params.bigscape_alignment_mode ?: "auto"
    def mibig_version = params.bigscape_mibig_version ? "--mibig-version ${params.bigscape_mibig_version}" : ""
    // opts_hash is a declared input purely so a changed BiG-SCAPE option
    // invalidates this task; interpolated values never would.
    def min_len = (params.bigscape_min_bgc_length ?: 0) as Integer
    def classify = params.bigscape_classify ? "--classify ${params.bigscape_classify}" : ""
    def include_singletons = params.bigscape_include_singletons ? "--include-singletons" : ""
    def mix = params.bigscape_mix ? "--mix" : ""
    """
    # BiG-SCAPE 2 requires path to Pfam-A.hmm file
    export PFAM_PATH=\$(readlink -f ${pfam_db})/Pfam-A.hmm
    # BiG-SCAPE loads GBKs in the iteration order of a set keyed on a hash STRING, so
    # without this the load order — and with it each pair's A/B orientation, which its
    # extended comparisons are not symmetric under — changes run to run. Measured:
    # 8-14% of distances differed between identical runs, and family centres moved.
    export PYTHONHASHSEED=0
    
    # BiG-SCAPE has NO --min-bgc-length option. The minimum lives in its config.yml
    # as MIN_BGC_LENGTH, applied via --config-file-path, and parse_config() reads
    # every key with config["KEY"] rather than .get(), so a partial file raises
    # KeyError -- the bundled config must be copied and patched, not replaced.
    #
    # Writing this as a CLI flag cost a validation run: `--min-bgc-length 2000` was
    # accepted by nextflow, reached bigscape, and died with "No such option" after
    # 134 tasks had completed. The flag was never checked against `bigscape cluster
    # --help`. Check the tool, not the database column it writes.
    CFG=""
    if [ "${min_len}" -gt 0 ]; then
        BUNDLED=\$(python -c "import big_scape, pathlib; print(pathlib.Path(big_scape.__file__).parent / 'config.yml')")
        sed "s/^MIN_BGC_LENGTH:.*/MIN_BGC_LENGTH: ${min_len}/" "\$BUNDLED" > bigscape_config.yml
        grep -q "^MIN_BGC_LENGTH: ${min_len}\$" bigscape_config.yml || {
            echo "ERROR: could not set MIN_BGC_LENGTH in the BiG-SCAPE config" >&2; exit 1; }
        CFG="--config-file-path bigscape_config.yml"
        echo "MIN_BGC_LENGTH set to ${min_len} bp via bigscape_config.yml"
    fi

    echo "Running BiG-SCAPE on ${taxon}"
    echo "Using Pfam database: \$PFAM_PATH"
    echo "Using ${task.cpus} cores"
    
    # Create input directory with all antiSMASH results
    mkdir -p antismash_input
    
    # Copy/link all antiSMASH result directories
    for dir in ${antismash_results}; do
        if [ -d "\$dir" ]; then
            echo "Linking \$dir"
            ln -s \$(readlink -f "\$dir") antismash_input/
        fi
    done
    
    # Count how many we have
    INPUT_COUNT=\$(ls -d antismash_input/*/ 2>/dev/null | wc -l)
    echo "Processing \$INPUT_COUNT antiSMASH results"
    
    if [ "\$INPUT_COUNT" -eq 0 ]; then
        echo "ERROR: No antiSMASH results found"
        echo "Listing what we received:"
        ls -la
        exit 1
    fi
    
    # Run BiG-SCAPE 2
    bigscape cluster \
        -i antismash_input \
        -o ${safe_taxon} \
        --pfam-path \$PFAM_PATH \
        --alignment-mode ${alignment_mode} \
        --gcf-cutoffs ${cutoffs} \
        \$CFG \
        ${include_singletons} \
        ${mix} \
        ${mibig_version} \
        --cores ${task.cpus} \
        ${classify}

    # Update index.html to use taxon-specific database name
    # BiG-SCAPE generates index.html with hardcoded 'data_sqlite.db' reference
    # Replace it with the actual taxon-specific database filename
    cd ${safe_taxon}
    sed -i "s|data_sqlite.db|${safe_taxon}.db|g" index.html
    cd ..

    echo "BiG-SCAPE analysis complete"
    """
}