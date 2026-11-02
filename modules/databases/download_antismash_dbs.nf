/**
 * Download antiSMASH databases (~50 GB)
 * Required for BGC detection with antiSMASH
 */
process DOWNLOAD_ANTISMASH_DBS {
    tag "antiSMASH databases"
    label 'process_low'
    storeDir "${params.outdir}/databases"

    output:
    path "antismash/", emit: db_dir

    script:
    """
    echo "=============================================="
    echo "Downloading antiSMASH Databases"
    echo "=============================================="
    echo "Size: ~50 GB"
    echo "Started at: \$(date)"
    echo ""

    # Set the database directory
    export ANTISMASH_DB_PATH=\$(pwd)/antismash
    mkdir -p \$ANTISMASH_DB_PATH

    # Same stale-CA problem as download_pfam.nf, one layer down: this downloader is
    # Python urllib, which fails with
    #     URLError: [SSL: CERTIFICATE_VERIFY_FAILED] unable to get local issuer certificate
    # and then reports "File not found on server. Please check your internet connection",
    # which points at the network rather than at the local CA store. Measured 2026-11-02:
    # it fails on the FIRST item and downloads nothing.
    for _ca in /etc/ssl/certs/ca-certificates.crt \\
               /etc/pki/tls/certs/ca-bundle.crt \\
               "\${CONDA_PREFIX:-}/ssl/cacert.pem"; do
        if [ -r "\$_ca" ]; then export SSL_CERT_FILE="\$_ca"; break; fi
    done
    echo "CA bundle: \${SSL_CERT_FILE:-<none found; using the default>}"

    # antiSMASH locates its own tools by probing PATH, so hmmpress must be reachable or
    # the final pre-build step dies with "cannot find executable for hmmpress" after
    # every archive has downloaded. Under Nextflow the conda env is activated and its
    # bin is already on PATH; this makes that explicit rather than incidental.
    export PATH="\${CONDA_PREFIX:-}/bin:\$PATH"

    # Download databases
    echo "Downloading databases..."
    download-antismash-databases --database-dir \$ANTISMASH_DB_PATH

    # Verify download
    echo ""
    echo "Verifying download..."
    if [ ! -d "\$ANTISMASH_DB_PATH" ] || [ -z "\$(ls -A \$ANTISMASH_DB_PATH)" ]; then
        echo "ERROR: Database directory is empty"
        exit 1
    fi

    echo ""
    echo "Database contents:"
    ls -lh \$ANTISMASH_DB_PATH/

    echo ""
    echo "=============================================="
    echo "antiSMASH databases ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
