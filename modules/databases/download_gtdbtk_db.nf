/**
 * Download GTDB-Tk reference database (~140 GB)
 * Required for phylogenetic classification with GTDB-Tk
 * Uses parallel downloads via aria2 for faster completion
 */
process DOWNLOAD_GTDBTK_DB {
    tag "GTDB-Tk database"
    label 'process_medium'
    storeDir "${params.outdir}/databases/gtdbtk"

    output:
    // Version-specific on purpose. storeDir skips the download when this path
    // already exists, so naming the release here is what makes changing
    // params.gtdb_release actually fetch new data instead of silently reusing
    // whatever was downloaded first. It also matches the existing on-disk
    // layout (databases/gtdbtk/release226), so pinning costs no re-download.
    path "release${params.gtdb_release}", emit: db_dir

    script:
    """
    echo "=============================================="
    echo "GTDB-Tk Database Download"
    echo "=============================================="
    echo "WARNING: this archive is ~132 GB and extracts to ~140 GB"
    echo ""

    # Mirror CHANGED 2026-11-02. This used to fetch a 14-part split package from
    # data.gtdb.aau.ecogenomic.org, whose TLS certificate EXPIRED on 2026-10-27:
    #
    #     SSL/TLS handshake failure: certificate has expired
    #
    # All 14 parts failed at 0 B/s. That mirror also publishes the only copy of the
    # split layout -- the same path 404s on the other two -- so moving mirrors means
    # moving to the single full package, which is why the concatenation step is gone.
    #
    # Certificate expiry dates checked the same day:
    #     data.gtdb.aau.ecogenomic.org   expired 2026-10-27   <- was hardcoded here
    #     data.gtdb.ecogenomic.org       2026-11-23           <- used now
    #     data.ace.uq.edu.au             2027-01-19           (404s on this path)
    #
    # A mirror's certificate will expire again. The symptom is the handshake error
    # above on every part at once; check the dates with
    #     echo | openssl s_client -connect <host>:443 -servername <host> \\
    #       | openssl x509 -noout -dates
    ARCHIVE_URL="https://data.gtdb.ecogenomic.org/releases/release${params.gtdb_release}/${params.gtdb_release}.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r${params.gtdb_release}_data.tar.gz"

    echo "Downloading from: \$ARCHIVE_URL"
    echo "Started at: \$(date)"
    echo ""

    # Same stale-CA guard as download_pfam.nf: a conda environment's CA bundle is
    # frozen on the day it was solved. aria2c reads SSL_CERT_FILE.
    for _ca in /etc/ssl/certs/ca-certificates.crt \\
               /etc/pki/tls/certs/ca-bundle.crt \\
               "\${CONDA_PREFIX:-}/ssl/cacert.pem"; do
        if [ -r "\$_ca" ]; then export SSL_CERT_FILE="\$_ca"; break; fi
    done
    echo "CA bundle: \${SSL_CERT_FILE:-<none found; using the default>}"

    # One file, so the parallelism is connections-per-file rather than files. The
    # mirror advertises accept-ranges: bytes, which is what -x/-s need. --continue
    # lets an interrupted 132 GB transfer resume instead of restarting.
    aria2c -x 16 -s 16 -k 20M \\
        --continue=true \\
        --max-tries=0 \\
        --retry-wait=10 \\
        --file-allocation=none \\
        --console-log-level=notice \\
        --summary-interval=60 \\
        -o gtdbtk_data.tar.gz \\
        "\$ARCHIVE_URL"

    echo ""
    echo "Download finished at: \$(date)"
    echo "Archive size: \$(du -h gtdbtk_data.tar.gz)"

    echo ""
    echo "Extracting database (this may take 30-60 minutes)..."
    tar -xzf gtdbtk_data.tar.gz

    # Clean up archive to save space
    rm gtdbtk_data.tar.gz

    # Find the release directory (e.g., release220, release226)
    RELEASE_DIR=\$(ls -d release* 2>/dev/null | head -1)

    if [ -z "\$RELEASE_DIR" ]; then
        echo "ERROR: Could not find extracted GTDB-Tk release directory"
        exit 1
    fi

    echo ""
    echo "GTDB-Tk database downloaded successfully to \$RELEASE_DIR"
    echo "Database contents:"
    ls -lh \$RELEASE_DIR/

    echo ""
    echo "=============================================="
    echo "GTDB-Tk database ready for use"
    echo "=============================================="
    """
}
