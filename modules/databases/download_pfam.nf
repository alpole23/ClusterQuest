/**
 * Download Pfam-A HMM database (~1.5 GB)
 * Required for BiG-SCAPE clustering
 */
process DOWNLOAD_PFAM {
    tag "Pfam database"
    label 'process_low'
    storeDir "${params.outdir}/databases"

    output:
    // Version in the path on purpose: storeDir skips the download whenever this
    // path exists, so without the release here a changed params.pfam_release
    // would be recorded in software_versions.json while the old data was
    // silently reused — a false provenance claim, worse than recording nothing.
    path "pfam${params.pfam_release}/", emit: pfam_db

    script:
    """
    echo "=============================================="
    echo "Downloading Pfam Database"
    echo "=============================================="
    echo "Size: ~1.5 GB"
    echo "Started at: \$(date)"
    echo ""

    mkdir -p pfam${params.pfam_release}
    cd pfam${params.pfam_release}

    # A conda environment's bundled CA store is fixed on the day the environment was
    # solved and is never refreshed. This one carried a 2026-07-22 bundle that does not
    # include the issuer ftp.ebi.ac.uk now presents, so wget died with
    #     ERROR: cannot verify ftp.ebi.ac.uk's certificate ... Self-signed certificate
    # and exit status 5 -- which `-q` SUPPRESSED, so the task failed with an empty
    # .command.err and a bare "exit status (5)". Prefer the system store, which the OS
    # keeps current; fall back to the environment's own. Note `wget --spider` does NOT
    # reproduce the failure: it verifies nothing, so a HEAD probe passes while the real
    # transfer fails.
    for _ca in /etc/ssl/certs/ca-certificates.crt \\
               /etc/pki/tls/certs/ca-bundle.crt \\
               "\${CONDA_PREFIX:-}/ssl/cacert.pem"; do
        if [ -r "\$_ca" ]; then export SSL_CERT_FILE="\$_ca"; break; fi
    done
    echo "CA bundle: \${SSL_CERT_FILE:-<none found; using the default>}"

    # Download Pfam-A HMM profiles
    if [ ! -f "Pfam-A.hmm.gz" ] && [ ! -f "Pfam-A.hmm" ]; then
        echo "Downloading Pfam-A.hmm..."
        # NOT -q: it hid the certificate error above. dot:giga keeps the log short
        # while leaving errors visible.
        wget --progress=dot:giga https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam${params.pfam_release}/Pfam-A.hmm.gz
        echo "Extracting..."
        gunzip Pfam-A.hmm.gz
    else
        echo "Pfam-A.hmm already exists"
    fi

    # Press the HMM database for faster searching
    if [ ! -f "Pfam-A.hmm.h3p" ]; then
        echo "Pressing Pfam database with hmmpress..."
        hmmpress Pfam-A.hmm
    else
        echo "Pfam database already pressed"
    fi

    # Nothing inside Pfam-A.hmm identifies its release, so stamp it.
    echo "pfam=${params.pfam_release}" > .db_version

    # Verify
    echo ""
    echo "Verifying download..."
    if [ ! -f "Pfam-A.hmm" ] || [ ! -f "Pfam-A.hmm.h3p" ]; then
        echo "ERROR: Pfam database files missing"
        exit 1
    fi

    echo ""
    echo "Database contents:"
    ls -lh

    echo ""
    echo "=============================================="
    echo "Pfam database ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
