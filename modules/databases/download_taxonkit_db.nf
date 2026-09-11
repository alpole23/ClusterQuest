/**
 * Download NCBI taxonomy database (~100 MB)
 * Required for TaxonKit taxonomy processing
 */
process DOWNLOAD_TAXONKIT_DB {
    tag "NCBI taxonomy"
    label 'process_low'
    storeDir "${params.outdir}/databases/taxonkit"

    output:
    // Version in the path, for the same reason as Pfam: the existing
    // databases/taxonkit/taxdump is a 2026-01-23 copy, and without the date here
    // storeDir would reuse it while the run reported the pinned date.
    path "taxdump_${params.taxdump_date}", emit: taxdump_dir

    script:
    """
    echo "=============================================="
    echo "Downloading NCBI Taxonomy Database"
    echo "=============================================="
    echo "Size: ~100 MB"
    echo "Started at: \$(date)"
    echo ""

    # Create directory for taxonomy dump
    mkdir -p taxdump_${params.taxdump_date}

    # A dated monthly snapshot, not the live taxdump.tar.gz. The live file is
    # rewritten daily, so with storeDir it would pin itself to whatever day the
    # pipeline first ran — this pipeline sat on a 2026-01-23 copy for seven
    # months that way. The archive is .zip; the live file is .tar.gz.
    echo "Downloading taxdmp_${params.taxdump_date}.zip from NCBI FTP..."
    wget -q --show-progress -O taxdmp.zip \\
        "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${params.taxdump_date}.zip"

    echo ""
    echo "Extracting taxonomy files..."
    unzip -q -o taxdmp.zip -d taxdump_${params.taxdump_date}
    rm taxdmp.zip

    echo "taxdump=${params.taxdump_date}" > taxdump_${params.taxdump_date}/.db_version

    # Verify required files
    echo ""
    echo "Verifying download..."
    REQUIRED_FILES="nodes.dmp names.dmp delnodes.dmp merged.dmp"
    for f in \$REQUIRED_FILES; do
        if [ ! -f "taxdump_${params.taxdump_date}/\$f" ]; then
            echo "ERROR: Required file taxdump_${params.taxdump_date}/\$f not found"
            exit 1
        fi
    done

    echo ""
    echo "Database contents:"
    ls -lh taxdump_${params.taxdump_date}/

    echo ""
    echo "=============================================="
    echo "NCBI taxonomy database ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
