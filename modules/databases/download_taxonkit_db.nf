/**
 * Download NCBI taxonomy database (~100 MB)
 * Required for TaxonKit taxonomy processing
 */
process DOWNLOAD_TAXONKIT_DB {
    tag "NCBI taxonomy"
    label 'process_low'
    storeDir "${params.outdir}/databases/taxonkit"

    output:
    path "taxdump", emit: taxdump_dir

    script:
    """
    echo "=============================================="
    echo "Downloading NCBI Taxonomy Database"
    echo "=============================================="
    echo "Size: ~100 MB"
    echo "Started at: \$(date)"
    echo ""

    # Create directory for taxonomy dump
    mkdir -p taxdump

    # A dated monthly snapshot, not the live taxdump.tar.gz. The live file is
    # rewritten daily, so with storeDir it would pin itself to whatever day the
    # pipeline first ran — this pipeline sat on a 2026-01-23 copy for seven
    # months that way. The archive is .zip; the live file is .tar.gz.
    echo "Downloading taxdmp_${params.taxdump_date}.zip from NCBI FTP..."
    wget -q --show-progress -O taxdmp.zip \\
        "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${params.taxdump_date}.zip"

    echo ""
    echo "Extracting taxonomy files..."
    unzip -q -o taxdmp.zip -d taxdump
    rm taxdmp.zip

    echo "taxdump=${params.taxdump_date}" > taxdump/.db_version

    # Verify required files
    echo ""
    echo "Verifying download..."
    REQUIRED_FILES="nodes.dmp names.dmp delnodes.dmp merged.dmp"
    for f in \$REQUIRED_FILES; do
        if [ ! -f "taxdump/\$f" ]; then
            echo "ERROR: Required file taxdump/\$f not found"
            exit 1
        fi
    done

    echo ""
    echo "Database contents:"
    ls -lh taxdump/

    echo ""
    echo "=============================================="
    echo "NCBI taxonomy database ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
