#!/bin/bash
# =============================================================================
# ClusterQuest SLURM Submission Script
# =============================================================================
#
# Usage:
#   sbatch submit_slurm.sh                    # Run with defaults
#   sbatch submit_slurm.sh "Streptomyces"     # Specify taxon
#   sbatch submit_slurm.sh "Pantoea" bigscape # Specify taxon and clustering
#
# =============================================================================

# ----------------SLURM Parameters----------------
#SBATCH -p normal                    # Partition (adjust for your cluster)
#SBATCH -n 2                         # CPUs for Nextflow head process
#SBATCH --mem=8g                     # Memory for head job
#SBATCH -N 1                         # Single node for head job
#SBATCH --time=48:00:00              # Max wall time
#SBATCH --mail-type=END,FAIL         # Email on completion or failure
#SBATCH -J clusterquest              # Job name
#SBATCH -o clusterquest_%j.out       # Standard output
#SBATCH -e clusterquest_%j.err       # Standard error

# ----------------Configuration-------------------
# Modify these defaults or pass as arguments

TAXON="${1:-Pantoea}"                # First argument or default
CLUSTERING="${2:-bigscape}"          # Second argument or default (none/bigscape)
OUTDIR="results"
RESUME="-resume"                     # Set to "" to start fresh

# ----------------Environment Setup---------------
# Uncomment and modify for your cluster

# Option 1: Load modules (cluster-specific)
# module load nextflow/24.10.4
# module load anaconda3

# Option 2: Use conda environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate nextflow

# Nextflow configuration
export NXF_HOME="${PWD}/.nextflow"
export NXF_DISABLE_CHECK_LATEST=true
export NXF_ANSI_LOG=false            # Disable ANSI colors in log files

# ----------------Validation----------------------
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: nextflow not found. Load the module or activate conda environment."
    exit 1
fi

# ----------------Run Pipeline--------------------
echo "=============================================="
echo "ClusterQuest Pipeline"
echo "=============================================="
echo "Start time:  $(date)"
echo "Taxon:       ${TAXON}"
echo "Clustering:  ${CLUSTERING}"
echo "Output:      ${OUTDIR}"
echo "Job ID:      ${SLURM_JOB_ID:-local}"
echo "Node:        ${SLURMD_NODENAME:-$(hostname)}"
echo "=============================================="
echo ""

nextflow run main.nf \
    -profile slurm \
    --taxon "${TAXON}" \
    --clustering "${CLUSTERING}" \
    --outdir "${OUTDIR}" \
    ${RESUME}

EXIT_CODE=$?

# ----------------Summary-------------------------
echo ""
echo "=============================================="
echo "Pipeline Complete"
echo "=============================================="
echo "End time:    $(date)"
echo "Exit code:   ${EXIT_CODE}"

# The directory the pipeline names for this taxon, sanitised exactly as
# nextflow.config does. Derived once and reused below, so one file cannot carry two
# different taxon-to-directory transforms (the report check used `tr \' \' \'_\'`,
# which differs for any other punctuation).
TAXON_DIR=$(printf '%s' "${TAXON}" \
            | sed -e 's/[^a-zA-Z0-9_]/_/g' -e 's/__*/_/g' -e 's/^_//' -e 's/_$//')

# The trace is at pipeline_info/<Taxon>/<run stamp>/pipeline_trace.tsv. This block
# tested ${OUTDIR}/pipeline_info/pipeline_trace.tsv, which stopped existing when the
# directory gained a taxon level -- so the summary has been omitting the task count
# ever since, as a quietly failing `if` rather than an error. Scoped to THIS job's
# taxon: an unscoped glob would report another taxon's counts under this header,
# newest-first across the whole outdir. Both layouts, so older outdirs still read.
TRACE=$(ls -1t "${OUTDIR}/pipeline_info/${TAXON_DIR}"/*/pipeline_trace.tsv \
               "${OUTDIR}/pipeline_info/${TAXON_DIR}"/pipeline_trace.tsv 2>/dev/null | head -1)
if [ -n "${TRACE}" ]; then
    # CACHED counts as completed. After a -resume nearly every row is CACHED, and
    # `grep -c COMPLETED` counted those as neither completed nor failed, so the
    # numbers did not add up. Status is field 4 of the configured trace fields.
    awk -F'\t' 'NR>1 {
        n++
        if ($4 == "COMPLETED" || $4 == "CACHED") ok++; else bad++
    } END {
        printf "Tasks:       %d/%d completed, %d not complete\n", ok, n, bad
    }' "${TRACE}"
    echo "Trace:       ${TRACE}"
fi

if [ -f "${OUTDIR}/main_analysis_results/${TAXON_DIR}/main_data_visualization/bgc_report.html" ]; then
    echo "Report:      ${OUTDIR}/main_analysis_results/*/main_data_visualization/bgc_report.html"
fi

echo "=============================================="

exit ${EXIT_CODE}
