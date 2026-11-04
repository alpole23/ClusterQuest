#!/usr/bin/env bash
# ClusterQuest demonstration runner.
#
#   bash demo/run.sh offline     # 24 genomes already on disk, no network
#   bash demo/run.sh live        # 58 Brenneria genomes fetched from NCBI
#
# Everything writes to results_demo/, which is also where the reference
# databases live, so the storeDir cache is hit instead of re-downloading.
set -euo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Nextflow needs a JVM, and the only one on this box is inside the conda env.
export JAVA_HOME="$HOME/miniconda3/envs/nextflow"
export PATH="$PATH:$HOME/miniconda3/envs/nextflow/bin"

# The antiSMASH conda env ships no CA bundle, so its Python cannot verify TLS
# and every database download dies with CERTIFICATE_VERIFY_FAILED. Point it at
# the system store. Harmless when the databases are already cached.
export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt

MODE="${1:-offline}"
case "$MODE" in
    offline) PARAMS=demo/params_offline.json ;;
    live)    PARAMS=demo/params_live_brenneria.json ;;
    *) echo "usage: bash demo/run.sh [offline|live]" >&2; exit 2 ;;
esac

echo "=== ClusterQuest demo: $MODE ==="
echo "params:  $PARAMS"
echo "outdir:  results_demo"
echo

# -resume, but NEVER bare. A bare `-resume` binds to the most recent session in
# .nextflow/history regardless of what that session was -- which on 2026-11-02
# meant this demo attached itself to a live 131 GB GTDB-Tk download and died with
#
#     ERROR ~ Unable to acquire lock on session with ID 5ee43dd8-...
#
# So resolve the last session that ran THIS demo's params file and resume that one
# by UUID. No match (a first run, or a pruned history) means no -resume at all,
# which is correct: a cold demo is ~5.5 min, and silently resuming the wrong
# session is far worse than running again.
SESSION=""
if [ -r .nextflow/history ]; then
    SESSION=$(awk -F'\t' -v p="$PARAMS" '$7 ~ p && $6 ~ /^[0-9a-f-]{36}$/ { s = $6 } END { print s }' \
        .nextflow/history)
fi

if [ -n "$SESSION" ]; then
    echo "resuming session $SESSION (last run of $PARAMS)"
    exec nextflow run main.nf \
        -profile local \
        -params-file "$PARAMS" \
        -resume "$SESSION" \
        "${@:2}"
else
    echo "no previous session for $PARAMS; running cold (~5.5 min)"
    exec nextflow run main.nf \
        -profile local \
        -params-file "$PARAMS" \
        "${@:2}"
fi
