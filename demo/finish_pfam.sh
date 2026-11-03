#!/usr/bin/env bash
# Unpack, press, and install Pfam 38.2 into the demo database cache.
# Separate from run.sh because hmmpress takes a few minutes and only needs doing once.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SRC=.demo_dbstage/pfam38.2
DST=results_demo/databases/pfam38.2
HMMPRESS=$(ls "$HOME"/.clusterquest/conda/env-*/bin/hmmpress 2>/dev/null | head -1)
[ -n "$HMMPRESS" ] || { echo "no hmmpress found in the conda cache" >&2; exit 1; }

if [ ! -f "$SRC/Pfam-A.hmm" ]; then
    [ -f "$SRC/Pfam-A.hmm.gz" ] || { echo "$SRC/Pfam-A.hmm.gz not found -- download it first" >&2; exit 1; }
    # A partially-downloaded file unzips to a truncated HMM that hmmpress accepts
    # and BiG-SCAPE then fails on, so verify before spending time on it.
    echo "verifying the archive..."
    gzip -t "$SRC/Pfam-A.hmm.gz"
    echo "decompressing..."
    gunzip -k "$SRC/Pfam-A.hmm.gz"
fi

echo "pressing (a few minutes)..."
rm -f "$SRC"/Pfam-A.hmm.h3?
"$HMMPRESS" "$SRC/Pfam-A.hmm"

mkdir -p "$DST"
cp -f "$SRC"/Pfam-A.hmm "$SRC"/Pfam-A.hmm.h3? "$DST"/
echo "pfam=38.2" > "$DST/.db_version"

echo
echo "installed into $DST:"
ls -lh "$DST"
