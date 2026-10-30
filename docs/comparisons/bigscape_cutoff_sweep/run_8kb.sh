#!/usr/bin/env bash
# Answer "what would the families look like at an 8 kb neighbourhood" without
# re-running antiSMASH. Needs Pfam-A.hmm, which lives on the external drive; set
# PFAM to it and run from this directory.
#
# Detection does not move with the neighbourhood -- the same 334 regions in the same
# 299 genomes were found at 5 kb and 10 kb on Erwiniaceae -- so re-cutting the regions
# antiSMASH already found is equivalent to re-running it, and costs minutes.
set -euo pipefail

AS=${AS:-/home/alexp/enterobacterales/results/antismash_results/Enterobacterales}
PFAM=${PFAM:?set PFAM to Pfam-A.hmm (release 38.2, as nextflow.config pins)}
KB=${KB:-8}
CORES=${CORES:-16}
OUT=${OUT:-nb${KB}kb}

export PYTHONHASHSEED=0          # see docs/comparisons/bigscape_determinism/

python3 trim_regions.py --in "$AS" --out "${OUT}_input" --kb "$KB"

bigscape cluster \
    -i "${OUT}_input" \
    -o "$OUT" \
    --pfam-path "$PFAM" \
    --alignment-mode auto \
    --gcf-cutoffs 0.2,0.3,0.4 \
    --include-singletons \
    --classify category \
    --cores "$CORES"

echo
echo "Compare against the 10 kb run with:"
echo "  python3 compare_runs.py $OUT/$OUT.db \\"
echo "      /home/alexp/enterobacterales/results/bigscape_results/Enterobacterales/Enterobacterales.db"
