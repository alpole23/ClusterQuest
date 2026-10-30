#!/usr/bin/env bash
# What do the families look like under a narrower comparison window?
#
# Three arms, because the interesting answer is not "8 kb" but "asymmetric":
#   sym8   core +/- 8 kb      -- the symmetric shrink
#   asym   core -0 / +10 kb   -- the whole cassette, none of the upstream flank
#   asym2  core -1 / +10 kb   -- as above with a little upstream slack
#
# Measured on three independently split family pairs, the phosphonate cassette lies
# ENTIRELY downstream of the rule core (median cross-family Jaccard 1.00 from +0 to
# +7 kb) while every upstream kilobase is divergent host context (0.00). So a
# symmetric shrink is the worst of both: it clips the terminal ATP-grasp, which sits
# 9,189-9,792 bp past the core, AND keeps 8 kb of the flank that splits on species.
#
# Needs Pfam-A.hmm, which lives on the external drive. Run from this directory.
set -euo pipefail

AS=${AS:-/home/alexp/enterobacterales/results/antismash_results/Enterobacterales}
PFAM=${PFAM:?set PFAM to Pfam-A.hmm (release 38.2, as nextflow.config pins)}
CORES=${CORES:-16}
REF=${REF:-/home/alexp/enterobacterales/results/bigscape_results/Enterobacterales/Enterobacterales.db}

export PYTHONHASHSEED=0          # see docs/comparisons/bigscape_determinism/

run_arm () {
    local name=$1; shift
    echo "=== $name ==="
    python3 trim_regions.py --in "$AS" --out "${name}_input" "$@"
    bigscape cluster -i "${name}_input" -o "$name" --pfam-path "$PFAM" \
        --alignment-mode auto --gcf-cutoffs 0.2,0.3,0.4 --include-singletons \
        --classify category --cores "$CORES"
    python3 compare_runs.py "$name/$name.db" "$REF"
}

run_arm sym8  --kb 8
run_arm asym  --up 0 --down 10
run_arm asym2 --up 1 --down 10

echo
echo "Read the asym arm against sym8: if asym merges the species-split families and"
echo "sym8 does not, the flank is upstream and the window should be asymmetric."
