#!/usr/bin/env bash
# Pipeline test suite. Run from anywhere:  bash tests/run_tests.sh
#
# Everything runs in a scratch directory, never in the project's results/, so the
# tests cannot clobber pipeline_info/ (trace, report, timeline all set overwrite=true).
set -uo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRATCH="${TEST_SCRATCH:-$(mktemp -d)}"
FIXTURES="$SCRATCH/fixtures"
BATCH_SIZE="${BATCH_SIZE:-3}"

# Pin the interpreter for the static checks before PATH is modified below, so they
# don't accidentally run under a borrowed conda python of a different version.
SYS_PY="$(command -v python3)"

PASS=0
FAIL=0

pass() { echo "  PASS  $1"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL  $1"; FAIL=$((FAIL + 1)); }
skip() { echo "  SKIP  $1 ($2)"; }
check() { if [ "$2" = "$3" ]; then pass "$1 ($2)"; else fail "$1 (got '$2', expected '$3')"; fi; }

echo "Project:  $PROJECT_DIR"
echo "Scratch:  $SCRATCH"

# GENBANK_TO_FASTA needs biopython. The tests run without conda (creating envs would
# make them slow), so borrow an interpreter that already has it: the current python,
# or one of the pipeline's cached conda envs. Without one, the conversion checks are
# skipped rather than reported as failures.
BIO_PY=""
if python3 -c "import Bio" 2>/dev/null; then
    BIO_PY="$(command -v python3)"
else
    for candidate in "$PROJECT_DIR"/work/conda/env-*/bin/python; do
        if [ -x "$candidate" ] && "$candidate" -c "import Bio" 2>/dev/null; then
            BIO_PY="$candidate"; break
        fi
    done
fi
if [ -n "$BIO_PY" ]; then
    export PATH="$(dirname "$BIO_PY"):$PATH"
    echo "Biopython: $BIO_PY"
else
    echo "Biopython: not found — GenBank→FASTA checks will be skipped"
fi

# Reports are disabled so a test run never truncates results/pipeline_info/.
cat > "$SCRATCH/no_reports.config" <<'EOF'
trace    { enabled = false }
report   { enabled = false }
timeline { enabled = false }
EOF

echo
echo "=== Building fixtures ==="
"$SYS_PY" "$PROJECT_DIR/tests/make_fixtures.py" --outdir "$FIXTURES" --n 7 --with-corrupt || exit 1

# Nextflow sets projectDir from the entry script's location, and the modules resolve
# helper scripts as ${projectDir}/scripts/... . Running the tests straight out of
# tests/ would therefore look for tests/scripts/. So stage each test script in the
# scratch dir (with its relative includes rewritten to absolute) and symlink
# scripts/ there — projectDir becomes the scratch dir and the paths resolve.
ln -sfn "$PROJECT_DIR/scripts" "$SCRATCH/scripts"

stage_nf() {
    sed -e "s|'\.\./|'$PROJECT_DIR/|g" "$PROJECT_DIR/tests/$1" > "$SCRATCH/$1"
}

run_nf() {  # run_nf <script> [extra args...]
    local script="$1"; shift
    stage_nf "$script"
    ( cd "$SCRATCH" && nextflow -quiet run "$SCRATCH/$script" \
        -c "$SCRATCH/no_reports.config" \
        -lib "$PROJECT_DIR/lib" \
        -w "$SCRATCH/work" \
        --fixtures "$FIXTURES" \
        --taxon "Test taxon" \
        --outdir "$SCRATCH/out" \
        --task_batch_size "$BATCH_SIZE" \
        "$@" 2>&1 )
}

echo
echo "=== Utils helpers ==="
UTILS_OUT="$(run_nf test_utils.nf)"
echo "$UTILS_OUT" | grep -E '^(PASS|FAIL) ' | sed 's/^/  /'
if echo "$UTILS_OUT" | grep -q '^FAIL '; then
    FAIL=$((FAIL + 1))
else
    PASS=$((PASS + 1))
fi

echo
echo "=== Task batching ==="
BATCH_OUT="$(run_nf test_batching.nf)"
if ! echo "$BATCH_OUT" | grep -q 'SUCCESS\|Succeeded'; then
    echo "$BATCH_OUT" | tail -20
fi

# 7 good genomes + 1 corrupt: all 8 rename, only the 7 parseable ones convert.
count() { echo "$BATCH_OUT" | grep -o "$1" | wc -l | tr -d ' '; }

check "renamed genomes"        "$(count 'RENAMED:')"        "8"
if [ -n "$BIO_PY" ]; then
    check "fasta conversions"      "$(count 'FASTA:')"        "7"
    check "corrupt genome skipped" "$(count 'FASTA: Broken')" "0"
else
    echo "  SKIP  fasta conversions (no biopython available)"
fi
check "reuse dirs copied"      "$(count 'COPIED:')"         "2"

# Assembly ID → mapped name pairing must survive batch staging (genome?.gbff)
MISPAIRED=0
for i in $(seq 1 7); do
    [ -f "$SCRATCH/out/ncbi_genomes/Test_taxon/renamed_genomes/Test_organism_${i}.gbff" ] || MISPAIRED=$((MISPAIRED + 1))
done
check "assembly ID pairing" "$MISPAIRED" "0"

# Copy fidelity: hidden file and nested directory must survive cp -rL
check "hidden file copied" \
    "$([ -f "$SCRATCH/out/antismash_results/Test_taxon/Genome_A/.antismash_meta" ] && echo yes || echo no)" "yes"
check "nested file copied" \
    "$([ -f "$SCRATCH/out/antismash_results/Test_taxon/Genome_A/nested/deep.txt" ] && echo yes || echo no)" "yes"


# --- batch composition is independent of arrival order ---------------------
echo ""
echo "=== Deterministic batching ==="
if command -v nextflow >/dev/null 2>&1; then
    DET=$(cd "$PROJECT_DIR" && nextflow run tests/test_batch_determinism.nf \
        -profile local 2>&1 || true)
    FWD=$(echo "$DET" | grep -oE "FWD<[^>]*>" | sed 's/FWD//' | sort | md5sum)
    REV=$(echo "$DET" | grep -oE "REV<[^>]*>" | sed 's/REV//' | sort | md5sum)
    TFWD=$(echo "$DET" | grep -oE "TFWD<[^>]*>" | sort | sed 's/TFWD//' | md5sum)
    TREV=$(echo "$DET" | grep -oE "TREV<[^>]*>" | sort | sed 's/TREV//' | md5sum)
    RAW=$(echo "$DET" | grep -oE "RAW<[^>]*>" | sed 's/RAW//' | sort | md5sum)
    [ -n "$(echo "$DET" | grep -oE 'FWD<')" ] && pass "batching ran" || fail "batching ran"
    [ "$FWD" = "$REV" ] && pass "file batches order-independent" \
                        || fail "file batches order-independent"
    [ "$TFWD" = "$TREV" ] && pass "tuple batches order-independent" \
                          || fail "tuple batches order-independent"
    # the guard is only meaningful if unsorted collate really does differ
    [ "$FWD" != "$RAW" ] && pass "unsorted collate does differ (guard is live)" \
                         || fail "unsorted collate does differ (guard is live)"
else
    skip "deterministic batching" "nextflow not available"
fi

# --- GFF3 pairing keeps dotted genome names --------------------------------
echo ""
echo "=== Recovered-ORF GFF3 pairing ==="
if command -v nextflow >/dev/null 2>&1; then
    PAIR_OUT=$(cd "$PROJECT_DIR" && nextflow run tests/test_pairing.nf -profile local 2>&1 || true)
    for NAME in "Pantoea_GCA_963520565.1" "Simple_name" "A.b.c" "Buchnera_B.tra"; do
        if echo "$PAIR_OUT" | grep -qF "PAIRED<${NAME}>"; then
            pass "paired ${NAME}"
        else
            fail "paired ${NAME}"
        fi
    done
else
    skip "GFF3 pairing" "nextflow not available"
fi

echo
echo "=== Python import / syntax check ($("$SYS_PY" --version 2>&1)) ==="
if "$SYS_PY" -m compileall -q "$PROJECT_DIR/scripts" > /dev/null 2>&1; then
    pass "all scripts compile"
else
    fail "compileall reported errors"
fi

UNDEF="$("$SYS_PY" "$PROJECT_DIR/tests/check_undefined.py" "$PROJECT_DIR/scripts")"
if [ -z "$UNDEF" ]; then
    pass "no undefined function calls"
else
    fail "undefined calls found:"; echo "$UNDEF" | sed 's/^/      /'
fi

# Declared scripts-version dependencies must still cover each script's real imports,
# otherwise -resume silently reuses output built from changed code.
if DEPS="$("$SYS_PY" "$PROJECT_DIR/tests/check_script_deps.py" 2>&1)"; then
    pass "script dependency declarations current"
else
    fail "script dependency declarations stale:"; echo "$DEPS" | sed 's/^/      /'
fi

# The report's JavaScript lives in Python string constants, so check_undefined.py
# cannot see it — that gap shipped a search box wired to an undefined function.
if JSCHK="$("$SYS_PY" "$PROJECT_DIR/tests/check_report_js.py" 2>&1)"; then
    pass "report JS handler checks"
else
    fail "report JS handler checks failed:"; echo "$JSCHK" | sed 's/^/      /'
fi

echo
echo "================================"
echo "  passed: $PASS   failed: $FAIL"
echo "================================"
if [ -z "${TEST_SCRATCH:-}" ]; then rm -rf "$SCRATCH"; fi
[ "$FAIL" -eq 0 ]
