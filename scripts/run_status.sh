#!/usr/bin/env bash
# Status of a long ClusterQuest run: is it alive, where is it, is anything wrong.
#
#   bash scripts/run_status.sh [bigdisk_root] [taxon] [run.log] [samples.tsv]
#
# Taxon, log and sample file all DEFAULT TO THE MOST RECENTLY MODIFIED ones rather
# than to a fixed name. They used to be hardcoded to the Enterobacterales run, which
# meant that during the 2026-10-19 Bacteria re-run this script reported that run's
# PID beside Enterobacterales' stages, trace and disk samples -- 1,349,525 genomes
# screened for a run that had screened 22,875. A status tool that mixes two runs is
# worse than none, which is the same lesson as the pgrep note below.
#
# Written because the ad-hoc checks kept lying. `pgrep -f nextflow` matches the
# shell that invokes it -- the pattern appears in that shell's own command line --
# so a check for a dead run reported ALIVE with the checking shell's PID, and a
# 1h-timeout failure went unnoticed for 24 minutes. The `[x]` bracket trick does
# not fix this either: when the surrounding command line mentions the real script
# name, as a launch command does, the bracketed pattern still matches it.
#
# So: identify processes by their java jar / script path via `ps -eo`, never pgrep,
# and never trust a PID recorded at launch time -- `pgrep | head -1` right after a
# setsid can capture a transient wrapper rather than the process that outlives it.
set -u
ROOT=${1:-/mnt/bigdisk/clusterquest}
newest() { ls -1t $1 2>/dev/null | head -1; }
TAXON=${2:-}
if [ -z "$TAXON" ]; then
    T=$(newest "$ROOT/results/pipeline_info/*/pipeline_trace.tsv")
    TAXON=$(basename "$(dirname "${T:-/none/unknown}")")
fi
# Newest log that is actually a PIPELINE log: every run prints [WORKDIR] at start,
# where a sampler's or archiver's log never does. Picking purely by mtime chose the
# trace archiver's and left the stage line blank.
nf_log() {
    for f in $(ls -1t $ROOT/*.log $ROOT/bench/*.log 2>/dev/null); do
        grep -ql '\[WORKDIR\]' "$f" 2>/dev/null && { echo "$f"; return; }
    done
}
LOG=${3:-$(nf_log)}
SAMPLES=${4:-$(newest "$ROOT/bench/*disk*.tsv")}
echo "run        taxon=$TAXON"
echo "           log=${LOG:-none}"
echo "           samples=${SAMPLES:-none}"
strip() { tr -d '\r' | sed 's/\x1b\[[0-9;]*[A-Za-z]//g'; }

NF_PID=$(ps -eo pid,cmd | grep '[n]extflow-.*-one\.jar run main\.nf' | awk '{print $1}' | head -1)
SA_PID=$(ps -eo pid,cmd | grep '[b]ench_disk_sampler\.sh' | awk '{print $1}' | head -1)

if [ -n "$NF_PID" ]; then
    echo "nextflow   ALIVE  pid $NF_PID  elapsed $(ps -o etime= -p "$NF_PID" | tr -d ' ')"
else
    echo "nextflow   *** NOT RUNNING ***"
fi
[ -n "$SA_PID" ] && echo "sampler    alive  pid $SA_PID" || echo "sampler    *** NOT RUNNING ***"
mountpoint -q "$ROOT" 2>/dev/null || mountpoint -q "$(dirname "$ROOT")" \
    && echo "mount      ok" || echo "mount      *** NOT MOUNTED ***"

echo "--- stage ---"
# grep -o, not an anchored line match: with -ansi-log false Nextflow concatenates
# several records onto one line ("[skipping] Stored process > ...[PROCESS ab/cdef]
# NAME (tag)"), so anchoring to ^ printed nothing while the run was progressing
# normally -- a status tool silently reporting no stage.
strip < "$LOG" 2>/dev/null | grep -oE '\[PROCESS [^]]*\] [A-Za-z_:]+( \([^)]*\))?' | tail -4

# A time-limit kill writes exit 143 and Nextflow may abort the session while
# reporting `failed=0`, so the exit codes on disk are checked rather than the
# summary line. 143 = SIGTERM, 137 = SIGKILL, 1 = ordinary failure.
#
# Read from the TRACE, not by walking the work directory. An earlier version ran
# `find work -maxdepth 3 -name .exitcode` and stat'd every task directory, which is
# O(tasks) against a disk the run itself is saturating: at 5,685 tasks it took over
# seven minutes and would only get worse, so the monitoring degraded exactly as the
# run grew big enough to need it. The trace is one sequential read, is maintained by
# Nextflow, and covers EVERY task rather than those a timestamp filter happened to
# keep -- which also retires the stale-task problem the filter existed for, since a
# previous run's tasks are not in this run's trace at all.
echo "--- task failures (from the trace) ---"
TRACE="$ROOT/results/pipeline_info/$TAXON/pipeline_trace.tsv"
if [ -f "$TRACE" ]; then
    awk -F'\t' 'NR>1 {
        n++
        if ($4 != "COMPLETED" && $4 != "CACHED") { bad++; if (bad<=5) print "  " $4 "  exit " $5 "  " $3 }
        if ($5 != "0" && $5 != "-" && $5 != "")  { nz++;  if (nz<=5)  print "  exit " $5 "  " $3 }
    } END {
        printf "  %d tasks, %d not complete, %d non-zero exit\n", n, bad, nz
    }' "$TRACE"
else
    echo "  (no trace at $TRACE)"
fi

echo "--- errors in log ---"
grep -iE 'ERROR|Exception|hasn.t exited|Terminated' "$LOG" 2>/dev/null | strip | head -5
echo "  (nothing above = clean)"

echo "--- latest sample ---"
if [ -f "$SAMPLES" ]; then
    head -1 "$SAMPLES" | tr '\t' '|'
    tail -1 "$SAMPLES" | tr '\t' '|'
    echo "  rows: $(($(wc -l < "$SAMPLES") - 1))"
fi
