#!/usr/bin/env bash
# Snapshot the Nextflow trace periodically, so a restart cannot take the only copy
# of a session's per-task timings with it.
#
# This is not hypothetical. The GenBank Enterobacterales run restarted five times;
# only 24,295 of 53,980 fetch tasks survived in ad-hoc archives, so its total CPU-h
# had to be scaled from 45% coverage. Snapshots are cheap (a few MB) and the
# alternative is remembering to copy the file before every unplanned restart.
#
#   bash scripts/bench_trace_archive.sh <trace-glob> <archive_dir> [interval_s]
#
# TAKES A GLOB, NOT A PATH, and re-resolves it on every tick. Since 2026-11-03 the
# trace lives at pipeline_info/<Taxon>/<UTC run stamp>/pipeline_trace.tsv, and the
# stamp is minted by Nextflow when the config is evaluated -- so no launcher can
# know the path before the run starts. The previous version took a fixed path and
# did `if [ -f "$TRACE" ]`, which on a miss looped forever in silence; the 0-byte
# archiver logs under bench/ are what that looks like from outside. Quote the glob
# so the shell does not expand it at launch:
#
#   bash scripts/bench_trace_archive.sh \
#        '/mnt/bigdisk/clusterquest/results/pipeline_info/Erwiniaceae/*/pipeline_trace.tsv' \
#        /mnt/bigdisk/clusterquest/bench/trace_archive 600
#
# A literal path still works -- it is a glob that matches one file.
#
# Snapshots only when the content has changed, so an idle run costs one md5 per
# interval. Deduplicate at analysis time by task hash -- the same task may appear in
# several snapshots, and a resumed task appears as CACHED in later ones.
#
# DO NOT SUM ACROSS SNAPSHOTS FROM DIFFERENT RUNS. A CACHED row carries the realtime
# of the session that actually ran the task, so the same task recurs in every later
# session's trace; adding them up multiplies the real work by the number of resumes.
set -u
GLOB=${1:?trace glob}; DIR=${2:?archive dir}; IV=${3:-600}
mkdir -p "$DIR"

# Identity of the run a resolved trace belongs to, so snapshots group by run and the
# retention cap below can be per-run. <Taxon>_<stamp> for the current layout,
# <Taxon> for one written before 2026-11-03.
src_id() {
    local p=$1 id
    case $p in
        */pipeline_info/*) id=${p#*/pipeline_info/}; id=${id%/pipeline_trace.tsv} ;;
        *)                 id=$(basename "$(dirname "$p")") ;;
    esac
    printf '%s' "${id//\//_}" | sed -e 's/[^a-zA-Z0-9_]/_/g'
}

last=""; misses=0
while true; do
    # Unquoted on purpose: $GLOB must glob. Newest first, so a resume is followed.
    TRACE=$(ls -1t $GLOB 2>/dev/null | head -1)
    if [ -n "$TRACE" ] && [ -f "$TRACE" ]; then
        misses=0
        cur=$(md5sum "$TRACE" | cut -d' ' -f1)
        if [ "$cur" != "$last" ]; then
            id=$(src_id "$TRACE")
            cp "$TRACE" "$DIR/trace_${id}_$(date -u +%Y%m%d_%H%M%S).tsv"
            last=$cur
            # Keep the 40 most recent FOR THIS RUN. A cap across all runs would let
            # an actively-snapshotting run evict an earlier run's only copy, which
            # is the same loss this script exists to prevent.
            ls -1t "$DIR/trace_${id}_"*.tsv 2>/dev/null | tail -n +41 | xargs -r rm -f
        fi
    else
        # Say so rather than looping in silence, but only once per dry spell: being
        # started before the run is the normal case, not an error.
        misses=$((misses + 1))
        if [ "$misses" = 3 ]; then
            echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) no trace matching $GLOB yet" >&2
        fi
    fi
    sleep "$IV"
done
