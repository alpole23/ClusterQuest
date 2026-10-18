#!/usr/bin/env bash
# Snapshot the Nextflow trace periodically, because the trace scope sets
# `overwrite = true` and every restart destroys the previous session's file.
#
# This is not hypothetical. The GenBank Enterobacterales run restarted five times;
# only 24,295 of 53,980 fetch tasks survived in ad-hoc archives, so its total CPU-h
# had to be scaled from 45% coverage. Snapshots are cheap (a few MB) and the
# alternative is remembering to copy the file before every unplanned restart.
#
#   bash scripts/bench_trace_archive.sh <trace.tsv> <archive_dir> [interval_s]
#
# Snapshots only when the content has changed, so an idle run costs one stat per
# interval. Deduplicate at analysis time by task hash -- the same task may appear in
# several snapshots, and a resumed task appears as CACHED in later ones.
set -u
TRACE=${1:?trace}; DIR=${2:?archive dir}; IV=${3:-600}
mkdir -p "$DIR"
last=""
while true; do
    if [ -f "$TRACE" ]; then
        cur=$(md5sum "$TRACE" | cut -d' ' -f1)
        if [ "$cur" != "$last" ]; then
            cp "$TRACE" "$DIR/trace_$(date +%Y%m%d_%H%M%S).tsv"
            last=$cur
            # keep the 40 most recent; they are supersets of each other within a
            # session, so old ones only matter across restarts
            ls -1t "$DIR"/trace_*.tsv 2>/dev/null | tail -n +41 | xargs -r rm -f
        fi
    fi
    sleep "$IV"
done
