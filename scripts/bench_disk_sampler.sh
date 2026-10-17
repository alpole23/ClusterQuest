#!/usr/bin/env bash
# Sample disk and memory while a run is in flight.
#
# PEAK is the number the storage work claims to have fixed, and `du` after a run
# measures the residue instead -- a different, smaller quantity. Peak only exists
# while the run is happening, so it has to be sampled. Every projection of it so
# far has been component arithmetic; this is the observation.
#
#   bash scripts/bench_disk_sampler.sh <outdir> <logfile> [interval_s] [workdir] [taxon] [batch_size]
#
# Three things here were wrong in the version that produced the RefSeq benchmark,
# and all three failed by reporting a plausible number rather than an error.
#
# The work directory is an ARGUMENT, not `./work`, and free space is read from the
# outdir's filesystem rather than `.`. Point `-w` at another drive and the old
# version sampled the repo's filesystem and a stale work dir from some earlier run.
#
# `tasks` counts at -maxdepth 3, because a task directory is work/XX/XXXXXX and
# .command.begin sits inside it. -maxdepth 2 matched nothing ever: the column read
# 0 for all 2,607 samples of the 43 h RefSeq run and was mistaken for a
# measurement. Treat that run's `tasks` column as absent, not as zero.
#
# The taxon scopes the genome counters. One outdir legitimately holds many taxa --
# that is the point of namespacing by taxon -- so counting across all of them adds
# a previous run's genomes to this one's. A 58-genome smoke run inflated the first
# samples taken here before this argument existed.
set -u
OUT=${1:?outdir}; LOG=${2:?logfile}; IV=${3:-60}; WORK=${4:-work}; TAXON=${5:-}
BATCH=${6:-25}   # download_batch_size, for converting batch count to genome count
[ -d "$OUT" ]  || { echo "outdir not found: $OUT" >&2; exit 1; }
[ -d "$WORK" ] || { echo "workdir not found: $WORK" >&2; exit 1; }
if [ -z "$TAXON" ]; then
    echo "warning: no taxon given; genome counters sum every taxon in $OUT" >&2
    TAXON='*'
fi

# genomes_screened counts prescreen verdict rows -- every genome the fetch stage
# resolved, pass or fail -- and genomes_kept counts the survivors antiSMASH will
# see. Together they give the retention rate and the sustained fetch throughput,
# which is what any wall-time projection actually rests on. Unquoted so a '*'
# taxon still globs.
PRESCREEN="$OUT/prescreen_results/$TAXON"
RENAMED="$OUT/ncbi_genomes/$TAXON"

printf 'epoch\telapsed_s\tfree_gb\tfs_used_gb\tmem_used_gb\tmem_avail_gb\tswap_used_gb\tbatches\tgenomes_screened\tgenomes_kept\n' > "$LOG"
START=$(date +%s)
while true; do
    now=$(date +%s)
    free_gb=$(df -BG --output=avail "$OUT" | tail -1 | tr -dc '0-9')
    # df, not du: whole-filesystem numbers are O(1) where `du` is O(files). Note this
    # is the WHOLE filesystem including the ~163 GB of reference databases, which is
    # why the column is named fs_used_gb and not out_gb -- it is not comparable to the
    # du-based out_gb in samples taken before 2026-10-12. Peak disk is still
    # recoverable from it, as a delta against the first row rather than an absolute.
    fs_used_gb=$(df -BG --output=used "$OUT" | tail -1 | tr -dc '0-9')
    read -r used avail < <(free -g | awk '/^Mem:/{print $3, $7}')
    swap=$(free -g | awk '/^Swap:/{print $3}')
    tasks=$(ls $PRESCREEN 2>/dev/null | wc -l)
    screened=$(( tasks * BATCH ))
    kept=$(ls $RENAMED/renamed_genomes 2>/dev/null | wc -l)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$now" "$((now-START))" "${free_gb:-0}" "${fs_used_gb:-0}" \
        "${used:-0}" "${avail:-0}" "${swap:-0}" "$tasks" "${screened:-0}" "${kept:-0}" >> "$LOG"
    sleep "$IV"
done
