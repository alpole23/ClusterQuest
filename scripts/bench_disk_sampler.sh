#!/usr/bin/env bash
# Sample disk and memory while a run is in flight.
#
# PEAK is the number the storage work claims to have fixed, and `du` after a run
# measures the residue instead -- a different, smaller quantity. Peak only exists
# while the run is happening, so it has to be sampled. Every projection of it so
# far has been component arithmetic; this is the observation.
#
#   bash scripts/bench_disk_sampler.sh <outdir> <logfile> [interval_s] [workdir] [taxon]
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

printf 'epoch\telapsed_s\tfree_gb\tout_gb\twork_gb\tmem_used_gb\tmem_avail_gb\tswap_used_gb\ttasks\tgenomes_screened\tgenomes_kept\n' > "$LOG"
START=$(date +%s)
while true; do
    now=$(date +%s)
    free_gb=$(df -BG --output=avail "$OUT" | tail -1 | tr -dc '0-9')
    out_gb=$(du -sBG --exclude=databases "$OUT" 2>/dev/null | cut -f1 | tr -dc '0-9')
    work_gb=$(du -sBG --exclude=conda "$WORK" 2>/dev/null | cut -f1 | tr -dc '0-9')
    read -r used avail < <(free -g | awk '/^Mem:/{print $3, $7}')
    swap=$(free -g | awk '/^Swap:/{print $3}')
    tasks=$(find "$WORK" -maxdepth 3 -name '.command.begin' 2>/dev/null | wc -l)
    screened=$(cat $PRESCREEN/prescreen_*.tsv 2>/dev/null | grep -cve '^genome' || true)
    kept=$(find $RENAMED -name '*.gbff' 2>/dev/null | wc -l)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$now" "$((now-START))" "${free_gb:-0}" "${out_gb:-0}" "${work_gb:-0}" \
        "${used:-0}" "${avail:-0}" "${swap:-0}" "$tasks" "${screened:-0}" "${kept:-0}" >> "$LOG"
    sleep "$IV"
done
