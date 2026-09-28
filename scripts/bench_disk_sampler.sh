#!/usr/bin/env bash
# Sample disk and memory while a run is in flight.
#
# PEAK is the number the storage work claims to have fixed, and `du` after a run
# measures the residue instead -- a different, smaller quantity. Peak only exists
# while the run is happening, so it has to be sampled. Every projection of it so
# far has been component arithmetic; this is the observation.
#
#   bash scripts/bench_disk_sampler.sh <outdir> <logfile> [interval_s]
set -u
OUT=${1:?outdir}; LOG=${2:?logfile}; IV=${3:-60}
printf 'epoch\telapsed_s\tfree_gb\tout_gb\twork_gb\tmem_used_gb\tmem_avail_gb\ttasks\n' > "$LOG"
START=$(date +%s)
while true; do
    now=$(date +%s)
    free_gb=$(df -BG --output=avail . | tail -1 | tr -dc '0-9')
    out_gb=$(du -sBG --exclude=databases "$OUT" 2>/dev/null | cut -f1 | tr -dc '0-9')
    work_gb=$(du -sBG --exclude=conda work 2>/dev/null | cut -f1 | tr -dc '0-9')
    read -r used avail < <(free -g | awk '/^Mem:/{print $3, $7}')
    tasks=$(find work -maxdepth 2 -name '.command.begin' 2>/dev/null | wc -l)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$now" "$((now-START))" "${free_gb:-0}" "${out_gb:-0}" "${work_gb:-0}" \
        "${used:-0}" "${avail:-0}" "$tasks" >> "$LOG"
    sleep "$IV"
done
