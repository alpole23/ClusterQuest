# GenBank Enterobacterales — predictions registered before the run finished

Written 2026-10-11, with the run ~20 min in and still in `NCBI_FETCH_METADATA`. Every
projection in `CLAUDE.md` so far was written up *after* the measurement that corrected it,
which makes it impossible to tell a model that worked from one that was retrofitted. These
are the numbers to be judged against.

**This run is larger than the 1 M-genome extrapolation point.** NCBI returned
**1,349,507** genome records for the order, so the "projected to 1M genomes" figures stop
being an extrapolation here and become a direct measurement past it.

## Measured, not predicted

| | |
|---|---:|
| genomes offered by NCBI | **1,349,507** |
| scale vs the RefSeq run (150,690) | 8.956x |
| box | 16 cores, 56 GB, ext4 on USB3 |
| link, measured 2026-10-11 | 51.6 MB/s single, 72.7 MB/s over 4 forks |

## Predictions

Stated with the reasoning, because a number without one cannot be wrong in an
interesting way.

| quantity | prediction | basis | confidence |
|---|---:|---|---|
| screen retention | **0.40-0.87%** | RefSeq was 0.87%; GenBank's extra 1.2 M are overwhelmingly more *E. coli* / *Salmonella* / *Klebsiella* clinical isolates, so this should fall rather than hold | firm on direction, wide on degree |
| genomes to antiSMASH | **5,400-11,700** | retention x 1,349,507 | follows from the above |
| total CPU-h | **630-920** | per-stage scaling below, not a linear 8.96x | moderate |
| wall clock | **2.5-4 d** | 44.9 h download-bound at 72.7 MB/s against 58.3 h CPU-bound at 16 cores | moderate |
| binding constraint | **CPU, not the link** | the hardwire moved download from ~130 h to ~45 h, below the CPU bound — this inverts the stated bottleneck | firm |
| peak disk | **150-250 GB** | GTDB-Tk's ~95 GB pplacer scratch is a floor independent of input, and dominated RefSeq's 150 GB peak | firm |
| peak memory | **~42 GB** | skani sketching GTDB's 143,614 references, measured at 41.7 GB PSS and input-independent | firm |
| BGC regions | **1,300-11,700** | genuinely unknown; RefSeq gave 1,303 and the clinical bulk may add almost none | **low — this is the open question** |
| BiG-SCAPE RAM | **2-30 GB** | corrected fit `GB = 1.80 + 2.04e-7*n^2` over that BGC range | firm given the BGC count |
| partitioning needed? | **no** | 30 GB at 11,707 BGCs fits 56 GB; `bigscape_partition` left `false` deliberately | moderate |
| task failures | **0** | RefSeq ran 6,098 tasks with 0 failures; the invalid-zip rate is absorbed by inline retries | firm |

### Per-stage CPU-h, scaled by what each stage actually tracks

Not a uniform 8.96x. `FETCH_RENAME_SCREEN` scales with *total* genomes; antiSMASH and
`RECOVER_ORFS` with *retained*; GTDB-Tk with *BGC-positive*; BiG-SCAPE with BGCs squared.

| stage | RefSeq | at 0.40% retention | at 0.87% | scales with |
|---|---:|---:|---:|---|
| `FETCH_RENAME_SCREEN` | 42.1 | 377 | 377 | total genomes |
| `ANTISMASH` | 30.5 | 126 | 274 | retained |
| `RECOVER_ORFS` | 20.9 | 86 | 187 | retained |
| `GTDBTK_CLASSIFY` | 9.1 | 38 | 82 | BGC-positive |
| `BIGSCAPE` | 0.7 | 1 | 1 | BGCs squared |
| **total** | **104.2** | **627** | **920** | |

A naive linear 8.96x gives 933 CPU-h and lands near the top of that range **for the wrong
reason** — it over-counts antiSMASH and under-counts nothing, so agreement would be
coincidence. If the measured total comes in near 930, check the per-stage split before
concluding the linear model worked.

## What will NOT be comparable

The `CLAUDE.md` 1 M projections (antiSMASH ~25,000 CPU-h, GTDB-Tk 563, BiG-SCAPE ~1,090,
~$1,320) assume **SLURM at 200 concurrent jobs and Erwiniaceae-like 11% prevalence**. Both
are wrong for this run: one 16-core box, and a clade at ~0.9%. The comparable axis is
**CPU-h per stage and per retained genome**, not wall clock and not those totals.

## Instrumentation

| what | where |
|---|---|
| per-task CPU, time, RSS, I/O | `results/pipeline_info/Enterobacterales/pipeline_trace.tsv` |
| peak disk, memory, swap, task count, genome counters | `bench/disk_samples.tsv`, 60 s interval |
| analysis | `scripts/benchmark_analysis.py --trace ... --targets ...` |

Three bugs in the sampler were fixed before it was trusted here, all of which reported a
plausible number rather than an error: it measured `./work` and `df .` regardless of where
`-w` pointed; `tasks` counted at `-maxdepth 2` where `.command.begin` sits at depth 3, so
that column read **0 for all 2,607 samples of the RefSeq run**; and the genome counters
summed every taxon in the outdir, picking up a 58-genome smoke run. Treat the RefSeq
`disk_samples.tsv` `tasks` column as absent, not as zero.

---

# Measured, 2026-10-17

Run completed `[SUCCESS] completed=22856 failed=0 cached=31202`. The predictions above
were registered on 2026-10-11 with the run 20 minutes in; nothing below has been
back-edited into them.

## Scorecard

| quantity | predicted | **measured** | |
|---|---|---:|---|
| screen retention | 0.40–0.87% | **0.113%** | ✗ **3.5x below the floor** |
| genomes to antiSMASH | 5,400–11,700 | **1,530** | ✗ follows from the above |
| BGC regions | 1,300–11,700 | **1,476** | ✓ bottom of range |
| total CPU-h | 630–920 | **~649** | ✓ |
| peak disk | 150–250 GB | **330 GB** | ✗ **1.3x over** |
| peak memory | ~42 GB | 8 GB system; 102 GB `peak_rss` on GTDB-Tk, which is the documented instrument artefact | ✓ in substance |
| BiG-SCAPE RAM | 2–30 GB | **19 GB reported, 0.8 CPU-h** | ✓ |
| partitioning needed? | no | **no** | ✓ |
| task failures | 0 | **0** of 54,058, and 0 non-zero exits | ✓ |
| binding constraint | CPU, not the link | **neither — per-task latency, then the USB HDD** | ✗ see below |

## The two misses

**Retention, 0.113% against a 0.40% floor.** The direction was right and argued for the
right reason — GenBank's extra 1.2 M genomes are overwhelmingly the same clinical isolates
that drove RefSeq to 0.87% — but the degree was underestimated by 3.5x. The estimate
stabilised early and never moved: 0.154% at 11 k genomes screened, 0.115% at 500 k,
0.113% at 1.35 M. Anyone sizing a clade like this should weight the clinical fraction far
more heavily than the RefSeq rate suggests.

**Peak disk, 330 GB against a 150–250 GB range.** The error is specific and avoidable:
the projection assumed one `GTDBTK_CLASSIFY` task at ~95 GB of pplacer scratch. **Two
shards ran concurrently**, and this repository's own benchmark notes already said "~95 GB
of scratch is per *concurrent* task, so two shards at once is ~190 GB". The number was on
the page and was not carried into the estimate.

## The binding constraint was neither prediction

The registered note said the hardwire had moved the bottleneck from the link to CPU. In
fact **nothing was saturated**: at 12 concurrent fetch tasks the network ran ~11 MB/s of
~56 available, disk 57.7% busy, CPU ~7 of 16 cores, memory 5 of 56 GB. The stage was
latency-bound — each task waiting in sequence on NCBI round trips, then `sync`, then the
corruption check, then diamond.

Four tuning rounds followed, and three of the four hypotheses were wrong (raw bandwidth,
fork count, the filesystem-wide `sync` — which was measured and exonerated). The real cap
was `memory = 8 GB` inherited from `process_medium`: the local executor admits a task only
while both cpus *and* memory remain free, so 56 GB of RAM pinned concurrency at **seven**
tasks whatever `download_max_forks` said. Measured `peak_rss` was p50 0.23 GB, max 1.20 GB
— 8.5x over-allocated. Fixing it gave 1.40x, after which the USB HDD saturated at 100%
and became the real ceiling. See PR #11.

## Per-stage cost

Actual CPU (`realtime x %cpu`), matching how the RefSeq baseline was computed. Allocation-
based figures (`realtime x cpus`) read 755 CPU-h and are **not** comparable, because the
allocation changed mid-run from 4 cpus to 1.

| stage | CPU-h | share | note |
|---|---:|---:|---|
| `FETCH_RENAME_SCREEN` | ~589 | 91% | scaled from 24,295 of 53,980 batches traced |
| `ANTISMASH` | 37.0 | 5.7% | 31 batches of 50 |
| `RECOVER_ORFS` | 19.3 | 3.0% | 8 batches of 200 |
| `GTDBTK_CLASSIFY` | 1.8 | 0.3% | single shard; wall time ~3 h, almost all pplacer I/O wait |
| `BIGSCAPE` | 0.8 | 0.1% | |

**Caveat on the fetch figure:** trace files are written with `overwrite = true`, so each
restart discarded the previous session's trace. Only 45% of fetch batches survive in the
archived traces and the total is scaled from those. Archive `pipeline_trace.tsv` before
every restart if a complete accounting matters.

## Result

| | **GenBank** | RefSeq (2026-09-27) |
|---|---:|---:|
| genomes offered | 1,349,507 | 150,690 |
| retained | 1,530 (0.113%) | 1,309 (0.87%) |
| BGC regions | **1,476** | 1,303 |
| GCF families | **112** | 72 |
| largest family | 186 | 186 |
| singletons | 43 | 16 |
| wall clock | ~60 h + 12 h outage | 43.1 h |

**8.9x the input gave 13% more BGCs but 56% more families.** The clinical bulk contributes
essentially nothing while the rare carriers are markedly more diverse than RefSeq showed —
40 new families from 173 extra BGCs. The largest family is 186 in both runs: that is
pantaphos, and RefSeq had already saturated it.

Do not read the singleton count (16 -> 43) as diversity without checking cohesion. The
dereplication argument that explained the 640-genome pilot does not apply here — this is
full depth — but `family_cohesion/` is the way to confirm rather than assume.

## Interruptions

One unplanned reboot: Windows Update (`TrustedInstaller.exe`, "Operating System: Upgrade
(Planned)") at 01:33 on 2026-10-14, with the run at 57.8%. Clean shutdown, no Event 41, no
ext4 errors; the journal replayed in 51 s. Cost was ~12 h of wall clock — WSL does not
auto-start — and roughly 100 in-flight genomes. All 31,187 completed batches resumed from
cache. Updates were then paused until Nov 1.
