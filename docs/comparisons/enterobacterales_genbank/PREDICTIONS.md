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
