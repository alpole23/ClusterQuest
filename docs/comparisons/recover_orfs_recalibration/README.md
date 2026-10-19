# Why `RECOVER_ORFS` is slow, and why better annotation makes it worse

Measured 2026-10-19 on the internal NVMe, after the USB dock failure of 10-18 made the
first attempt at this question worthless.

## The retraction first

On 2026-10-18 this project recorded `RECOVER_ORFS` at **"~7.5 min/genome, 10x slower on
diverse input"**, measured during the whole-Bacteria reference pilot, and a revised 3-5
week whole-Bacteria projection was built on it.

**That measurement is void.** All seven `diamond` processes were in uninterruptible-sleep
(D state) on a failing USB-SATA bridge at the time — `UASPStor` Event 129 bus resets plus
41 x `disk` Event 153 "I/O retried", all on 10-18. The number measures a dying disk.

It is withdrawn rather than corrected, and nothing should be rebuilt on it. The Bacteria
run's trace carries **no `RECOVER_ORFS` row at all** — the session aborted on Nextflow's
`IllegalThreadStateException` before writing one — so there is no record to salvage even
if one wanted to.

A clean measurement now gives **11.8x**, which is close to the withdrawn figure. That is a
coincidence worth naming rather than a vindication: the withdrawn number was attributed to
"diverse input being harder to annotate", and the real mechanism is the opposite (below).
Right magnitude, wrong cause, invalid measurement.

## The question

Enterobacterales measured **45.3 CPU-s per retained genome**. The Bacteria pilot's genomes
are NCBI *reference* genomes — RefSeq, PGAP-annotated, about as good as bacterial
annotation gets. Recovery exists to rebuild genes a submission left out, so a well-annotated
genome should give it *less* to do.

It got slower instead. By 11.8x.

## What the stage actually spends its time on

`scripts/genome/recover_orfs.py` runs three external things per genome. Timed by patching
`subprocess.run`, so the script itself is unmodified and this is exactly what the pipeline
runs. One genome (`Abyssisolibacter_fermentans_MCWD3`, 4,216 annotated CDS), 4 threads,
both pools staged on NVMe so the disk cannot confound it:

| phase | Enterobacterales pool (447 MiB) | Bacteria pool (6.89 GiB) | ratio |
|---|---:|---:|---:|
| `prodigal` | 2.35 s | 2.45 s | 1.0x |
| **`diamond makedb`** | 3.21 s | **52.29 s** | **16.3x** |
| `diamond blastx` | 2.51 s | 45.15 s | 18.0x |
| python (parse / gaps / merge / write) | 0.46 s | 0.64 s | 1.4x |
| **total** | **8.53 s · 40.7 CPU-s** | **100.5 s · 600.5 CPU-s** | **11.8x · 14.8x** |

Same genome, same code, same hardware. **The only variable is the protein pool.**

The 40.7 CPU-s reproduces the Enterobacterales run's traced 45.3 CPU-s/genome from first
principles, which is the control: the clean rig agrees with the real run.

## `diamond makedb` runs once per genome

`recover_orfs.py:257`, inside `main()`, inside the per-genome loop that
`modules/genome/recover_orfs.nf` drives with `for GENOME in ${genomes}`:

```python
if n:                      # n = number of intergenic gaps found
    db = a.workdir / 'pool'
    subprocess.run([a.diamond, 'makedb', '--in', str(a.pool), '-d', str(db), '--quiet'],
                   check=True, capture_output=True)
    hom = homology_orfs(a.diamond, gf, db, a.threads, records)
```

The pool is one file for the whole run — `BUILD_PROTEIN_POOL` is a gather step producing a
single `protein_pool.faa`. So every genome rebuilds **the identical 7.19 GB database from
scratch**. A 200-genome batch builds it 200 times.

Two aggravating details:

- **The gate never fires.** `if n:` asks only whether the genome has at least one
  unannotated stretch of `MIN_GAP = 300` bp. Every bacterial genome has intergenic DNA at
  that length — rRNA operons, prophage, repeats, ordinary inter-operon spacing. A pristine
  reference genome pays the full database build to search a handful of small gaps.
- **It takes every core.** `makedb` is called with no `--threads`, unlike the `blastx`
  below it which correctly passes `--threads ${task.cpus}`. Standalone on the Bacteria
  pool: 84.35 s / 268 CPU-s at 4 threads, **58.95 s / 500 CPU-s at 16**. Under Nextflow
  that is 8 concurrent tasks each grabbing 16 cores on a 16-core box, which is what the
  traced `%cpu` of 211-256% against `cpus = 4` is showing.

## Why good annotation makes the pool bigger

This is the answer to the question as asked, and it runs opposite to intuition.

**The pool is built from the retained genomes themselves.** `BUILD_PROTEIN_POOL` pools
every already-annotated protein across everything the screen kept — the design premise
being that screened genomes are each other's best reference for exactly the genes recovery
is trying to rebuild.

So annotation quality is an input to pool size, not a discount on it:

| | Enterobacterales (GenBank) | Bacteria pilot (RefSeq reference) |
|---|---:|---:|
| genomes retained | 1,530 | 4,193 |
| pool | 447 MiB | **6.89 GiB** |
| proteins | 1,177,087 | **19,818,123** |
| **proteins per retained genome** | **769** | **4,726** |

2.74x the genomes, **16.8x the proteins.**

A PGAP-annotated genome carries ~4,700 proteins, and the pilot's figure is 4,726 — every
retained genome donating its full complement. Enterobacterales averages 769, about 16% of
one genome's worth, because most GenBank submissions carry no CDS translations at all.
Measured directly on one recovery batch: **83 of 130 genomes (64%) had no existing
annotation**, logged by the script as `no existing annotation, leaving it to antiSMASH`.

So the Enterobacterales run was fast *because its input was badly annotated*. Three
quarters of its retained genomes contributed nothing to the pool, the pool stayed at
447 MiB, and the per-genome database rebuild stayed cheap enough to hide. Feeding the
stage well-annotated genomes removed the thing that was accidentally protecting it.

The cost is not pure waste. The bigger pool recovered **46 genes against 23** on the same
genome (28 by homology against 5) — a broader reference set genuinely finds more. The
`makedb` half of the cost is waste; the `blastx` half buys something.

## Recalibrated cost

Per retained genome, Bacteria-pilot pool, measured:

| | wall | CPU |
|---|---:|---:|
| as written | 100.5 s | 600.5 CPU-s |
| with the database prebuilt | **46 s** | **126 CPU-s** |

(prebuilt: `blastx` 42.88 s / 123.2 CPU-s + `prodigal` 2.46 s / 2.6 CPU-s + ~0.6 s python)

Over the pilot's 4,193 retained genomes:

| | CPU-h | wall per 200-genome batch |
|---|---:|---:|
| as written | **699** | 5.6 h |
| database built once | **147** | 2.6 h |

For scale, the Enterobacterales run's entire `RECOVER_ORFS` stage was **19.25 CPU-h** over
1,530 genomes (8 tasks, traced), and the whole GenBank run was 649 CPU-h.

The 5.6 h/batch also explains the kill: `process_medium` allows 2 h, so the task died at
exactly 2 h having completed roughly a third of its batch. The 12 h `withName` override
now in `conf/labels.config` is sufficient, but it is treating the symptom.

## The scaling hazard, which is the real finding

**As written the stage is quadratic in retained genomes.** The pool grows with N and the
database is rebuilt N times, so `makedb` alone is O(N x pool) = O(N²). Ten times the
retained genomes is a hundred times the database-building cost.

This was invisible at Enterobacterales scale only because poor GenBank annotation held the
pool to 447 MiB. On well-annotated input the term is fully exposed.

**Building the database once removes the `makedb` term but not the quadratic.** `blastx`
still searches every genome against the whole pool, which still grows with N: 2.51 s
against 1.18 M proteins, 45.15 s against 19.8 M — linear in pool size, N times over.

So a whole-Bacteria run needs the pool bounded, not just the database cached. Options, in
rough order of appeal, none of them measured yet:

- **Dereplicate the pool** (`diamond cluster`, `mmseqs linclust`). 19.8 M proteins drawn
  from 4,193 genomes is enormously redundant, and the homology pass needs one
  representative per protein family, not one per genome.
- **Scope the pool** to a clade or GCF rather than the whole run.
- **Cap it**, sampling representatives once it exceeds a size.

## Fix

Build the DIAMOND database once, in `BUILD_PROTEIN_POOL`, and pass the `.dmnd` to
`RECOVER_ORFS` alongside (or instead of) the FASTA:

1. `BUILD_PROTEIN_POOL` emits `protein_pool.dmnd`; add `bioconda::diamond` to its conda
   spec in `conf/conda.config` (it currently has python + biopython only).
2. `RECOVER_ORFS` takes it as an input and passes `--pool_db`.
3. `recover_orfs.py` skips `makedb` when `--pool_db` is given.
4. Wire the new output at `subworkflows/antismash_analysis.nf:76`.

Nextflow stages inputs as symlinks, so the 7.19 GB database costs one copy on disk however
many tasks consume it.

Measured payoff on the pilot: **2.2x wall, 4.8x CPU, ~552 CPU-h saved.**

## Reproducing

```bash
# the per-call timing harness used above
python docs/comparisons/recover_orfs_recalibration/timed_recover.py --genome <g>.gbff --out /tmp/o.gff3 \
    --pool protein_pool.faa --prodigal $(which prodigal) \
    --diamond $(which diamond) --threads 4 --workdir /tmp/wd
```

Stage both the pool and the genome on local NVMe. The 10-18 measurement was taken with the
pool on a failing USB disk, and the resulting number was wrong by an unknown factor in an
unknown direction — `/proc/diskstats` and `ps -o stat` (look for `D`) before trusting any
timing on external storage.
