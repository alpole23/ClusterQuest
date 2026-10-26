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

### Measured 2026-10-19: dereplication is NOT the answer, and neither is scoping

Three measurements, in the order they were taken, each one killing the proposal before it.

**1. Dereplicating the pool buys 1.37x.** `diamond linclust --approx-id 90` over the
19,818,123-protein pool: **14,496,200 clusters**, a 27% reduction, for 880 s wall /
1.6 CPU-h / 9.86 GB peak RSS.

The reason is the input. The pilot ran `assembly_reference_only: true` -- ONE genome per
species -- so the set is already dereplicated at the organism level, which is exactly the
redundancy 90% clustering removes. The proposal was to collapse "the same protein from 400
strains of one species" in a set containing no such strains. A full RefSeq run would be far
more strain-redundant and would dereplicate better, but by how much is unmeasured, and what
matters is redundancy among *retained* genomes rather than among all of them.

**2. It does not saturate**, which was the hypothesis that would have justified it anyway:

| genomes | clusters | new clusters per added genome |
|---:|---:|---:|
| 10 | 46,907 | 4,391 |
| 100 | 439,542 | 4,610 |
| 1,000 | 4,037,228 | 3,855 |
| 2,000 | 7,633,160 | 3,457 |
| 4,193 | 14,496,200 | **2,971** |

A 4x increase in genomes drops the marginal cost only 23%. Every new species still brings
~3,000 genuinely novel proteins at 90% identity, so on taxonomically diverse input the
dereplicated pool stays essentially linear in N and the quadratic survives. A lower
threshold might saturate, but 50% sits at the `MIN_PCT = 40.0` search floor and would be
unsafe; 7.0% of multi-member clusters already carry disagreeing product names at 90%
(160,338 of 2,282,841), and the representative's name is what lands in the recovered
gene's GFF3 note.

**3. Scoping the query to pepM neighbourhoods buys 11%.** A well-annotated 4.79 Mb genome
emits **701,071 bp across 1,125 gaps -- 14.6% of itself**. Only genes near pepM can reach a
region (rule core plus a 10 kb neighbourhood), so a +/-50 kb window should have cut the
query ~46x. Measured against the full pool:

| query | wall | CPU |
|---|---:|---:|
| 701 kb | 40.39 s | 118.8 CPU-s |
| 15 kb (2.2%) | 35.98 s | 99.7 CPU-s |

**blastx here is almost all fixed cost: ~99 CPU-s to load a 7.19 GB database, ~0.028
CPU-s per kb of query.** Shrinking the query 46x touches the small term.

### What the lever actually is: one blastx per BATCH, not per genome

diamond pays that ~99 CPU-s once per invocation, and `modules/genome/recover_orfs.nf`
invokes it once per genome -- 200 times a batch, against the same database. Concatenating
every genome's gaps into one query and splitting the output by query id amortises the load
across the batch instead of paying it 200 times.

Query ids are already `{rid}|{start}-{end}` and carry the source record; prefixing the
genome name makes them unambiguous across a batch. prodigal stays per-genome -- it is 2.4%
of the cost and trains its model on the whole genome.

Projected per genome, 4,193 retained:

| | CPU-s/genome | pilot CPU-h |
|---|---:|---:|
| as originally written | 600.5 | 699 |
| database built once (shipped) | 126 | 147 |
| + one blastx per batch | ~23 | ~27 |
| + pepM-window scoping on top | ~3.5 | ~4 |

Note the ordering. Scoping alone is 11%; scoping *after* batching is worth ~6x, because
once the fixed cost is amortised the marginal term is all that is left. That inversion is
only visible once the fixed/marginal split is measured rather than assumed.

**Caveat on those projections: the split is fitted on TWO points.** Two points cannot
separate a fixed cost from a linear one -- this repo has made exactly that error before,
with GTDB-Tk, where 8 and 64 genomes were used to argue O(1) and the two candidate models
differed by 94x at a million. Take a third and fourth query size before building a plan on
the batched figure.

### Measured 2026-10-25, and shipped

The caveat above was taken up rather than ignored. Three sizes, same box, same pool,
4 threads, idle:

| genomes per blastx | CPU-s/genome | peak RSS |
|---:|---:|---:|
| 1 (what the pipeline did) | 126.7 | 3.46 GB |
| 12 | 16.8 | — |
| 200 (`pepm_prescreen_batch_size`) | **13.4** | 7.26 GB |

**9.5x, against the ~5.5x the two-point fit projected** -- conservative this time, where
the previous projection in this document was optimistic by 2.23x. Solving the three
points gives ~120 CPU-s fixed per invocation, which is the database load, and the rest
is marginal.

Contention-adjusted, this projected the stage at **~35 CPU-h**: 13.4 scaled by the 2.23x
isolation-to-production factor measured earlier, over 4,193 genomes.

**Measured in production 2026-10-26: 47.93 CPU-h, 41.2 CPU-s/genome** -- 37% over that
projection. The 9.5x isolated ratio held exactly; what did not was the contention factor,
which came in at **3.07x** rather than 2.23x, because 4 concurrent tasks x 4 CPUs now
fully subscribes all 16 cores. That is the THIRD time in this document an idle-box
measurement has under-predicted production. **Multiply an isolated figure by ~3 before
quoting it for this machine**, or measure under load.

The pilot comparison that used to sit here -- "327.6 CPU-h -> ~35" -- is withdrawn with
the rest of the pilot's cost table; see `bacteria_reference_pilot/`. The controlled
before/after above is the defensible evidence for what batching bought, because it holds
the genomes, pool, box and load constant.

**Equivalence, which is the point rather than the speed.** This stage feeds antiSMASH's
gene calls, so a silent difference propagates into every BGC boundary downstream:
13 genomes gave **13/13 byte-identical GFF3s** against the per-genome path, and the same
12 genomes compared against the **200-genome** batch gave 12/12. 713 recovered genes
either way, same homology/prodigal split.

**Memory was the thing that nearly stopped it, and the obvious fixes were wrong.** The
batch peaks at 7.26 GB against `process_medium`'s 8 GB. It is diamond, driven by QUERY
size and not the database -- 0.35 MB of gaps 3.46 GB, 16 MB 7.16, 64 MB 7.26 -- so it
**plateaus** well before a full batch. An earlier reading of the first two of those
points called it "~60 MB of RAM per MB of gaps", which is the two-point error again, in
the same document that warns about it.

Because it plateaus, both ways of bounding it by reshaping the query fail, and both
**change the hits**:

| lever | MAXRSS | CPU | hits vs unsplit |
|---|---:|---:|---|
| split the batch into 4 | 7.16 GB -- **no saving** | 2.1x | +106 / -47, 78 above the 40% floor |
| `diamond -b 0.02` | 0.33 GB | **5.5x** | +13 |
| `diamond -b 0.01` | 0.45 GB | 4.8x | +72 |

So the allocation moves instead: `RECOVER_ORFS` is declared at 12 GB. One unsplit call
is the only variant verified to reproduce the per-genome path byte for byte, and a
memory fix that changes the answer is not a fix.

Note `-M/--memory-limit` is **not available for blastx** in diamond 2.2.6 -- it is a
cluster-workflow option and errors with `Option is not permitted for this workflow`.

**Separate defect found while verifying, not fixed with it.** `run_blastx` asks diamond
for `sseqid`, which BLAST and DIAMOND truncate at the first whitespace, while the pool's
headers are `>locus_tag|product with spaces`. Every homology-recovered gene therefore
carries a one-word product: `DNA polymerase III subunit beta` is written as `DNA`,
`IS3 family transposase` as `IS3`. 290 of 713 genes in the equivalence test.

Fixed in [#14](https://github.com/alpole23/ClusterQuest/pull/14) by asking for `stitle`
instead -- one field, kept out of the batching change because that change's acceptance
test was byte-identical output. Verified on the same 12 genomes against the pre-fix
output: **12/12 identical ignoring `product=`**, 713 recovered genes either way,
products containing a space **0 of 290 -> 280 of 290** (the 10 that stay single-word
genuinely are, e.g. `transposase`).

**Anything produced before that carries the truncated names**, the pilot's recovered
ORFs included -- which is why it is worth landing before a re-run rather than after.

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
