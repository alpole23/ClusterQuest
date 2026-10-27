# Whole-Bacteria reference sweep

One NCBI reference genome per species across all of Bacteria, to answer whether a
domain-wide phosphonate sweep is tractable on one workstation and what it would find.
It is.

Two runs, and the distinction matters for every number below:

| | pilot, 2026-10-19 → 10-25 | **rerun, 2026-10-25 → 10-26** |
|---|---|---|
| outcome | completed, but produced **no report** | completed with the report |
| cost accounting | **withdrawn — see below** | clean |
| code | pre-#12/#13/#14 | post |

**Read the rerun column.** The pilot's published stage costs could not be reproduced
and are retracted; the section "The pilot's cost table is withdrawn" says exactly which
numbers were wrong and how that was established.

## Configuration

```json
{"taxon": "Bacteria", "assembly_source": "RefSeq",
 "assembly_reference_only": true, "outdir": "/mnt/bigdisk/clusterquest/results"}
```

`assembly_reference_only` **requires** `assembly_source = "RefSeq"` — reference
designation exists only on GCF\_ accessions, and the combination with GenBank silently
returns zero genomes. Measured: Enterobacterales GenBank 0, RefSeq 643, Bacteria 22,871.

## Result (rerun, measured)

| | |
|---|---:|
| genomes fetched | 22,875 |
| retained by the pepM screen | **4,193 (18.3%)** |
| BGC regions | **1,345** |
| BGC-positive genomes | 1,264 |
| GCF families | 1,022 (938 singletons, largest 94) |
| GTDB-Tk classified | 1,263 genomes, 24 phyla |
| tasks / failures | 1,056 / **0** |
| wall | 12 h 51 m |
| **CPU, work actually performed** | **173.8 CPU-h** |
| CPU, incl. tasks resumed from the pilot | 266.5 CPU-h |
| peak disk / memory | 512 GB / 20 GB (0 swap) |
| `bgc_report.html` | **47.7 MB**, 11 panes, 4,193 genome pages, 1,022 family pages |

Taxonomic spread of carriers: Pseudomonadota 582, Actinomycetota 271, Bacillota 230,
Bacteroidota 43, Campylobacterota 30, across 24 phyla.

### Cost by stage

The rerun resumed the pilot's session, so the fetch and GTDB-Tk were cached. Separating
the two is not bookkeeping pedantry — **the cached column is the only surviving record
of what the pilot actually cost**, and it is what exposes the retraction below.

| stage | ran here | CPU-h | resumed | CPU-h |
|---|---:|---:|---:|---:|
| `ANTISMASH` | 84 | **120.46** | — | — |
| `RECOVER_ORFS` | 21 | **47.93** | — | — |
| `GTDBTK_CLASSIFY` | — | — | 1 | **69.85** |
| `FETCH_RENAME_SCREEN` | — | — | 915 | 19.44 |
| `BUILD_PROTEIN_POOL` | — | — | 1 | 2.47 |
| `BIGSCAPE` | 1 | 1.61 | — | — |
| `GCF_BIOSYNTHETIC_TREE` | 1 | 1.58 | — | — |
| `BIGSCAPE_REFERENCES` | 1 | 1.01 | — | — |
| everything else (19 stages) | 12 | 1.23 | 19 | 0.93 |
| | **120** | **173.82** | **936** | **92.69** |

### Unit costs, for sizing anything else

| | |
|---|---:|
| fetch + pepM screen | **3.06 CPU-s per genome fetched** |
| `RECOVER_ORFS` | **41.2 CPU-s per genome retained** |
| `ANTISMASH` | **103.4 CPU-s per genome retained** |
| `GTDBTK_CLASSIFY` | **199 CPU-s per genome classified** |
| `BIGSCAPE` | 1.61 CPU-h at 1,345 BGCs |
| protein pool | 19,818,123 unique proteins from 20,310,261 CDS; 6.89 GB FASTA, 7.19 GB database |

## The pilot's cost table is withdrawn

The earlier version of this document reported **428.5 CPU-h** split as `RECOVER_ORFS`
327.6, `ANTISMASH` 67.7, `AGGREGATE_TAXONOMY` 11.0, `FETCH_RENAME_SCREEN` 8.6,
`GTDBTK_CLASSIFY` 6.0, other 7.6. **None of that should be cited.**

Nextflow writes a resumed task into the trace carrying its *original* execution metrics,
so the rerun's cached rows are the pilot's own measurements. They contradict the
published table directly:

| stage | published for the pilot | pilot's actual, from its cached row |
|---|---:|---:|
| `GTDBTK_CLASSIFY` | 6.0 | **69.85** — 8 h 43 m × 8 cpus |
| `FETCH_RENAME_SCREEN` | 8.6 | **19.44** |
| `AGGREGATE_TAXONOMY` | 11.0 | **it never ran** |

The third is self-refuting and should have been caught when written: `AGGREGATE_TAXONOMY`
being skipped *was failure 4*. A stage cannot be both the documented example of a process
that never executed and 2.6% of the bill.

**The pilot's own trace no longer exists.** `trace`, `report` and `timeline` use
`overwrite = true` under a per-taxon path, so rerunning the same taxon destroyed it —
behaviour this repo documents as intended, which is correct for the common case and
costly here. The cached rows survive only because they were copied forward into the new
trace.

What can still be said about the pilot: it completed 1,052 tasks with no failures and
produced no report. Its totals are unrecoverable.

**This also means the headline "428.5 → 266.5" is not a result.** The defensible
statement about the ORF-recovery fix is the controlled comparison in
`recover_orfs_recalibration/`: same genomes, same pool, same box, **126.7 → 13.4
CPU-s/genome isolated, 9.5×**, with output verified byte-identical. In production this
run measured 41.2 CPU-s/genome.

## Finding 1 — good annotation makes ORF recovery *more* expensive

`RECOVER_ORFS` is costly here *because* the input is well annotated, which is backwards
from intuition. The pool is built **from the retained genomes**, so annotation quality is
an input to its size rather than a discount on it:

| | Enterobacterales (GenBank) | this sweep (RefSeq reference) |
|---|---:|---:|
| genomes retained | 1,530 | 4,193 |
| pool | 447 MiB | **6.89 GiB** |
| proteins | 1,177,087 | **19,818,123** |
| **proteins per retained genome** | **769** | **4,726** |

2.74× the genomes, **16.8× the proteins**. A PGAP-annotated genome carries ~4,700
proteins and this figure is 4,726 — every genome donating its full complement.
Enterobacterales averages 769 because most GenBank submissions carry no CDS translations
at all: **83 of 130 genomes (64%)** in one recovery batch had none.

So Enterobacterales was fast *because its input was badly annotated*.

### Retracted measurements, kept on the record

Three in this area, each a different failure of method:

- The 2026-10-18 figure of "~7.5 min/genome, 10× slower on diverse input" — **withdrawn**:
  every diamond process was in D state on a failing USB-SATA bridge. Re-measured clean the
  ratio is 11.8×, close to the withdrawn number by a completely different mechanism. Right
  magnitude, wrong cause, invalid measurement.
- The projection replacing it (126 CPU-s/genome, 147 CPU-h) — **optimistic**, because it
  timed one genome in isolation on an idle box while production shares 16 cores.
- The projection for the batched version (~35 CPU-h) — **optimistic again, same reason**,
  against 47.93 measured. The isolated ratio held exactly; the contention factor was 3.07×
  where 2.23× was assumed.

Three times now, an isolated benchmark has under-predicted production. **Multiply an
idle-box measurement by ~3 before quoting it as a production figure on this machine**, or
measure under load.

## Finding 2 — 938 singletons is dereplication, not diversity

One genome per species means a family can only form where **different species** share a
cluster. Enterobacterales gave 72 families from 1,303 regions with 16 singletons; this
gives 1,022 from 1,345 with 938.

**Do not quote 1,022 families as a diversity result.** The comparable prior is the
640-genome Enterobacterales pilot, which gave 40 GCFs from 50 BGCs and was initially
misread the same way.

What is solid: **1,345 regions across 4,193 phylogenetically distinct species**, and a
largest family of **94 species sharing one cluster**.

Note the pilot reported 1,020 / 935 / 93 against this run's 1,022 / 938 / 94 on an
identical 1,345 regions. The difference is small and **unexplained**: `PYTHONHASHSEED` is
pinned in all four BiG-SCAPE processes and membership at the 0.30 cutoff has been stable
under that pin before. Given the pilot's other numbers did not survive scrutiny, treat
this run's as the measured ones rather than assuming the clustering moved.

See `hgt_candidates/` for what the family structure does support.

## Finding 3 — four failures, three of them silent

The pilot's most useful output. All four are fixed and **all four are confirmed fixed in
production** — verified by the stages being *present* in the rerun's trace, not by the
absence of errors.

| # | failure | how it presented | fixed in |
|---|---|---|---|
| 1 | `RECOVER_ORFS` killed at exactly 2 h | exit 143; `process_medium`'s limit, no retry label | labels override |
| 2 | GCF × species heatmap, 225,819 px | `ValueError` from `savefig` — the last call, after 430 CPU-h | [#13](https://github.com/alpole23/ClusterQuest/pull/13) |
| 3 | GCF tree, 82,620 px | same bug, sibling script, after the first was "fixed" | [#13](https://github.com/alpole23/ClusterQuest/pull/13) |
| 4 | `VISUALIZE_RESULTS` never ran | **`[SUCCESS] completed=7 failed=0`** | [#13](https://github.com/alpole23/ClusterQuest/pull/13) |

**Failure 4 is the one to remember.** `NCBI_FETCH_METADATA` emits `taxonomy_report.jsonl`
as `optional: true`, NCBI returns none for "Bacteria", and `EXTRACT_TAXONOMY` takes it as
a plain `path` input — so the empty channel skipped that process, then
`AGGREGATE_TAXONOMY`, then the report. A process that does not run is not a process that
failed, so Nextflow printed success.

The file is **never read**: `extract_taxonomy.py` takes `assembly_report`, `taxdump_dir`,
`output_file`. A file nothing opens gated the pipeline's deliverable.

**`failed=0` says nothing about what did not run.** Check the trace for the processes you
*expected*, not only for the ones that broke. The rerun was checked that way:
`EXTRACT_TAXONOMY`, `AGGREGATE_TAXONOMY`, `GCF_BIOSYNTHETIC_TREE` and `VISUALIZE_RESULTS`
all present.

Failure 5, nearly: GTDB-Tk refuses to start if any genome id contains `= ( ) [ ] ; ,` or a
space. Five of 4,193 carry a semicolon from NCBI strain fields — one bad name kills the
whole shard, 40 s in, after everything upstream.

## Going wider: full RefSeq Bacteria

Counts queried 2026-10-26:

| Bacteria at NCBI | genomes |
|---|---:|
| RefSeq, reference only | 22,893 ← **this sweep** |
| **RefSeq, all** | **544,787** |
| GenBank, all | 3,419,995 |

### Cost, from this run's unit costs

**Retention is the one unknown, and it dominates.** Full depth dilutes: Enterobacterales
went 8.8% on its reference set to **0.87%** on full RefSeq, a 10× drop, because the extra
genomes are redundant clinical strains. Bacteria should dilute less, since the carriers
are Pseudomonadota / Actinomycetota / Bacillota rather than the *E. coli* / *Salmonella* /
*Klebsiella* bulk. Bracket **3–6%**; the table below assumes **30,000 retained (5.5%)**.

| stage | scaling | CPU-h |
|---|---|---:|
| fetch + screen | linear in fetched | 463 |
| `RECOVER_ORFS` | **quadratic** — pool grows with retained | **~1,380** |
| `ANTISMASH` | linear in retained | 862 |
| `GTDBTK_CLASSIFY` | linear in BGC-positive (~9,000) | ~500 |
| `BIGSCAPE` | quadratic in BGCs (~9,600) | ~83 |
| **total** | | **~3,300 CPU-h** |

At 16 cores and realistic utilisation that is **~12 days of compute**, plus ~3.9 days of
download wall that partly overlaps: call it **2–3 weeks**. Disk ~1–1.5 TB against 12 TB
available.

**This supersedes an estimate of ~1,250 CPU-h made on 2026-10-25**, which was built from
the pilot's withdrawn stage costs — it under-counted antiSMASH and GTDB-Tk by using
figures that were never real.

### The memory blocker, found and removed

`BUILD_PROTEIN_POOL` kept the deduplicated pool resident in a Python dict keyed on the
protein **sequence** — ~550 bytes per unique protein, **measured 10.2 GB over 19.8 M**,
linear in the pool. At ~30,000 retained genomes that is ~68 GB, which does not fit a
56 GB box: the binding constraint on going wider was memory, not CPU.

Fixed in [#16](https://github.com/alpole23/ClusterQuest/pull/16) by remembering a sha256
digest and streaming each protein straight to the output, so the pool is never resident.
Measured on 400 genomes: **516 → 108 bytes per protein**, output byte-identical.

| retained genomes | unique proteins | before | after |
|---:|---:|---:|---:|
| 4,193 *(measured)* | 19.8 M | 10.2 GB | **~2.2 GB** |
| 30,000 | ~142 M | ~68 GB | **~14.5 GB** |

**Resource declarations now scale with the input rather than with this workstation**,
also in #16 — `BUILD_PROTEIN_POOL` from the genome count, `RECOVER_ORFS` from the
database size. A literal fitted to one box wastes a large node and silently blocks
scheduling on a small one, and this pipeline has to run on better hardware than the
machine it was developed on.

### What else to distrust at that scale

- **A 51 GB DIAMOND database.** `RECOVER_ORFS` was measured against 7.2 GB and peaks at
  7.26 GB of RAM, driven by query size. Neither the time nor the memory of a 7× larger
  database has been measured. The declaration now scales with the database (#16), but
  that rule is fitted to a single point: diamond chunks its database, so peak may be
  *bounded* rather than linear, in which case it over-requests. **This is the biggest
  remaining unknown in the table above**, and it sits on the stage that dominates.
- **BiG-SCAPE fits, which was not expected.** ~9,600 BGCs gives 20.7 GB on the corrected
  fit, against a ~15,000-BGC ceiling for this box — so it runs unpartitioned and the
  partitioning machinery stays insurance.
- **GenBank at 3.42 M is out.** The download alone is ~24 days of wall time.

### The run gates itself

Fetch and screen cost ~463 CPU-h and ~4 days, and report the **real** retained count,
from which every other figure here is recomputable. The decision point is four days in,
not three weeks in.

## Data

Trace at `results/pipeline_info/Bacteria/pipeline_trace.tsv` (rerun; the pilot's was
overwritten). Disk samples and logs under `/mnt/bigdisk/clusterquest/bench/`. Related:
`recover_orfs_recalibration/`, `hgt_candidates/`,
`enterobacterales_genbank/PREDICTIONS.md`.
