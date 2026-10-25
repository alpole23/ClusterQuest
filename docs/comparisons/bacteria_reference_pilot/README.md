# Whole-Bacteria reference pilot (2026-10-19 → 10-25)

One NCBI reference genome per species across all of Bacteria, to answer whether
a domain-wide phosphonate sweep is tractable on one workstation and what it
would find. It is. The run cost **428.5 CPU-h** over 1,052 tasks and returned
1,345 BGC regions.

It also broke four times, and the failures are the more useful half of this
document: three of them produced no error a reader would notice.

## Configuration

```json
{"taxon": "Bacteria", "assembly_source": "RefSeq",
 "assembly_reference_only": true, "outdir": "/mnt/bigdisk/clusterquest/results"}
```

`assembly_reference_only` **requires** `assembly_source = "RefSeq"` — reference
designation exists only on GCF\_ accessions, and the combination with GenBank
silently returns zero genomes. Measured: Enterobacterales GenBank 0, RefSeq 643,
Bacteria 22,871.

## Result

| | |
|---|---:|
| genomes fetched | 22,875 |
| retained by the pepM screen | **4,193 (18.3%)** |
| BGC regions | **1,345** |
| GCF families | 1,020 |
| — singletons | 935 |
| — multi-member | 85 (largest 93, 27 with ≥4) |
| GTDB-Tk classified | 1,263 genomes, 24 phyla |
| tasks / failures | 1,052 / **0** |
| total CPU | **428.5 CPU-h** |
| peak disk / memory | 568 GB / 35 GB (2 GB swap) |

Taxonomic spread of carriers: Pseudomonadota 582, Actinomycetota 271, Bacillota
230, Bacteroidota 43, Campylobacterota 30, across 24 phyla.

### Cost by stage

| stage | CPU-h | share |
|---|---:|---:|
| `RECOVER_ORFS` | **327.6** | **76.4%** |
| `ANTISMASH` | 67.7 | 15.8% |
| `AGGREGATE_TAXONOMY` | 11.0 | 2.6% |
| `FETCH_RENAME_SCREEN` | 8.6 | 2.0% |
| `GTDBTK_CLASSIFY` | 6.0 | 1.4% |
| everything else | 7.6 | 1.8% |

**The shape inverts again.** On RefSeq Enterobacterales the fetch was 40.4% and
recovery 20%; here recovery is 76% and the fetch 2%. Part of that is the fetch
being cached from an earlier session, but the dominant cause is real and is the
pilot's main operational finding.

## Finding 1 — good annotation makes ORF recovery *more* expensive

`RECOVER_ORFS` ran at **281 CPU-s/genome**, against 45.3 on GenBank
Enterobacterales. The input here is better annotated, which should give recovery
*less* to do.

The pool is built **from the retained genomes**, so annotation quality is an
input to its size rather than a discount on it:

| | Enterobacterales (GenBank) | this pilot (RefSeq reference) |
|---|---:|---:|
| genomes retained | 1,530 | 4,193 |
| pool | 447 MiB | **6.89 GiB** |
| proteins | 1,177,087 | **19,818,123** |
| **proteins per retained genome** | **769** | **4,726** |

2.74× the genomes, **16.8× the proteins**. A PGAP-annotated genome carries
~4,700 proteins and the pilot's figure is 4,726 — every genome donating its full
complement. Enterobacterales averages 769 because most GenBank submissions carry
no CDS translations at all: **83 of 130 genomes (64%)** in one recovery batch had
no annotation, logged as `no existing annotation, leaving it to antiSMASH`.

So Enterobacterales was fast *because its input was badly annotated*. Full
analysis and the fix in `docs/comparisons/recover_orfs_recalibration/`.

### A retraction and a correction

The 2026-10-18 figure of "~7.5 min/genome, 10× slower on diverse input" is
**withdrawn**: every diamond process was in D state on a failing USB-SATA bridge
at the time. Re-measured clean, the ratio is 11.8× — close to the withdrawn
number, by a completely different mechanism. Right magnitude, wrong cause,
invalid measurement.

Separately, the projection *replacing* it was also wrong, by 2.23×. It said 126
CPU-s/genome and 147 CPU-h; production gave **281 and 327.6**. The benchmark
timed one genome in isolation with four uncontended threads and a warm page
cache, while production runs 7–8 tasks sharing 16 cores at ~2.3 effective CPUs
each. The ratio the fix delivers held; its absolute did not.

## Finding 2 — 935 singletons is dereplication, not diversity

One genome per species means a family can only form where **different species**
share a cluster. Enterobacterales gave 72 families from 1,303 regions with 16
singletons; this gives 1,020 from 1,345 with 935.

**Do not quote 1,020 families as a diversity result.** The comparable prior is
the 640-genome Enterobacterales pilot, which gave 40 GCFs from 50 BGCs and was
initially misread the same way.

What is solid: **1,345 regions across 4,193 phylogenetically distinct species**,
and a largest family of **93 species sharing one cluster**.

This artefact surfaced three separate ways — as the singleton count, as a figure
that could not render, and as a 98%-empty matrix. See
`docs/comparisons/hgt_candidates/` for what the family structure does support.

## Finding 3 — four failures, three of them silent

Worth reading even if the numbers above are all you want.

| # | failure | how it presented |
|---|---|---|
| 1 | `RECOVER_ORFS` killed at exactly 2 h | exit 143; `process_medium`'s limit, no retry label |
| 2 | GCF × species heatmap, 225,819 px | `ValueError` from `savefig` — the last call, after 430 CPU-h |
| 3 | GCF tree, 82,620 px | same bug, sibling script, after the first was "fixed" |
| 4 | `VISUALIZE_RESULTS` never ran | **`[SUCCESS] completed=7 failed=0`** |

**Failure 4 is the one to remember.** `NCBI_FETCH_METADATA` emits
`taxonomy_report.jsonl` as `optional: true`, NCBI returns none for "Bacteria",
and `EXTRACT_TAXONOMY` takes it as a plain `path` input — so the empty channel
skipped that process, then `AGGREGATE_TAXONOMY`, then the report. A process that
does not run is not a process that failed, so Nextflow printed success.

The file is **never read**: `extract_taxonomy.py` takes `assembly_report`,
`taxdump_dir`, `output_file`. A file nothing opens gated the pipeline's
deliverable, and the Enterobacterales copy that let the branch run is 104 bytes.

**`failed=0` says nothing about what did not run.** Check the trace for the
processes you *expected*, not only for the ones that broke.

Failure 5, nearly: GTDB-Tk refuses to start if any genome id contains
`= ( ) [ ] ; ,` or a space. Five of 4,193 carry a semicolon from NCBI strain
fields — one bad name kills the whole shard, 40 s in, after everything upstream.

All five fixed in [#13](https://github.com/alpole23/ClusterQuest/pull/13).

## What this says about going wider

**Fetch and screen scale fine.** 1.36 CPU-s/genome here against 1.57 on
Enterobacterales — clade-insensitive, and 2.0% of the bill.

**Recovery is the constraint, and it is quadratic.** The pool grows with
retained genomes and every genome is searched against all of it. Building the
DIAMOND database once per run instead of once per genome removes one factor
([#12](https://github.com/alpole23/ClusterQuest/pull/12)) but not the quadratic:
blastx still searches each genome against a pool that grows with N.

**Dereplicating that pool does not fix it**, measured rather than assumed:
`diamond linclust` at 90% gives 14,496,200 clusters from 19,818,123 — 27% — and
it does **not saturate**: a 4× increase in genomes drops novel clusters per
genome only 23%, from 3,855 to 2,971. On taxonomically diverse input every new
species brings ~3,000 genuinely novel proteins. Scoping the query to pepM
neighbourhoods fails too, for a different reason: a 46× smaller query saves 11%,
because blastx against a 7.19 GB database is ~99 CPU-s of fixed load plus
~0.028 CPU-s per kb searched.

**The lever was the invocation count, and it has been measured** (2026-10-25).
diamond pays that fixed cost once per call and the module called it once per
genome. Batching every genome's gaps into one query takes this stage from
**126.7 to 13.4 CPU-s/genome** isolated — 9.5× — which projects the pilot's
327.6 CPU-h to **~35**, and `RECOVER_ORFS` from 76.4% of the bill to ~14%.

Verified byte-identical rather than assumed: 13/13 GFF3s against the per-genome
path, and 12/12 for the same genomes against a 200-genome batch. The memory
allocation moves 8 → 12 GB, because diamond peaks at 7.26 GB on a full batch and
neither way of bounding that by reshaping the query is safe — both change the
hit set. Detail in `recover_orfs_recalibration/`;
[#12](https://github.com/alpole23/ClusterQuest/pull/12).

Note this removes a constant factor, not the quadratic: blastx still searches
each genome against a pool that grows with N.

## Data

Trace, disk samples and logs under `/mnt/bigdisk/clusterquest/bench/`. Related:
`recover_orfs_recalibration/`, `hgt_candidates/`,
`enterobacterales_genbank/PREDICTIONS.md`.
