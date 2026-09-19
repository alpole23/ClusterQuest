# CLAUDE.md

Nextflow pipeline for analyzing phosphonate biosynthetic gene clusters (BGCs) in bacterial genomes using antiSMASH with optional BiG-SCAPE clustering and GTDB-Tk phylogenetic analysis.

## Environment Setup

```bash
conda activate nextflow    # Activate conda environment before running
```

## Quick Start

```bash
nextflow run main.nf --taxon "Pantoea ananatis"           # Full pipeline
nextflow run main.nf -resume                               # Resume previous run
nextflow run main.nf --workflow download --taxon "Streptomyces coelicolor"
nextflow run main.nf --clustering bigscape                 # With clustering
nextflow run main.nf -profile slurm                        # HPC execution
```

## Cross-Taxon Result Reuse

When analyzing a taxon that's a subset of a previously analyzed taxon, you can reuse existing antiSMASH and GTDB-Tk results to avoid redundant computation.

### Usage

```bash
# First run on broad taxon (e.g., family level)
nextflow run main.nf --taxon "Erwiniaceae"

# Later, run on subset taxon, reusing results from the broader run
nextflow run main.nf --taxon "Pantoea" --reuse_antismash_from "Erwiniaceae" --reuse_gtdbtk_from "Erwiniaceae"
```

**Direction matters:** reuse only works broad → narrow (superset → subset). Running Erwiniaceae after Pantoea can reuse antiSMASH results (`--reuse_antismash_from "Pantoea"`) since the check is per-genome and Pantoea genomes will get cache hits while other genera run fresh. But GTDB-Tk reuse will not work in this direction — it requires *all* genomes in the current run to exist in the prior results, which fails when the new taxon is broader.

### antiSMASH Reuse

1. For each genome in the current run, the pipeline checks if results exist in the reuse directory
2. Results are reused if:
   - The antiSMASH version matches (major.minor, e.g., 7.1.x matches 7.1.y)
   - The parameter configuration matches (tracked via hash)
3. Genomes without existing results are processed normally
4. Each antiSMASH result includes a `.antismash_meta` file that stores version and params_hash

### GTDB-Tk Reuse

1. The pipeline checks if ALL genomes in the current run exist in the reuse results
2. If yes: The summary TSV is filtered and the phylogenetic tree is pruned to only include current genomes
3. If no: GTDB-Tk runs fresh on all current genomes (all-or-nothing approach)

This is more efficient than re-running GTDB-Tk, especially since the classify step (pplacer) is memory-intensive.

### Clustering

BiG-SCAPE always runs fresh for the current genome set, as clustering depends on the complete set of BGCs being analyzed together.

### When to Use

- Running on a genus after analyzing the family (e.g., Pantoea after Erwiniaceae)
- Re-running with different clustering parameters (antiSMASH results unchanged)
- Adding new genomes to a previous analysis

## Configuration

Parameters in `nextflow.config` are organized by subworkflow to make it easy to find relevant settings.

### Global

| Parameter | Default | Description |
|-----------|---------|-------------|
| `workflow` | "full" | Pipeline mode: `download`, `bgc_analysis`, or `full` |
| `outdir` | "results" | Output directory for all results |
| `task_batch_size` | 100 | Genomes per task for the short per-genome steps (rename, GenBank→FASTA, antiSMASH result copy) |

### DOWNLOAD_GENOMES

Downloads and prepares bacterial genomes from NCBI.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `taxon` | "Erwiniaceae" | NCBI taxon (species, genus, family, order, etc.) |

### ANTISMASH_ANALYSIS

BGC detection using antiSMASH.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_genomes` | null | Path to pre-downloaded genomes (for `bgc_analysis` workflow) |
| `reuse_antismash_from` | null | Taxon name to reuse antiSMASH results from |
| `antismash_minimal` | false | Minimal mode (faster, skips domain analysis) |
| `antismash_cb_general` | false | ClusterBlast: Compare vs antiSMASH DB |
| `antismash_cc_mibig` | false | ClusterCompare: Advanced MIBiG scoring |
| `antismash_smcog_trees` | false | Phylogenetic trees for BGC genes |
| `antismash_phosphonate_neighbourhood` | 10 | kb of flank kept around the rule core (antiSMASH's own value is 5); `null` restores it |

**Note:** Detection is hardcoded to phosphonate rule only (`--hmmdetection-limit-to-rule-names phosphonate`). `--cb-knownclusters`, `--clusterhmmer`, and `--tigrfam` are always enabled. `--no-zip-output` is also always passed: the `{genome}.zip` antiSMASH writes by default is an archive of its own output directory (~6 MB/genome) that nothing downstream reads. The whole-genome summary GenBank is off by default too and gated behind `--antismash_summary_gbk` (~11 MB/genome, also unread by any step) — worth enabling on small sets, not on a genus. Neither flag is in `Utils.antismashParamsHash`: they change packaging, not results, so toggling them does not invalidate `--reuse_antismash_from`.

### Region Analysis

BGC counting, tabulation, and statistics.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_analysis` | true | Enable region analysis and visualization |
| `count_per_contig` | false | Count per contig (true) or per genome (false) |
| `split_hybrids` | false | Split hybrid types (T1PKS-NRPS → T1PKS + NRPS) |

### CLUSTERING

Gene Cluster Family (GCF) clustering.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `clustering` | "bigscape" | `none` or `bigscape` |

**BiG-SCAPE Options** (when `clustering = "bigscape"`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bigscape_cutoffs` | "0.30" | GCF distance threshold(s), comma-separated |
| `bigscape_alignment_mode` | "auto" | `auto`, `global`, or `glocal` |
| `bigscape_mibig_version` | "" | MIBiG version (e.g., "3.1") or "" to exclude |
| `bigscape_classify` | "category" | `""`, `category`, `class`, or `legacy` |
| `bigscape_include_singletons` | true | Include unclustered BGCs |
| `bigscape_mix` | false | Allow mixing BGC classes in same GCF |
| `bigscape_reference_dir` | `assets/phosphonate_reference_bgcs` | Characterised clusters measured by `BIGSCAPE_REFERENCES`; `""` disables |

| Parameter | Default | Description |
|-----------|---------|-------------|

### PHYLOGENY

Phylogenetic placement using GTDB-Tk. ⚠️ Requires ~140 GB disk and ~56-64 GB RAM.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_gtdbtk` | true | Enable phylogenetic analysis |
| `reuse_gtdbtk_from` | null | Taxon name to reuse GTDB-Tk results from |
| `gtdbtk_bgc_genomes_only` | true | Only analyze genomes with detected BGCs |
| `gtdbtk_cpus` | 8 | CPUs for identify/align steps |
| `gtdbtk_pplacer_cpus` | 1 | CPUs for pplacer (keep low, memory-bound) |
| `gtdbtk_min_perc_aa` | 10 | Minimum % amino acids in MSA |
| `gtdbtk_outgroup` | null | Outgroup pattern for tree rooting (e.g., "g__Escherichia") |

## Output Structure

```
results/
├── databases/                    # Cached databases (preserve these)
├── ncbi_genomes/${taxon}/
│   ├── ncbi_dataset/            # Raw NCBI download
│   ├── renamed_genomes/         # Standardized genome files
│   └── name_map.json            # Assembly ID to genome name mapping
├── antismash_results/${taxon}/  # Per-genome antiSMASH results
├── bigscape_results/${taxon}/   # BiG-SCAPE clustering output
│   ├── ${taxon}.db              # SQLite database with clustering results
│   ├── gcf_representatives.json # GCF data with KCB hits and gene diagrams
│   ├── reference_distances.tsv  # every BGC's distance to each characterised cluster
│   └── reference_summary.json   # per-GCF nearest reference, per-reference nearest GCF
├── gtdbtk_results/${taxon}/     # GTDB-Tk phylogenetic placement
├── pipeline_info/
│   ├── pipeline_trace.tsv       # Per-task timing and resource data
│   ├── pipeline_report.html     # Nextflow execution report
│   ├── pipeline_timeline.html   # Visual timeline
│   └── software_versions.json   # Tool versions
└── main_analysis_results/${taxon}/
    ├── region_counts.tsv        # BGC counts per genome
    ├── region_tabulation.tsv    # Detailed BGC information
    ├── taxonomy_map.json        # Genome taxonomy mapping
    ├── bgc_report.html          # Interactive HTML report
    ├── genomes/                 # Per-genome HTML pages
    ├── kcb_identification_chart.png
    ├── rarefaction_curve.png
    ├── gcf_heatmap/             # GCF biosynthetic tree and heatmap outputs
    │   ├── gcf_biosynthetic_tree.png
    │   ├── gcf_biosynthetic_tree.svg
    │   ├── gcf_species_heatmap.png
    │   ├── gcf_species_heatmap.svg
    │   └── phosphonate_metadata.json
    └── coupling_enzyme_trees/   # Coupling enzyme phylogenetic trees
        ├── coupling_trees_manifest.json
        ├── tree_A/              # Combined pepM tree (all BGCs + references)
        │   ├── pepm_tree.nwk
        │   ├── itol_coupling_class.txt
        │   ├── itol_gcf.txt
        │   ├── itol_source.txt
        │   ├── itol_organism.txt
        │   └── missing_pepm.txt    # BGCs with no extractable pepM (if any)
        └── tree_B/              # Per-class coupling enzyme trees
            ├── Synthase/
            │   ├── Synthase_tree.nwk
            │   └── itol_*.txt
            ├── Decarboxylase/
            │   ├── Decarboxylase_tree.nwk
            │   └── itol_*.txt
            ├── Reductase/
            │   └── ...
            └── Transaminase/
                └── ...
```

## SLURM HPC Execution

### Configuration

The pipeline uses process labels from `conf/labels.config` with SLURM-specific overrides in `nextflow.config`. The SLURM profile (`-profile slurm`) provides:

- Queue management: 200 concurrent jobs, 20 submissions/min
- High-memory queue for `process_high_memory` label (GTDB-Tk)
- BiG-SCAPE memory override: 64 GB (scales O(n²))
- Extended time limits for network-bound processes

Default resources are defined by labels and can be overridden in the SLURM profile.

**Submission rate is the scaling constraint, not compute.** At 20 submissions/min, job
*count* matters more than job duration for large taxa — which is why the short per-genome
steps are batched (see Task Batching) and why `ANTISMASH` is capped at 2 CPUs so more
genomes run concurrently per allocation. Before adding a new per-genome process, check
whether it should be batched instead.

### Running on SLURM

```bash
# Direct execution
nextflow run main.nf -profile slurm --taxon "Pantoea"

# Via sbatch (recommended for long runs)
sbatch submit_slurm.sh
```

## Benchmarking

### Purpose

Benchmark the pipeline to estimate runtime for large-scale analyses (e.g., 3 million genomes).

### Recommended Benchmark Size

| Sample Size | Use Case |
|-------------|----------|
| 100-200 | Good statistical power, captures variance |
| 300-500 | Excellent confidence intervals, tests SLURM scaling |
| 1000+ | Very accurate but excessive for benchmarking |

### Running a Benchmark

1. **Pre-download databases** (one-time cost, not included in benchmark):
   ```bash
   # Databases are cached in results/databases/
   # antiSMASH (~50GB), GTDB-Tk (~140GB), Pfam, TaxonKit
   ```

2. **Run benchmark with trace enabled** (already configured in nextflow.config):
   ```bash
   nextflow run main.nf -profile slurm --taxon "Pantoea"
   ```

3. **Analyze trace file**:
   ```bash
   # Trace file location: results/pipeline_info/pipeline_trace.tsv
   # Contains: task_id, name, status, realtime, cpus, memory, peak_rss, peak_vmem
   ```

### Extrapolation Formula

```
Total time = (mean_time_per_genome × total_genomes) / max_parallel_jobs

# With 95% confidence interval:
CI = mean ± (1.96 × std_dev / sqrt(n))
```

### Key Metrics

- **Wall time per genome**: Primary scaling factor (antiSMASH is bottleneck)
- **CPU hours per genome**: For HPC allocation requests
- **Peak memory per genome**: For SLURM memory allocation
- **Success/failure rate**: For planning retries

### Example: Pantoea Benchmark

- **Taxon**: Pantoea (~1500 genomes)
- **Configuration**: Full analysis (antiSMASH + BiG-SCAPE + GTDB-Tk)
- **Expected output**: Per-genome timing data for 3M genome extrapolation

### Measured: Pantoea and Erwiniaceae (2026-08-29)

Two full local runs, no `--reuse_antismash_from`, on one 56 GB box. Reference figures for
sizing resources — re-measure after any change to the antiSMASH parameter set. Full
write-up with the cost model: `docs/benchmark_2026-08-29.html`.

| | *Pantoea* (genus) | *Erwiniaceae* (family) | ratio |
|---|---|---|---|
| genomes | 1,736 | 2,758 | 1.59x |
| tasks / failures | 1,771 / 0 | 2,803 / 0 | — |
| CPU-hours | 60.43 | 75.55 | **1.25x** |
| wall clock | 4 h 50 m | 5 h 56 m | 1.23x |
| BGCs / GCFs | 320 / 13 | 333 / 19 | 1.04x / 1.46x |
| BGC-positive genomes | 285 (16.4%) | 298 (10.8%) | — |
| I/O (rchar+wchar) | 7,254 GB | 8,877 GB | 1.22x |

| Process | CPU-h (Pantoea → Erwiniaceae) | scaled | share |
|---------|-------------------------------|--------|-------|
| `ANTISMASH` | 47.28 → 61.58 | 1.30x | 81.5% |
| `GTDBTK_CLASSIFY` | 12.08 → 12.42 | 1.03x | 16.4% |
| `NCBI_DATASETS_DOWNLOAD` | 0.51 → 0.86 | 1.68x | 1.1% |
| `BIGSCAPE` | 0.14 → 0.16 | 1.14x | 0.2% |

**antiSMASH is O(n) and dominant today.** BiG-SCAPE was measured separately (below).

**GTDB-Tk's scaling is linear, measured 2026-09-18.** An earlier version of this file
claimed it was fixed-dominated ("59% more genomes cost 3% more, therefore O(1)"), which
compared *total* genomes (1,736 → 2,758) while `gtdbtk_bgc_genomes_only` means GTDB-Tk
only ever saw the BGC-positive ones: 285 → 298, **4.6% apart**. Two points that close
cannot separate fixed cost from linear cost, and the two models differed by 94x at a
million genomes.

Measured directly instead, single shard, module flags, 8 vs 64 genomes:

| N | wall | CPU | peak RSS |
|---:|---:|---:|---:|
| 8 | 31.5 min | 28.2 min | 55.7 GB |
| 64 | 53.7 min | 63.5 min | 56.3 GB |

**slope 0.630 CPU-min/genome, intercept 23.2 CPU-min.** The fixed cost is real but small;
per-genome cost dominates above ~40 genomes. At the ~121,000 BGC-positive genomes a
million-genome run yields, that is **~1,271 CPU-h**, against the 48 CPU-h the cost model
had assumed — so the million-genome total moves from 22,868 to roughly 24,100 CPU-h, ~5%.

Peak memory is set by the reference data, not by genome count (55.7 vs 56.3 GB), and sits
close enough to this box's 56 GB that the figures may be memory-bound.

`GTDBTK_CLASSIFY` is sharded (`gtdbtk_shard_size`, default 5,000). Each shard re-pays only
the 23.2 CPU-min intercept, so 20 shards cost ~7.7 CPU-h on ~1,271 — **0.6%**, which makes
the sharding cheap insurance rather than a gamble. Data:
`docs/comparisons/gtdbtk_scaling/`.

### Measured: BiG-SCAPE scaling (2026-08-30)

`scripts/bench_bigscape_scaling.py`, ten sizes from 333 to 10,000 BGCs, built by
replicating real Erwiniaceae region GBKs with unique identities (BiG-SCAPE dedupes on the
sha256 of raw file bytes).

| BGCs | comparisons | wall | CPU-s | us/pair | peak RSS | db | B/pair |
|-----:|------------:|-----:|------:|--------:|---------:|---:|-------:|
| 333 | 55,278 | 45.1 s | 277 | 816 | 1.62 GB | 15 MB | 267 |
| 1,000 | 499,500 | 148.1 s | 911 | 297 | 2.01 GB | 75 MB | 151 |
| 2,000 | 1,999,000 | 350.4 s | 2,217 | 175 | 2.54 GB | 242 MB | 121 |
| 4,000 | 7,998,000 | 991.3 s | 6,077 | 124 | 3.40 GB | 856 MB | 107 |
| 6,000 | 17,997,000 | 32.4 m | 11,148 | 108 | 5.68 GB | 1.86 GB | 103 |
| 8,000 | 31,996,000 | 55.4 m | 18,729 | 104 | 9.08 GB | 3.25 GB | 102 |
| 10,000 | 49,995,000 | 85.6 m | 28,535 | 103 | 14.21 GB | 5.06 GB | 101 |

**Comparisons are all-pairs at every size** — the `distance` table holds exactly n(n-1)/2
rows, n^2.000, no pruning. That is the check the benchmark rests on; if it ever stops
holding, the curve is measuring something else.

**CPU: ~1,090 CPU-h at 121,000 BGCs, about $54.** Cost per comparison falls 816 -> 103 us
and has flattened; the local exponent climbs 1.05 -> 1.89 and is converging on 2. Fit
`cpu_s = 1643 + 2.68e-4*n^2` on the top four sizes (worst residual 2.4%). The harness also
prints a power-law fit (n^1.241, R2=0.99) — **do not extrapolate that one**, it understates
by 10x. An earlier 21,000 CPU-h estimate overstated by 20x for the mirror reason: it scaled
the whole 333-BGC runtime quadratically, fixed Pfam and database cost included.

**Memory is the wall, and it is invisible below 4,000 BGCs.** Peak RSS looks flat across
the first seven sizes (1.62 -> 3.40 GB) and an earlier draft concluded memory was not a
constraint. It was measured over a range where the growth had not started. The local
exponent runs 0.11 -> 0.29 -> 0.59 -> 1.27 -> 1.63 -> **2.01**: above ~4,000 BGCs it is
exactly quadratic. Fit `GB = 1.14 + 1.29e-7*n^2` gives **~1.9 TB at 121,000 BGCs**, 4.3 TB
at *Pantoea*-level BGC yield.

| RAM budget | largest single BiG-SCAPE job |
|-----------:|-----------------------------:|
| 48 GB (`process_high` label) | 19,100 BGCs |
| 64 GB (SLURM `withName: BIGSCAPE`) | 22,100 BGCs |
| 128 GB | 31,400 BGCs |
| 1 TB | 89,100 BGCs |

So a million genomes cannot be clustered in one pass. Partition by BGC class, or
pre-cluster, before that point. **The 64 GB SLURM override caps the pipeline at ~22,000
BGCs** — roughly 180,000 genomes at Erwiniaceae-level prevalence, which is where this will
first bite.

**Database: ~0.74 TB at 121,000 BGCs.** Bytes per comparison converged to 101, so 7.3
billion rows in one SQLite file (16.9 billion / 1.71 TB at the higher yield). Real, but a
smaller problem than the RAM.

One trap: BiG-SCAPE shells out to `fasttree` for its GCF trees. Nextflow supplies it via
the activated conda env; invoking the binary by path does not, and the run dies *after*
computing every distance. The harness prepends the executable's own bin to `PATH`.

**antiSMASH cost tracks base pairs, not genome count.** Per-genome cost *fell* 98.0 →
80.4 CPU-s when the taxon widened, because Erwiniaceae drags in 213 *Buchnera*
endosymbionts at ~0.6 Mb. Smallest to largest bin is a 6.2x runtime spread (6.9 s vs
42.6 s). Do not carry either figure to an arbitrary genome set without checking size
distribution.

**GTDB-Tk does not need 104 GB.** Nextflow's `peak_rss` reads ~104 GB and is an artefact:
it sums RSS across pplacer's forked children, which all map the same reference database,
and it did not move between the two runs (104.3 → 104.0 GB) despite 59% more queries.
A 30-second sampler over the Erwiniaceae run measured the truth — largest single process
**pplacer at 47.0 GB**, peak system memory in use **19.2 GB**, minimum available 37.6 GB,
peak swap 2.0 GB. Most of pplacer's 47 GB is file-backed `mmap`, reclaimable under cgroup
pressure. The SLURM profile still requests 128 GB as deliberate insurance (an OOM kill
costs a 90-minute task and its queue slot; the over-request costs ~$0.05) — but note it
tips the billed dimension from cores to memory at 8 cores. Verify with `sacct -o MaxRSS`
before paying that premium across sharded jobs.

**Storage, complete accounting (Erwiniaceae, 2,758 genomes).** Earlier figures counted
only `antismash_results/` and understated peak by 6.7x. `publishDir` used `mode: 'copy'`,
so `work/` held a second copy of everything.

| | GB | MB/genome |
|---|---:|---:|
| `results/ncbi_dataset` (raw download) | 24.0 | 8.91 |
| `results/renamed_genomes` | 24.0 | 8.91 |
| `results/antismash_results` | 21.0 | 7.79 |
| `results/` other | 0.4 | 0.13 |
| `work/` (download, rename, antiSMASH, GTDB-Tk) | 70.7 | 26.25 |
| **peak during a run** | **140.0** | **51.99** |

Every genome's sequence was stored **four times** — raw in `work/`, raw published, renamed
in `work/`, renamed published — 35.3 MB/genome, 68% of all storage. Two changes on
2026-08-30 address it:

- **`params.publish_mode = 'link'`** (was a hardcoded `mode: 'copy'` in all 17 modules).
  Hard links mean a published file and its `work/` counterpart share an inode, so
  `results/` costs no extra disk. Verified: the published copy survives
  `nextflow clean -f`, the inode's link count simply drops to 1. Needs `outdir` and
  `workDir` on one filesystem — set `'copy'` if they are not, or if anything edits
  published files in place, since a write through a hard link also rewrites the cached
  task output and corrupts `-resume`.
- **`NCBI_DATASETS_DOWNLOAD` no longer publishes the genomes**, only the metadata that
  `main.nf`'s `bgc_analysis` entry reads. The `*.gbff` payload is republished by
  `RENAME_GENOMES` anyway.

Together: ~52 -> ~16 MB/genome peak, and the 1M-genome projection goes 49.6 TB -> ~16 TB.

Two Nextflow traps found doing this, both verified against the real output declarations:
`publishDir`'s `pattern:` publishes **nothing at all** when any output is declared with a
`**` glob, and a directory output (`path "ncbi_dataset/"`) is published as a single item
that `saveAs` cannot filter inside — so the fix needed `saveAs` returning null *and*
removal of the directory output, which was emitted but never consumed.

**antiSMASH output after `--no-zip-output` / `--no-summary-gbk`:** 24.63 -> 8.55 MB per
genome, **-65.3%**, measured like-for-like on 40 genomes present in both runs. I/O did not
improve — ~3.2 GB per genome of small-file read+write is the shared-filesystem risk.

**Projected to 1M genomes** at UIUC internal rates (`max(cores x $1.19, GB x $0.08)`
per day, storage $8.75/TB/month): antiSMASH ~25,000 CPU-h / $1,240; GTDB-Tk in 50k
batches 563 CPU-h / $30; BiG-SCAPE ~1,090 CPU-h / $54 (measured, above); storage ~16 TB
with hard links and the raw-download publish dropped, plus 0.74 TB of BiG-SCAPE
database / $145 per month. Compute is not the constraint
at ~$1,320 — BiG-SCAPE's ~1.9 TB memory requirement is, and it forces partitioning.

## Module & Script Organization

```
main.nf                 # Parameter validation + entry point (~140 lines)
nextflow.config         # Parameters and SLURM profile

tests/                  # Test suite — run with: bash tests/run_tests.sh
├── run_tests.sh        # Entry point; runs everything in a scratch dir
├── test_utils.nf       # Utils.groovy helpers (optArg, isValidInput, sanitizeTaxon)
├── test_batching.nf    # Batched per-genome processes end to end
├── make_fixtures.py    # Synthetic genomes (incl. one corrupt) + fake reuse results
└── check_undefined.py  # Static scan for calls to undefined/unimported names

subworkflows/           # Workflow composition (one file per subworkflow)
├── helpers.nf          # placeholder(), clusteringEnabled(), batchSize() — included like processes
├── download_genomes.nf # DOWNLOAD_GENOMES
├── antismash_analysis.nf # ANTISMASH_ANALYSIS (with reuse)
├── clustering.nf       # CLUSTERING
├── phylogeny.nf        # PHYLOGENY (with reuse)
└── bgc_analysis.nf     # BGC_ANALYSIS (composes the four above)

assets/
└── reference_sequences/
    ├── reference_pepM.faa             # 7 PEP mutase references (Tree A anchors)
    └── reference_coupling_enzymes.faa # 7 coupling enzyme references (Tree B anchors)

conf/
├── conda.config        # Centralized conda environments by process
└── labels.config       # Process labels (resource allocations, error handling)

lib/
└── Utils.groovy        # Shared Groovy utilities (sanitizeTaxon, antismashParamsHash, buildReusePath)

modules/
├── databases/          # Database download processes (antiSMASH, GTDB-Tk, Pfam, etc.)
├── genome/             # Genome processing (NCBI download, rename, GenBank→FASTA)
├── analysis/           # BGC analysis (antiSMASH, counting, tabulation, reuse)
├── clustering/         # BiG-SCAPE clustering and stats extraction
├── phylogeny/          # GTDB-Tk classification (with reuse support)
├── visualization/      # HTML report generation
└── utilities/          # Version collection

scripts/
├── utils/              # Shared Python utilities
│   ├── constants.py      # BGC_COLORS (comprehensive ~110-entry palette), GENE_COLORS, KCB_THRESHOLDS,
│   │                     # COUPLING_COLORS, COUPLING_ORDER, LEGACY_CLASS_NAMES, GCF_PALETTE,
│   │                     # DOMAIN_NAMES + domain_name(accession),
│   │                     # load_coupling_classes(path, region_only=False)
│   ├── tree_building.py  # build_nj_tree(labels, dm_rows) — shared NJ tree construction
│   ├── tree_layout.py    # Cladogram layout/drawing (linear + circular) for Bio.Phylo trees
│   ├── bigscape_db.py    # BiG-SCAPE SQLite queries (BGC records, families, domains)
│   ├── itol.py           # iTOL dataset writers (colorstrip, binary, simplebar, text)
│   ├── bgc_labels.py     # label_from_path, make_labels_unique
│   ├── colors.py         # family_color, genus_color
│   ├── gene_diagram.py   # Gene arrow SVG generation
│   ├── parsers.py        # Duration, memory, timestamp parsing; format_bytes(v, precision=1)
│   ├── trace.py          # Nextflow trace aggregation + resource/Gantt HTML (single implementation;
│   │                     # re-exported by viz/ for convenience)
│   ├── plotting.py       # Deterministic matplotlib output: pins svg.hashsalt + Agg backend,
│   │                     # exports SVG_METADATA for timestamp-free savefig
│   └── antismash_parser.py  # antiSMASH JSON parsing + BGC label/region helpers
├── viz/                # Visualization modules — every report section lives here;
│   │                     # visualize_results.py only orchestrates them
│   ├── charts.py         # KCB pie charts, BGC color utilities
│   ├── tree_viz.py       # Prunes the GTDB-Tk tree to the analysed genomes (Bio.Phylo)
│   ├── tables.py         # Genome tables, summary statistics, distribution table
│   ├── clustering.py     # BiG-SCAPE stats + GCF representative HTML
│   ├── taxonomy.py       # Interactive taxonomy tree
│   ├── distribution.py   # GCF × genus heatmap, genus-specific / widespread tables
│   ├── genome_pages.py   # Per-genome HTML pages
│   ├── rarefaction.py    # GCF rarefaction curve
│   ├── report_assets.py  # REPORT_CSS / REPORT_JS for bgc_report.html
│   └── report_sections.py # _build_* section builders + build_coupling_table_rows
├── taxonomy/           # Taxonomy processing scripts
├── genome/             # Genome processing scripts
├── clustering/         # Clustering statistics and GCF representative extraction
├── analysis/           # BGC counting and tabulation
├── phylogeny/          # GTDB-Tk result filtering
├── bgc_coupling_tree.py  # Coupling enzyme NJ trees (Tree A: pepM, Tree B: per-class)
└── visualize_results.py  # Main visualization entry point
```

### Workflow Structure

The pipeline uses DSL2 subworkflows for modularity. Parameters in `nextflow.config` are organized to mirror this structure.

```
workflow (entry point)
│
├── DOWNLOAD_GENOMES          # Download and prepare genomes from NCBI
│   ├── NCBI_DATASETS_DOWNLOAD
│   ├── CREATE_NAME_MAP
│   ├── RENAME_GENOMES
│   └── EXTRACT_TAXONOMY
│
└── BGC_ANALYSIS              # Main analysis pipeline
    │
    ├── ANTISMASH_ANALYSIS    # BGC detection (with reuse support)
    │   ├── CHECK_ANTISMASH_REUSE
    │   ├── ANTISMASH
    │   └── COPY_ANTISMASH_RESULT
    │
    ├── Region Analysis       # BGC statistics
    │   ├── COUNT_REGIONS
    │   ├── TABULATE_REGIONS
    │   └── AGGREGATE_TAXONOMY
    │
    ├── CLUSTERING            # GCF clustering
    │   ├── BIGSCAPE
    │   ├── EXTRACT_CLUSTERING_STATS
    │   └── EXTRACT_GCF_REPRESENTATIVES
    │
    ├── PHYLOGENY             # GTDB-Tk (with reuse support)
    │   ├── CHECK_GTDBTK_REUSE
    │   ├── GTDBTK_CLASSIFY
    │   └── FILTER_GTDBTK_RESULTS
    │
    ├── GCF_BIOSYNTHETIC_TREE # GCF biosynthetic NJ tree (when bigscape enabled, runs before visualization)
    │                         # Also outputs phosphonate_itol_coupling.txt (coupling annotation)
    └── VISUALIZE_RESULTS     # HTML report generation (receives GCF tree PNG + coupling annotation)
```

**Running a single stage:**

```bash
# Run only download
nextflow run main.nf --workflow download --taxon "Pantoea"

# Run full pipeline (default)
nextflow run main.nf --taxon "Pantoea"
```

Note: Nextflow 26's strict parser rejects `-entry <WORKFLOW>`; select the stage with
`--workflow` instead.

### Process Labels

Processes use labels for resource allocation and error handling:

| Label | CPUs | Memory | Description |
|-------|------|--------|-------------|
| `process_local` | 1 | 1 GB | Runs on head node |
| `process_low` | 1 | 2 GB | Light scripts |
| `process_medium` | 4 | 8 GB | Visualization, GCF trees |
| `process_high` | 8 | 32 GB | BiG-SCAPE |
| `process_high_memory` | 8 | 128 GB | GTDB-Tk pplacer |

Per-process override in `conf/labels.config`:

| Process | CPUs | Memory | Rationale |
|---------|------|--------|-----------|
| `ANTISMASH` | 2 | 6 GB × attempt | Measured over 950 genomes: %cpu p99 = 145% (never 2 full cores), peak RSS p99 = 4.7 GB. Escalates once on an OOM kill, then skips the genome. |

Error handling labels:
- `tolerant`: Individual failures don't stop pipeline (per-genome processes)
- `retry_on_error`: Retry on transient errors (network downloads)

### Cache Invalidation: Nextflow Hashes Source, Not Rendered Script

**Nextflow hashes the *unevaluated* script block plus the declared input values.
It never hashes the rendered command.** Interpolating a value into the script
body therefore does not invalidate anything.

Measured on Nextflow 26.04 with a two-process probe:

| Where the changed value lives | Task re-runs? |
|---|---|
| Interpolated into the `script:` body | **No** — `cached=1` |
| Declared as a `val` input | **Yes** — `completed=1` |

This matters because scripts are invoked as `python ${projectDir}/scripts/foo.py`
— an interpolated path, not a declared input — so Nextflow cannot see the Python
file at all. `Utils.scriptsHash` exists to close that gap, but it was embedded as
a `# scripts-version:` **comment inside the script block**, which lands on the
wrong side of the table above. It invalidated nothing.

The consequence is silent: `-resume` reuses output produced by Python code that
has since changed, with no error and no warning. Found when a partitioned run
failed with `KeyError: 'genome'` because `PARTITION_BGCS` served a manifest built
by the previous version of `partition_bgcs.py` — the digest had changed from
`8cc0a42b6b97` to `cc4bf61ea8c5` and the task was reused regardless.

**The fix is to pass the digest as a `val` input**, as `PARTITION_BGCS` now does:

```groovy
// module
input:
val scripts_version

// call site
PARTITION_BGCS(taxon, antismash_results, pfam_db_ch,
               Utils.scriptsHash(projectDir,
                   ['clustering/partition_bgcs.py', 'utils']))
```

`ANTISMASH` was never affected: it already passes `antismash_version` and
`antismash_params_hash` as `val` inputs, which is the pattern that works.

**All 17 processes now take the digest as a `val` input.** The last two —
`create_name_map` and `rename_genomes_parallel` — were deferred because they sit
upstream of antiSMASH and converting them invalidates every antiSMASH task
(~2,807, ~6 h). They were converted once that cache had been cleared for other
reasons, when the re-run was already unavoidable and the conversion therefore
free.

`tests/check_script_deps.py` **rejects** the `# scripts-version:` comment spelling
outright, so it cannot reappear silently. It reads the digest only from the call
site.

`CHECK_GTDBTK_REUSE` needs no digest at all: it declares `cache false`. Its
sibling `FILTER_GTDBTK_RESULTS` in the same file is the process that runs Python,
and that is where the digest belongs — a reminder that one `.nf` file can hold
more than one process.

## Reference Database Versions

**Databases are pinned, and the pins are recorded in the output.** Every database
downloads through `storeDir`, which skips the process whenever its output path already
exists. An unpinned URL therefore tracks nothing: `releases/latest/` and
`current_release/` resolve exactly once, on the first run ever, and are never
re-checked. This pipeline classified against a **2026-01-23 NCBI taxdump for seven
months** that way, with no record anywhere of which version produced any result.

| Parameter | Default | Notes |
|-----------|---------|-------|
| `gtdb_release` | `226` | GTDB-Tk enforces `MIN_REF_DATA_VERSION`; 2.6.1 wants r226. **Bump the `gtdbtk` conda pin before bumping this.** |
| `pfam_release` | `38.2` | Current as of 2026-08 |
| `taxdump_date` | `2026-08-01` | NCBI *monthly archive* (`taxdump_archive/taxdmp_<date>.zip`), not the live `taxdump.tar.gz`, which is rewritten daily and would pin itself to whatever day you first ran |
| `check_db_updates` | `true` | Warn when upstream moves past a pin. One HTTP request per database at startup; set `false` for offline runs |

antiSMASH is not listed because its databases ship keyed to the tool release, so the
conda pin on `antismash` already pins them.

At startup the run prints what it is classifying against:

```
Reference databases:
  GTDB     226  <- upstream now 232; edit params.gtdb_release to upgrade
  Pfam     38.2  (current)
  taxdump  2026-08-01  (current)
```

`lib/DbVersions.groovy` does the check. It **never fails a run and never changes what is
downloaded** — an offline machine still runs, and every failure path returns null. The
pins are also written into `software_versions.json` as `db_gtdb_release`,
`db_pfam_release` and `db_taxdump_date`, so a published result can say what it was
produced against.

**Do not upgrade casually.** Two runs on different GTDB releases are not directly
comparable — genera get reclassified between releases, which moves the taxonomy tree and
the per-clade BGC prevalence the whole analysis rests on. Changing a pin changes the
`storeDir` output path, which is what triggers the fresh download; for GTDB that is
~140 GB.

**Release numbers are not sequential counters.** GTDB has run 202, 207, 214, 220, 226,
232 — r226 is one release behind r232, not six.

## Software Versions

Versions are dynamically collected from installed tools. Most use `--version` flag, but TaxonKit uses `version` subcommand.

| Tool | Conda Spec | Purpose |
|------|------------|---------|
| antiSMASH | bioconda::antismash | BGC detection |
| BiG-SCAPE | bioconda::bigscape | GCF clustering |
| GTDB-Tk | bioconda::gtdbtk | Phylogenetic placement |
| TaxonKit | bioconda::taxonkit | Taxonomy processing |
| HMMER | bioconda::hmmer | Coupling enzyme HMM alignment and search |

Version information is output to `results/pipeline_info/software_versions.json`.

## HTML Report Features

The interactive HTML report (`bgc_report.html`) includes:

### Tabs
The report uses 7 tabs:
- **Overview**: Summary statistics grid, rarefaction curve, pipeline resource usage (collapsible) and software versions
- **Phylogeny**: NCBI taxonomy tree + GTDB-Tk phylogenetic tree and BGC distribution
- **Genomes**: Searchable genome table with links to individual genome pages
- **GCF Analysis**: GCF biosynthetic NJ tree (embedded as base64), dynamic coupling enzyme class table, BiG-SCAPE clustering statistics, GCF visualization, and the **pepM identity vs gene-cluster similarity** figure with its correlation table — the Yu et al. replication that is the evidence the GCF assignments above it can be trusted
- **Pipeline Info**: resource usage, **BiG-SCAPE partitioning feasibility**, software versions
- **Novel BGCs**: BGC regions without KnownClusterBlast matches
- **KCB Hits**: Known cluster matches grouped by MIBiG entry

### The All-BGCs Tree Was Removed

The report used to carry a global all-BGCs circular tree beside the family-centre
tree. It is gone, on both partitioned and unpartitioned runs.

**On a partitioned run it was quietly wrong.** `bgc_all_bgcs_tree.py` fills
unmeasured pairs with a constant:

```python
row.append(distances.get(key, 1.0))
```

The merged `distance` table holds only within-partition comparisons by design, so
on Erwiniaceae 23,712 of 55,278 pairs — **42.9%** — were that constant, with no
warning printed (unlike `bgc_gcf_tree.py`, which reports its missing pairs). This
is the same failure that `BIGSCAPE_CENTERS` was built to fix for the GCF tree; the
all-BGCs tree never got the equivalent treatment.

It was still valid unpartitioned, where BiG-SCAPE measures every pair, so gating
it on `params.bigscape_partition` was an option. It was removed outright instead,
so that a figure means the same thing in every run mode.

What replaces it is the **family-centre tree**, whose every centre-to-centre
distance is measured via `BIGSCAPE_CENTERS`, so it means the same thing in both
run modes.

`scripts/bgc_all_bgcs_tree.py` has no caller since `PARTITION_TREES` was removed.
It is kept as a standalone tool and now counts substituted distances, warns, and
refuses above 5% — which is exactly the protection an ad-hoc run against a merged
partitioned database needs.

### Report Tabs

Eight, in this order: Overview, Phylogeny, Genomes, GCF Analysis, **GCF Trees**,
Novel BGCs, KCB Hits, Pipeline. Tabs are pure CSS radio buttons, so adding one
means an `#tabN:checked ~ #contentN` rule in `viz/report_assets.py` alongside the
markup — there is no JavaScript involved in tab switching.

The GCF Trees tab holds the family-centre tree. The
coupling-enzyme class table stays in GCF Analysis: the tree figures carry their
own colour legends, and the table is a classification reference rather than a
tree legend.

### The GTDB-Tk Tree Is Not Drawn in the Report

`prepare_phylo_tree_for_js` prunes the GTDB-Tk tree to the analysed genomes and writes
`pruned_phylo_tree.nwk`, which is the useful output — it opens in iTOL, FigTree or
Dendroscope. Its returned dict reaches the report generator, but **nothing renders it**:
there is no tree renderer in `REPORT_JS`. The Phylogeny tab is titled for the GCF
distribution it actually shows and points at the Newick files.

Adding a rendered tree means writing a renderer, not re-enabling one. `tree_viz.py` was
cut from 1,137 lines to 220 on 2026-08-28; the three circular-tree plotters it used to
hold were dead (exported from `viz/__init__`, called from nowhere) and
`plot_circular_phylogenetic_tree` timed out after two minutes on a real GTDB-Tk tree
because it drove a hand-rolled Newick parser rather than Bio.Phylo. That parser is gone
too — it survived only as an `except` fallback, unreachable in practice since biopython
is a hard dependency of every environment that runs this code. The Bio.Phylo path prunes
20,051 terminals to 285 in about two seconds.

### The Taxonomy Tree Summarises; It Does Not Embed Genomes

The tree exists to answer "how much of this clade carries a BGC" at each rank. Two things
worked against that.

**The prevalence figure was computed and discarded.** `genomes_with_bgcs` was read from the
node stats and never rendered; the headline number was `avg 0.18`, an average over every
genome including the many with none — the same distortion that made "avg BGCs/genome 0.2"
useless on the Overview. Nodes now lead with `285 of 1735 genomes (16.4%)`, then total
BGCs, then a per-BGC-positive average.

**The per-species genome tables were a second copy of the Genomes tab** — 1,735 genome
links, 351 tables, roughly 90% of the tree's 657 KB. They are kept, because they are
useful, but rendered from a JSON payload on first expand rather than inlined.
`renderTaxonomyGenomes()` in `viz/report_assets.py` mirrors `render_genome_list()` in
`viz/taxonomy.py`, including the `greenBg`/`redFont` thresholds — **change both together.**

Note the tree is collapsed on load by a `DOMContentLoaded` handler that hides every
`.node-children` except the first, so inspecting the static HTML is misleading: it shows
the pre-JavaScript state.

Measured on Pantoea: tree block 657,575 to 306,829 chars, inline tables 351 to 0, inline
genome links 1,735 to 0, payload 73 KB for all 1,735 genomes.

For a single-genus run the upper ranks are redundant — domain through family all report
the same 285 of 1,735, because every genome sits in all of them. The tree earns its keep
at broader scope, where those ranks differentiate.

### The Genome Table Renders From JSON

Fully rendered, 1,735 genome rows were 612 KB — the largest single element in the report
— parsed and painted on load although almost nobody scrolls past the first screenful.

`viz/tables.generate_genome_table_html` now returns a dict: `initial_rows` (the first
`INITIAL_GENOME_ROWS`, currently 100, rendered as HTML) and `data_json` (all rows as a
compact array-of-arrays, not objects — repeating six keys 1,735 times is pure overhead).
The page carries the array in a `<script type="application/json">` island and renders
rows from it on demand.

**Search runs over the array, not the DOM**, so it still covers every genome. That is the
point rather than a side effect: `ATCC 35400` and `GCA_963520565.1` each match exactly one
genome, and neither is in the first 100 rows — a DOM-based filter over a truncated table
would silently find nothing.

`renderGenomeRows()` in `viz/report_assets.py` mirrors `_genome_row_html()` in
`viz/tables.py`; **change both together**. A test asserts the server-rendered first row and
the JS-rendered equivalent are identical.

Measured on Pantoea: Genomes panel 612,422 to 194,518 chars, rows in the DOM 1,735 to 100,
whole report 3.62 to 3.23 MB.

### Rarefaction Curve
- Shows GCF discovery saturation across sampled genomes
- Generated from BiG-SCAPE SQLite database (`{taxon}.db`), plus `region_counts.tsv`
- **Coverage is reported as Chao2**, the incidence-based asymptotic richness estimator
  (`chao2()` in `viz/rarefaction.py`): `S_est = S_obs + ((m-1)/m)·Q1²/(2·Q2)`, with Q1/Q2
  the GCFs seen in exactly one/two genomes. Coverage = `S_obs/S_est`. The report shows
  S_obs, S_est, coverage, and Q1/Q2
- The older `saturation` value — a ratio of the curve's early slope to its late slope —
  is still returned and plotted for continuity, but it is **not** a calibrated coverage
  estimate. Benchmarked against known ground truth it is directionally right and accurate
  when genuinely saturated, but optimistic when undersampled (reported 71% against 61%
  true; 39% against 19% true). Cite Chao2, not this
- **The x-axis is every analysed genome, not just BGC-carrying ones.** Genomes with no
  phosphonate BGC produce no region GBK, so they never reach the BiG-SCAPE DB — the
  denominator has to come from `region_counts.tsv`, passed as `--counts`. Without it the
  function falls back to BGC-positive genomes only and labels the axis and report text
  accordingly instead of silently overstating coverage
- Resampling is seeded (`--seed`, default 0) — see the reproducibility note in Known Issues

**Measured on *P. ananatis* (343 genomes, 2026-08-25):** 192 genomes (56.0%) carry a
phosphonate BGC, giving 225 regions in 6 GCFs with incidence 180/23/16/4/1/1. Chao2
estimates 7.00 families → 85.8% coverage (Q1=2, Q2=0, so the no-doubleton branch fires).
The slope-ratio saturation reads 92-97% on the same data — a 6-11 point overstatement,
matching the simulated bias. Note the denominator change moved the slope ratio (96.8% →
92.0%) but left Chao2 unchanged (85.8%): Chao2 depends on the incidence distribution
rather than the axis length, which is a further reason to prefer it.

### GCF Visualization
- Shows representative BGCs for each Gene Cluster Family
- Includes gene arrows with functional annotations
- Color-coded by gene function (core biosynthetic, transport, regulatory, etc.)
- Links to antiSMASH results for detailed analysis (paths relative to `main_analysis_results/{taxon}/`)
- Displays KCB hit or "Potentially Novel" designation for each GCF representative
- Novel BGCs tab shows GCF family assignment when clustering is enabled

### GCF Biosynthetic Tree
- `GCF_BIOSYNTHETIC_TREE` runs **before** `VISUALIZE_RESULTS` — its PNG output is passed as `gcf_tree_png` input to create an explicit Nextflow data dependency
- The tree PNG is embedded as base64 in `bgc_report.html`, making the report self-contained
- Published copies also exist in `gcf_heatmap/` for standalone use
- `conf/conda.config` uses `withName: 'GCF_BIOSYNTHETIC_TREE'` for the conda environment
- Also outputs `phosphonate_itol_coupling.txt` (coupling enzyme class colorstrip) and
  `phosphonate_coupling_support.tsv` (reference support per assignment), both passed to `VISUALIZE_RESULTS`

### Dynamic Coupling Enzyme Class Table
- The coupling enzyme table in the GCF Analysis tab is built dynamically at report-generation time — not hardcoded
- `build_coupling_table_rows()` in `viz/report_sections.py` queries the BiG-SCAPE SQLite DB to determine the dominant coupling class per GCF family, so the table remains accurate regardless of GCF family ID shifts between runs
- Input: `phosphonate_itol_coupling.txt` (from `GCF_BIOSYNTHETIC_TREE`) + BiG-SCAPE DB
- `load_coupling_classes(path, region_only=False)` is the shared parser for this file, defined in `utils/constants.py` and imported by all scripts that read the iTOL coupling colorstrip

### `generate_html_report` Structure

`scripts/visualize_results.py` is an orchestrator (~555 lines, two functions:
`generate_html_report` and `main`). It parses arguments, calls the section builders in
`viz/`, and assembles the page.

**Report content belongs in a `viz/` module, not here.** The file previously grew to
3,300 lines by keeping its own copies of functions that also existed in `viz/`; because
nothing imported the `viz/` versions, the two drifted apart and the package looked live
while being dead. Every section builder is now imported from `viz/`:

| Import from | Provides |
|-------------|----------|
| `viz.tables` | genome table, summary statistics, BGC distribution table |
| `viz.clustering` | BiG-SCAPE statistics, GCF representative visualisation |
| `viz.taxonomy` | interactive taxonomy tree |
| `viz.tree_viz` | `prepare_phylo_tree_for_js` (+ the Newick parse/prune helpers behind it) |
| `viz.distribution` | GCF × genus heatmap, genus-specific and widespread tables |
| `viz.genome_pages` | per-genome HTML pages |
| `viz.rarefaction` | rarefaction curve |
| `viz.report_assets` | `REPORT_CSS`, `REPORT_JS` |
| `viz.report_sections` | `_build_*` blocks, `build_coupling_table_rows` |
| `utils.trace` | trace aggregation + resource-usage HTML |

- CSS and JS live in `viz/report_assets.py` as the module-level string constants `REPORT_CSS` / `REPORT_JS` (plain strings, not f-strings, so braces need no escaping)
- Five helper functions in `viz/report_sections.py` handle the large content blocks, keeping `generate_html_report` to ~260 lines:
  - `_build_kcb_content(kcb_stats, taxon_clean, gcf_data)` → `{kcb_mapping_section, novel_bgcs_tab_content, kcb_hits_tab_content}`
  - `_build_bigscape_overview_cards(gcf_data)` → overview grid HTML
  - `_build_bigscape_section_html(bigscape_stats_html, gcf_visualization_html, taxon_clean)` → GCF Analysis tab section
  - `_build_versions_html(versions_data)` → software versions table
  - `_build_rarefaction_section(rarefaction_stats)` → rarefaction curve block

### BGC Distribution Analysis
- GCF × Genus heatmap showing BGC distribution across taxonomic groups
- Genus-specific GCFs table (potential taxon markers)
- Widespread GCFs table (found in 5+ genera, conserved or HGT)
- Uses GTDB-Tk taxonomy when available, falls back to NCBI taxonomy
- Phylogenetic tree files available in `results/gtdbtk_results/` for external viewers (iTOL, FigTree)

## Post-Pipeline BGC Analysis Scripts

These standalone scripts (in `scripts/`) perform additional analyses after the main pipeline completes. They operate on the BiG-SCAPE SQLite database and antiSMASH outputs.

**Shared building blocks** — extend these rather than re-implementing per script:

| Module | Provides | Used by |
|--------|----------|---------|
| `utils/bigscape_db.py` | `connect`, `fetch_bgc_records`, `fetch_families`, `fetch_domain_set`, `fetch_best_domain_per_cds`, `record_metadata` | pfam / synteny / architecture trees |
| `utils/itol.py` | `write_dataset` plus `write_colorstrip` / `write_binary` / `write_simplebar` / `write_text`, `simple_legend`, `dataset_path` | every script that emits iTOL files |
| `utils/tree_layout.py` | `assign_layout`, `max_depth`, `draw_cladogram` and their circular counterparts | gcf / all-BGCs trees |
| `utils/tree_building.py` | `build_nj_tree(labels, dm_rows)` | all NJ trees |
| `utils/antismash_parser.py` | JSON loading, `build_json_index`, `parse_bgc_label`, `parse_location_bounds`, `find_region_feature`, `region_bounds`, `cds_in_region` | coupling annotation / coupling trees |
| `utils/constants.py` | palettes, `DOMAIN_NAMES` + `domain_name(accession)`, `load_coupling_classes` | synteny / architecture trees, every coupling script |
| `utils/trace.py` | `parse_timestamp`, `aggregate_trace_by_process`, `generate_resource_usage_html`, `generate_gantt_chart_html` | HTML report |
| `utils/plotting.py` | `SVG_METADATA`; pins `svg.hashsalt` and the Agg backend on import | every script that writes an SVG |

`bgc_all_bgcs_tree.py` and `bgc_gcf_tree.py` still query SQLite directly — they read the
`distance` table and family centers, which the shared query set does not cover.

### `scripts/prune_antismash_results.py` — reclaim antiSMASH disk

At 1,735 genomes `antismash_results/` was 41.8 GB, and a BGC-negative genome's directory
is the same size as a BGC-positive one (~23 MB either way). The bulk is not regions:

| file | size | note |
|------|------|------|
| `{genome}.gbk` | 8.6 MB | annotated genome |
| `{genome}.json` | 6.7 MB | **required by `CHECK_ANTISMASH_REUSE`** |
| `{genome}.zip` | 5.6 MB | archive of the very same directory |
| `js`/`images`/`css` | 708 KB | byte-identical in every genome's directory |

So the largest *safe* win is not deleting BGC-negative genomes — it is dropping the
redundant `.zip` from every directory, which loses nothing.

| tier | frees (Pantoea) | reuse still works? |
|------|-----------------|--------------------|
| `archives` | 9.9 GB (24%) | yes — lossless |
| `strip` (default) | 25.4 GB (61%) | **yes** |
| `purge` | more | **no** — antiSMASH re-runs on those genomes |

`strip` keeps `{genome}.json` and `.antismash_meta` for BGC-negative genomes, which is
exactly what `CHECK_ANTISMASH_REUSE` tests for, so `--reuse_antismash_from` still skips
them. `purge` deletes the directory outright: for a low-prevalence taxon that forfeits
most of the compute bill on any later reuse run.

Dry-run by default; `--apply` deletes. It refuses to `--apply` while a Nextflow run is
active (exit 2), because pruning published output races with `publishDir`.

**It prunes only the published copy.** `publishDir` uses `mode: 'copy'`, so `work/` holds
another copy of every result and accumulates one per run — a single genome was found in
`work/` twice plus `results/` once, three copies of 22 MB. Reclaim those with
`nextflow clean -f` once you no longer need `-resume`.

**Not a pipeline stage, deliberately.** `publishDir` re-publishes from `work/` on
`-resume`, silently undoing a prune, and a stage that deletes published output races with
other processes still publishing. It is a post-run tool, like the other `bgc_*.py`
scripts.

### `scripts/bgc_pfam_tree.py` — Jaccard-distance NJ tree of BGCs

Builds a Neighbor-Joining tree based on Pfam domain presence/absence (Jaccard distance).

```bash
python scripts/bgc_pfam_tree.py \
    --db results/bigscape_results/Pantoea/Pantoea.db \
    --bgc_type phosphonate \
    --outdir results/bgc_trees/Pantoea
    [--family_id 2]    # Optional: restrict to a single GCF
```

Outputs: `_pfam_tree.nwk`, `_jaccard_distances.tsv`, `_domain_matrix.tsv`, `_metadata.json`

### `scripts/bgc_synteny_tree.py` — LCS-based gene-order tree

Builds a tree based on ordered domain sequences (one domain per CDS, sorted by genomic position). Uses normalized LCS distance. **Note:** Can be confused by strand orientation.

```bash
python scripts/bgc_synteny_tree.py \
    --db results/bigscape_results/Pantoea/Pantoea.db \
    --bgc_type phosphonate \
    --outdir results/bgc_trees/Pantoea/GCF2_synteny \
    [--family_id 2]
```

Outputs: `_synteny_tree.nwk`, `_domain_sequences.tsv`, `_lcs_distances.tsv`, `_metadata.json`

### `scripts/bgc_architecture_tree.py` — Architecture deduplication tree

Groups BGCs by exact domain multiset (orientation-independent), then builds a generalized Jaccard NJ tree of the unique architectures. Best for within-GCF comparison.

```bash
python scripts/bgc_architecture_tree.py \
    --db results/bigscape_results/Pantoea/Pantoea.db \
    --bgc_type phosphonate \
    --family_id 2 \
    --outdir results/bgc_trees/Pantoea/GCF2_arch
```

Outputs: NJ tree + five iTOL annotation files (count bar, domain binary, genus colorstrip, arch label, genome list).

Architecture labels: `arch_001_n138` (arch rank, count). GCF2 phosphonate → 205 BGCs → 21 unique architectures; arch_001 (n=138) is the dominant core.

### `scripts/bgc_coupling_annotation.py` — iTOL coupling enzyme colorstrip

Classifies each phosphonate BGC by the coupling enzyme acting on phosphonopyruvate (the branching step immediately downstream of PEP mutase). Reads antiSMASH JSON files for rich SMCOG and rule-based-cluster annotations, then outputs an iTOL DATASET_COLORSTRIP file.

```bash
python scripts/bgc_coupling_annotation.py \
    --antismash_dir results/antismash_results/Pantoea \
    --metadata results/bgc_trees/Pantoea/phosphonate_metadata.json \
    --outfile results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \
    --bgc_type phosphonate
```

**⚠️ The GCF numbers below are run-specific and will not match your output.**
Counts are also from that older run. For reference, the Pantoea genus run of
2026-08-25 (1,735 genomes, 320 regions) gave: Synthase 236, Reductase 29,
Decarboxylase 49 (27 of which were then split off as a separate
Decarboxylase-Nucleotidyltransferase class — see below), Transaminase 6,
Unknown 0 — the zero being the segment-based membership fix.
`family.id` in the BiG-SCAPE database is an `INTEGER PRIMARY KEY AUTOINCREMENT` — it
records the order families happened to be written, not a stable biological identity.
It shifts between runs, between taxa, and whenever the genome set changes: the table
below is from a Pantoea run with at least 8 families, while a *P. ananatis* run
(2026-08-25, 343 genomes) produced only 6. Treat the **class → pathway → marker**
columns as the durable content and the GCF column as an example. Within a single run
`family.center_id` is a better handle, since it points at the BGC record serving as
the family center rather than at write order.

The rendered report does not have this problem: `build_coupling_table_rows()` re-derives
the dominant coupling class per family from the database at report time, so the HTML
table self-corrects. Only the hard-coded numbers in this file go stale.

**Coupling enzyme classes detected (Pantoea, n=1212 BGCs):**

| Class | Marker | Pathway | GCF | Count |
|-------|--------|---------|-----|-------|
| FrbC | SMCOG1271 (HMGL-like) | → phosphonomethylmalate → phosphinothricin-type | GCF-2/3 | 920 |
| Fe-ADH | Fe-ADH rule | → phosphonolactate (reductase route) | GCF-4/6 | 112 |
| Ppd | SMCOG1055 (ThDP) or TPP_enzyme_C | → 2-phosphonoacetaldehyde → 2-AEP | GCF-1/5/8 | 156 |
| PalB | SMCOG1019 (Aminotran_1_2/PF00155) | → phosphonoalanine | GCF-7 | 20 |
| Unknown | — | — | — | 4 |

**Region boundary reading:** CDS scanning is limited to the region feature's extent, read
with `parse_location_bounds(..., span=False)` — the first coordinate pair only. This is
deliberately narrower than what `bgc_coupling_tree.py` uses (see Known Issues); widening
it reclassifies BGCs whose region feature has a compound location.

**Key insights:**
- Classification maps almost perfectly onto BiG-SCAPE GCF families — coupling enzyme type is the primary determinant of GCF membership.
- The `Fe-ADH` rule-based marker (iron-containing alcohol dehydrogenase / 2-Hacid_dh_C) is antiSMASH's marker for the phosphonopyruvate reductase (→ phosphonolactate) pathway.
- AEP-pathway BGCs (GCF-1/8 in that run) use Ppd as coupling enzyme regardless of tailoring enzymes downstream.

**The Decarboxylase-Nucleotidyltransferase class was removed on 2026-09-11.** From
2026-08 to 2026-09 a fifth class split off the decarboxylases on `TPP_enzyme_C` plus an
`NTP_transf_3`/`NTP_transf_2` hit, on the theory that a cytidylyltransferase in the BGC
marked the CDP-activated phosphonolipid route. Checked against the two clusters whose
chemistry is known from lab work, it came out **inverted**:

| Cluster | Lab truth | Class assigned | Margin |
|---|---|---|---|
| *P. ananatis* LMG 5342 `HE617160.1.region002` | phosphonolipid | Decarboxylase | 0.0 |
| *Winslowiella iniecta* B149 `JRXF01000012.1.region001` | **not** a lipid | Decarb-Nucleotidyltransferase | 0.0 |

The confirmed lipid carries no NTP_transf at all — several of its biosynthetic CDS are
unannotated in that assembly, so the marker simply is not visible — while the confirmed
non-lipid carries one, activating a substrate for some other energetically unfavourable
step. NTP transfer is generic activation chemistry, not a lipid signature.

The margin is structural, not incidental: both classes scored against the same DhpF /
Fom2 / Ppd reference set (via the `_SHARED_REFS` mechanism, now also gone), so percent
identity was identical for both by construction and `margin` was always exactly 0.0. No
amount of added sequence evidence could have separated them.

`LEGACY_CLASS_NAMES` in `utils/constants.py` still maps the old `TPP+NTP` and `Ppd-CDP`
spellings, now onto `Decarboxylase`, so annotation files from those runs still load.

**Predicting phosphonolipid vs. small molecule: the failures were measurement, the
open question is real.** Four signals were tried — coupling class, `NTP_transf_3` copy
number, CDP-alcohol phosphatidyltransferase proximity, and TIGRFAM — and none separated
the two characterised examples. **All four were measured on a cluster the pipeline could
not see.** *P. ananatis* LMG 5342's 2012 deposit annotates 6 of its 15 genes, and
antiSMASH runs gene finding only on records with ZERO CDS, so the other 9 — including
the class-V transaminase and both CDP-alcohol phosphatidyltransferases — were invisible
to every metric built on top. See `RECOVER_ORFS`.

With the genes restored the question is still open, but for a better reason. The two
clusters are now **nearly identical in domain content**:

| | LMG 5342 r2 (**lipid**) | *Winslowiella* B149 (**not** lipid) |
|---|---:|---:|
| aepZ-family transaminase (PF00266) | 1 | 2 |
| NTP_transf_3 (PF12804) | 1 | 1 |
| CDP-alcohol phosphatidyltransferase (PF01066) | **2** | 1 |
| pepM / Ppd | 1 / 1 | 1 / 1 |

Copy number is the only domain-level difference, and at one example per class that is
not signal. So whatever distinguishes them is **not visible in Pfam content** — the
acceptor specificity of the CDP-alcohol phosphatidyltransferase, substrate availability
and regulation are the places left to look.

Note what this means for the original hypothesis that NTP_transf plus a CDP-alcohol
phosphatidyltransferase marks a phosphonolipid: it is **supported** by the confirmed
lipid, which carries both. It is simply not discriminating. It was rejected three times
on evidence that could not see those genes at all.

**PalB detection was corrected on 2026-08-25.** It previously used SMCOG1013
(Aminotran_3, fold type IV PLP), which is a different aminotransferase class from
PalB — an AAT superfamily enzyme (fold type I PLP), annotated Aminotran_1_2 / PF00155
/ SMCOG1019. The marker is now SMCOG1019 in both `bgc_coupling_annotation.py` and
`bgc_coupling_tree.py`'s `CLASS_MARKERS`.

**Measured impact on the Pantoea genus run: 5 of 6 Transaminase calls were false
positives.** Counts went Transaminase 6 → 1, Unknown 0 → 5; no other class moved. The
one surviving call is `CEUYZP010000005.1.region001` (*Pantoea* sp. E956-1_S3, locus
`ctg5_2`), annotated `SMCOG1019: aminotransferase`. The five reclassified BGCs carry
an Aminotran_3 enzyme that is not the coupling enzyme, so `Unknown` is the honest
label; the coupling enzyme trees are the way to resolve what they actually are.

Note `PF00155` never appears literally in antiSMASH JSON — it writes the domain *name*
`Aminotran_1_2`. That domain was tested as an additional fallback and made no
difference to the outcome (the five reclassified BGCs do not carry it), so it was left
out: `Aminotran_1_2` hits 82 of 313 phosphonate regions and is too promiscuous to use
as a coupling-enzyme marker on its own.

The remaining caveat from the original note still stands: coupling enzymes are not
always adjacent to pepM (the phosphonoalamide BGC places PalB far from it). Region
membership is now segment-based, so the whole region is scanned regardless of distance.

**Note on PalA:** PalA (phosphonopyruvate hydrolase, a phosphonate degradation/resistance gene) does not confound the classification — all GCF types show clear biosynthetic markers.

### `scripts/bgc_coupling_tree.py` — Coupling enzyme phylogenetic trees

Builds FastTree ML trees (LG model) for pepM (Tree A) and per-class coupling enzymes (Tree B), with characterized reference sequences as phylogenetic anchors. **Standalone only — no longer a pipeline stage** (removed 2026-08-27; see "Why There Is No pepM Tree"). Consumes `phosphonate_metadata.json` and `phosphonate_itol_coupling.txt` from `GCF_BIOSYNTHETIC_TREE`, and uses `load_coupling_classes` from `utils/constants.py`. `hmmer` and `fasttree` are no longer in any pipeline conda environment, so install them yourself and pass `--hmmbuild`/`--hmmalign`/`--hmmsearch`/`--fasttree`. Use it to place an individual ambiguous enzyme phylogenetically — something the scalar reference-support score cannot do.

**HMM strategy (4 steps per class):**
1. `hmmbuild` from a single seed reference → initial HMM
2. `hmmalign` all references to initial HMM → aligned references
3. `hmmbuild` from aligned references → refined HMM
4. `hmmalign` refs + query sequences → final alignment → FastTree ML tree (LG model)

**Sequence extraction (annotation-first with HMM fallback):**
- Primary: antiSMASH annotation markers (SMCOG/domain hits from `gene_functions` / `sec_met_domain`) — zero extra compute, already in JSON
- Fallback: for any BGCs the annotation missed, extract all CDS from the region and run `hmmsearch` against the class reference HMM; select the highest-scoring hit per BGC
- This handles divergent sequences that escape SMCOG thresholds (ThDP decarboxylases that carry `TPP_enzyme_C` but not `SMCOG1055`)

**CLASS_MARKERS — extraction markers per class:**

| Class | Marker type | Marker | Rationale |
|-------|-------------|--------|-----------|
| Synthase (FrbC-like) | smcog | SMCOG1271 | HMGL-like phosphonomethylmalate synthase |
| Decarboxylase (Ppd-like) | domain | TPP_enzyme_C | Divergent ThDP decarboxylases carry this but not SMCOG1055 |
| Reductase (VlpB-like) | domain | Fe-ADH | Phosphonopyruvate reductase (iron-containing ADH) |
| Transaminase (PalB-like) | smcog | SMCOG1019 | Aminotran_1_2/PF00155, AAT superfamily (corrected 2026-08-25) |

**Tree outputs per class:** `{class}_tree.nwk` + four iTOL annotation files (coupling class colorstrip, GCF colorstrip, source colorstrip, organism text labels).

**Reference sequences** (`assets/reference_sequences/`):

| FASTA ID | Protein | Function | Source |
|----------|---------|----------|--------|
| `BGC0000904\|ABB90393\|FrbD` | FrbD | PEP mutase | *Streptomyces rubellomurinus* (FR-900098) |
| `BGC0000904\|ABB90392\|FrbC` | FrbC | phosphonomethylmalate synthase | *Streptomyces rubellomurinus* |
| `BGC0000897\|ACZ13456\|DhpE` | DhpE | PEP mutase | *Streptomyces luridus* (Dehydrophos) |
| `BGC0000897\|ACZ13457\|DhpF` | DhpF | phosphonopyruvate decarboxylase | *Streptomyces luridus* |
| `BGC0000938\|ACG70831\|Fom1` | Fom1 | PEP mutase | *Streptomyces fradiae* (Fosfomycin) |
| `BGC0000938\|ACG70832\|Fom2` | Fom2 | phosphonopyruvate decarboxylase | *Streptomyces fradiae* |
| `BGC0000806\|AHL24479\|PepM` | PepM | PEP mutase | *Glycomyces* sp. NRRL B-16210 |
| `BGC0000806\|AHL24480\|Ppd` | Ppd | phosphonopyruvate decarboxylase | *Glycomyces* sp. |
| `Phosphonoalamide_BGC\|WP_030764868\|PnaD` | PnaD | PEP mutase | *Streptomyces* sp. NRRL B-2790 |
| `Phosphonoalamide_BGC\|WP_051781701\|PnaA` | PnaA | phosphonopyruvate transaminase | *Streptomyces* sp. NRRL B-2790 |
| `Valinophos_BGC\|WP_031174023\|VlpA` | VlpA | PEP mutase | *Streptomyces durhamensis* NRRL B-3309 |
| `Valinophos_BGC\|WP_063765859\|VlpB` | VlpB | phosphonopyruvate reductase | *Streptomyces durhamensis* |
| `Pantaphos_BGC\|WP_013027161\|HvrA` | HvrA | PEP mutase | *Pantoea ananatis* LMG 5342 |
| `Pantaphos_BGC\|WP_013027159\|HvrC` | HvrC | phosphonomethylmalate synthase | *Pantoea ananatis* LMG 5342 |

**PalB-like detection uses SMCOG1019** (Aminotran_1_2 / PF00155, AAT superfamily,
fold type I PLP) as of 2026-08-25. It previously used SMCOG1013 (Aminotran_3, fold
type IV) — see the correction note above.

**Standalone usage:**
```bash
python scripts/bgc_coupling_tree.py \
    --antismash_dir  results/antismash_results/Pantoea \
    --metadata       results/main_analysis_results/Pantoea/gcf_heatmap/phosphonate_metadata.json \
    --coupling_annotation <itol_coupling_colorstrip.txt> \
    --ref_pepm_faa   assets/reference_sequences/reference_pepM.faa \
    --ref_coupling_faa assets/reference_sequences/reference_coupling_enzymes.faa \
    --outdir         results/main_analysis_results/Pantoea/coupling_enzyme_trees \
    --hmmbuild       $(which hmmbuild) \
    --hmmalign       $(which hmmalign) \
    --hmmsearch      $(which hmmsearch) \
    --fasttree       $(which FastTree) \
    --tree           both    # "A", "B", or "both"
```

### `scripts/bgc_gcf_heatmap.py` — GCF × Species presence/absence heatmap

Generates a heatmap of GCF membership across organism groups, with a GTDB-Tk phylogenetic tree as column ordering and a Jaccard/complete-linkage row dendrogram matching BiG-SCAPE's clustering algorithm. Uses `load_coupling_classes` from `utils/constants.py`.

```bash
python scripts/bgc_gcf_heatmap.py \
    --metadata            results/bgc_trees/Pantoea/phosphonate_metadata.json \
    --coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \
    --gtdbtk_tree         results/gtdbtk_results/Pantoea/gtdbtk_output/classify/gtdbtk.bac120.classify.tree.1.tree \
    --gtdbtk_summary      results/gtdbtk_results/Pantoea/gtdbtk_output/gtdbtk.bac120.summary.tsv \
    --outdir              results/bgc_trees/Pantoea
```

- **Data source**: Only region-level BGC records with GCF assignments at `cutoff=0.3` (303 phosphonate BGCs; sub-records like cand_cluster/protocluster are excluded)
- **True singletons**: Single-member GCFs (size=1), not unassigned records
- **Row dendrogram**: `scipy.spatial.distance.pdist(metric='jaccard')` + `linkage(method='complete')` — matches BiG-SCAPE's clustering algorithm
- **Column tree**: GTDB-Tk phylogenetic tree pruned to representative genomes per organism group, rendered as a cladogram
- **Outputs**: `gcf_species_heatmap.png` and `.svg`

### `scripts/bgc_itol_annotations.py` — iTOL annotation files

Generates iTOL annotation files from bgc_pfam_tree.py or bgc_synteny_tree.py outputs.

```bash
python scripts/bgc_itol_annotations.py \
    --treedir results/bgc_trees/Pantoea \
    --bgc_type phosphonate
```

Outputs: `_itol_gcf.txt` (color strip), `_itol_domains.txt` (binary), `_itol_domaincount.txt` (bar chart).

### Key Pfam accessions for phosphonate BGCs

Verified from antiSMASH clusterhmmer output on Pantoea phosphonate clusters:

| Pfam | Name | Function |
|------|------|----------|
| PF13714 | PEP_mutase | PEP mutase (pepM/aepX) — hallmark gene |
| PF00296 | HMGL-like (HEPD) | 2-hydroxyethylphosphonate dioxygenase |
| PF00682 | FrbC-like (PmmS) | Phosphonomethylmalate synthase (HMGL superfamily) |
| PF02775 | ThDP_C | Phosphonopyruvate decarboxylase |
| PF00266 | Aminotrans_V | 2-AEP transaminase |
| PF13649 | Radical_SAM | Radical C–P chemistry |

**Note on HMGL annotation:** AntiSMASH/BiG-SCAPE annotates phosphonomethylmalate synthase as `PF00682 (HMGL-like)` because it structurally belongs to the HMGL superfamily. The antiSMASH JSON provides richer context via `gene_functions: biosynthetic-additional (smcogs) SMCOG1271: 2-isopropylmalate synthase` and `sec_met_domain: HMGL-like`. BiG-SCAPE only stores the Pfam accession and bit score — no SMCOG or functional description.

### Data Sources

- BiG-SCAPE DB `hsp` table: Pfam accession + bit_score per CDS (populated by antiSMASH clusterhmmer)
- AntiSMASH JSON: richer annotations including `gene_functions`, `sec_met_domain` (SMCOG hits, TIGRFAM), and `product`
- Domain sequences in TSV come from the BiG-SCAPE DB (best Pfam hit per CDS, ordered by `nt_start`)

## Development Notes

### Configuration

- **Conda environments**: Defined centrally in `conf/conda.config` (not in individual modules)
- **Resource labels**: Defined in `conf/labels.config`, applied via process labels in modules
- **SLURM overrides**: Profile-specific adjustments in `nextflow.config`

### Utilities

- `Utils.sanitizeTaxon(name)`: Sanitize taxon for filesystem paths (removes special chars)
- `Utils.antismashParamsHash(params)`: Generate MD5 hash of antiSMASH parameters for reuse tracking
- `Utils.buildReusePath(params, projectDir, tool, taxon, subPath)`: Build absolute path for result reuse
- `Utils.isValidInput(input)`: Check if input is valid (not a placeholder)
- `Utils.optArg(flag, input)`: Build `"--flag path"` for a real input, `""` for a placeholder — use this instead of comparing against a specific `NO_*` name, which silently passes the sentinel through when the names drift apart

Workflow-level helpers live in `subworkflows/helpers.nf` and are included like processes
(`include { placeholder } from './helpers'`), since Nextflow functions are file-scoped:

- `placeholder(name)`: Value channel holding a sentinel file for an optional input
- `clusteringEnabled(method)`: `params.clustering == method`
- `batchSize()`: `params.task_batch_size` coerced to Integer

### Module Guidelines

- Use `publishDir` for outputs, `storeDir` for database downloads
- Use appropriate labels: `process_low`, `process_medium`, `process_high`, `process_high_memory`
- Use `tolerant` label for per-genome processes where individual failures are acceptable
- antiSMASH uses `cache 'lenient'` for directory inputs
- COLLECT_VERSIONS searches `work/conda/` for installed tool versions
- BiG-SCAPE database (`bigscape_db`) is passed explicitly through pipeline for rarefaction curve generation
- Build optional arguments with `Utils.optArg('--flag', input)` rather than comparing a
  staged file against a specific sentinel name. Placeholders are only guaranteed to start
  with `NO_`; hardcoding `!= 'NO_FILE'` passes the sentinel through as a real path once
  the workflow emits a differently-named one

### Testing

```bash
bash tests/run_tests.sh          # full suite (~1 min)
TEST_SCRATCH=/tmp/t bash tests/run_tests.sh   # keep the scratch dir for debugging
BATCH_SIZE=2 bash tests/run_tests.sh          # exercise a different batch size
```

Everything runs in a scratch directory, never in `results/`, so a test run cannot
truncate `pipeline_info/`. The suite covers:

- `Utils` helpers — `optArg` across real/placeholder/list/empty/null inputs, `isValidInput`, `sanitizeTaxon`, and that `batchSize()` yields a real Integer (`collate()` silently fails otherwise)
- Batched per-genome processes — assembly-ID pairing across a batch, one output per genome after `.flatten()`, a corrupt genome skipped without losing its batch, and reuse-copy fidelity for hidden and nested files
- Static checks — every script compiles, and `check_undefined.py` finds calls to names that are never defined or imported (this is what surfaced the phylo-fallback `NameError`)

Notes on the runner: Nextflow derives `projectDir` from the entry script's location, so
the test scripts are staged into the scratch dir with a `scripts/` symlink — otherwise
modules would look for `tests/scripts/...`. GenBank→FASTA needs biopython; the runner
borrows an interpreter that has it (system python or a cached conda env) and skips those
two assertions if none is available.

### Parameter validation: the schema is the contract

`nextflow_schema.json` declares every parameter, and `main.nf` calls nf-schema's
`validateParameters()` before anything else. Two silent failures become errors:

| you type | before | now |
|---|---|---|
| `--taxn Pantoea` | runs *Erwiniaceae*, reports success | `* --taxn (Pantoea): False schema always fails` |
| `--run_gtdbtk false` | **enables** GTDB-Tk | `Value is [string] but should be [boolean]` |

The second is the nastier one. A command-line param arrives as a **string**, and every
non-empty string is true in Groovy, so `--run_gtdbtk false` reads as enabled. Measured on
Nextflow 26.04.3: `--run_gtdbtk false` gives `String "false"` and takes the TRUE branch,
while `-params-file {"run_gtdbtk": false}` gives `Boolean false`. Every gate here is
`if (params.x)`, so it cost two real mistakes — a run that spent 373 CPU-min and 93 GB on
GTDB-Tk after being told not to, and an A/B that would have screened both arms while
reporting them as screen-on against screen-off. Nextflow does not catch either on its own,
and `NXF_ENABLE_STRICT=true` does not change that.

**Adding a param means adding it to the schema**, or the pipeline rejects it at runtime
for everyone. `tests/check_schema.py` compares the schema against `nextflow config -flat`
and fails on drift in either direction; `--write` regenerates it. The generator reads the
param list from Nextflow rather than by parsing the config text — a regex missed one
param, and because the same regex checked its own output the gap stayed invisible until
a run failed.

Three things worth knowing about the schema's shape. All 53 properties sit at the **root**,
not in `$defs` groups: `additionalProperties` only sees properties declared in the same
schema object, so an `allOf`/`$defs` layout rejects every grouped param instead of only
unknown ones. And `"False schema always fails"` is what an unknown parameter looks like —
the message comes from the JSON-schema library, and it names the offending flag. And a param
whose documented "off" value is `null` needs `["integer","null"]` rather than the type its
default implies — the generator cannot infer that from a default of `10`, so
`antismash_phosphonate_neighbourhood` is listed in `NULLABLE` in `tests/check_schema.py`
and `--write` preserves it. Without that, the one setting `nextflow.config` tells you to
use would be rejected.

Disabling something on the command line is no longer possible; use a params file:

```bash
echo '{ "run_gtdbtk": false, "pepm_prescreen": false }' > off.json
nextflow run main.nf -params-file off.json --taxon "Pantoea ananatis"
```

### Task Batching

Steps whose per-genome work is under a couple of seconds are batched — one job per
genome is almost entirely scheduler overhead, and at 3M genomes the submission rate
limit becomes the bottleneck rather than the compute.

Batched processes: `RENAME_GENOMES`, `GENBANK_TO_FASTA`, `COPY_ANTISMASH_RESULT`
(`params.task_batch_size`, default 100). `CHECK_ANTISMASH_REUSE` is not batched because
it already runs with `executor 'local'`.

When adding or changing a batched process:

- **Flatten downstream.** A batched process emits one list per task; consumers that work
  per genome need `.flatten()` (see `DOWNLOAD_GENOMES.out.renamed_genomes`)
- **Keep failures per-genome.** The batch script must catch per-item errors and continue,
  exiting non-zero only if every item failed — otherwise batching turns one bad genome
  into 100 lost ones
- **Watch for input name collisions.** NCBI genomes are all named `genomic.gbff`, so
  `RENAME_GENOMES` stages them as `genome?.gbff` and pairs staging order against the
  assembly-ID list via a manifest
- **`collate()` needs a real Integer.** Params given on the command line arrive as
  strings, which silently fail to dispatch — always go through `batchSize()`

### pepM Pre-Screen (`--pepm_prescreen`)

**Off by default.** It *removes* genomes from the analysis, so turn it on
deliberately — the same posture as `--bigscape_partition`.

Every phosphonate BGC carries a PEP mutase, so a genome without one cannot hold
what this pipeline looks for. Establishing that costs **~0.9 CPU-s** against
antiSMASH's **41.4**, which is what makes an order-scale run tractable.

```
RENAME_GENOMES -> PEPM_PRESCREEN -> ANTISMASH   (only genomes that pass)
```

**Two modes, one decision.** NCBI GenBank annotation is inconsistent — 794 of
2,771 Erwiniaceae genomes (28.7%) carry no CDS translations at all:

| genome | mode | cost |
|---|---|---:|
| annotated | `diamond blastp` over its proteins | 0.234 CPU-s |
| unannotated | `diamond blastx` over its contigs | 2.030 CPU-s |

Both search the same seven references at the same threshold. On the 341 genomes
where both could run they agreed on **100.0%** of calls at bitscore 100 — same
tool, same references, same cutoff, only the input representation differs. That is
why this branch is safe where the others in this codebase were not: the two paths
are verified to make identical decisions, rather than merely intended to.

**blastp is 9.2x cheaper** (measured; an earlier estimate of 14x was optimistic),
which takes the screen from 2.11 to 0.88 days at a million genomes.

**Validated over all 2,771 Erwiniaceae genomes against the real BGC calls:**

| | |
|---|---:|
| sensitivity | **298 / 298** |
| false positives | 8 of 2,473 |
| retained | 306 (11.0%) |
| true-positive bitscore | 154-552 |

**The separation is wide but not absolute, and an earlier version of this table said
otherwise.** It quoted "background bitscore <= 51", which cannot be true of the full
validation: the 8 false positives *are* negatives scoring above the cut of 100. That
figure came from a 60-genome unannotated sample (true positives 154-552, negatives
topping out at 51) and was wrongly carried over as a property of all 2,473 negatives.

Those 8 are not misfires. They are genomes carrying a credible PEP mutase with no
assembled cluster around it — antiSMASH examines them and correctly reports nothing,
which is the behaviour a screen should have at its margin. The bound that matters is
the one on the other side: **no true positive scored below 154**, so the cut at 100 has
54 points of headroom against a miss, which is the direction that loses data.

A `--min_density` guard (500 proteins/Mb) routes partially-annotated genomes to
blastx, closing the one failure mode the two modes do not share. Observed density
was 576-1,015 with nothing below 500, so it costs nothing today.

**A pseudogene-flagged pepM was invisible to the screen, and is no longer.** NCBI
withholds `/translation` from any CDS it flags `/pseudo`, and `parse_genome` collected
only CDS that had one — so a genome whose pepM is annotated `phosphoenolpyruvate mutase`
*and* `/pseudo` reached diamond with no pepM in its protein set and scored **0.0**. Found
on the held-out clade (below), where it cost 2 of 98 true positives; `--min_density` does
not catch it, because both genomes run 726-763 CDS/Mb and density is a whole-genome proxy
for a single-gene problem. Such a CDS is now translated from its own coordinates, internal
stops kept as `X` rather than truncating, since a pseudogene spreads its signal across the
frameshift. Erwiniaceae is unchanged in every field.

### Held-out clade: *Bacteroides fragilis* (2026-09-19)

Every earlier validation was on a clade the reference set draws from — Erwiniaceae
supplies HvrA (1 of 7), the actinomycete set supplies the other 6 and contains the source
strain of one, which self-matched at 828. **Bacteroidota supplies none**, and is held out
in BGC space too: all 5 reference clusters sit 0.82-0.93 from their nearest *B. fragilis*
BGC, none inside the 0.30 cutoff.

136 genomes, ground truth from an unscreened arm: **98 BGC-positive (72%), 143 regions,
19 GCFs**. After the pseudogene fix:

| | |
|---|---:|
| sensitivity | **98 / 98** |
| false positives | 2 of 38 |
| true-positive bitscore | 342-567 |
| top negative | 330 |

**The classes separate completely** — any cut in (330, 342] gives 98/98 with zero false
positives, where the actinomycete set had no such cut. The caveat is diversity rather
than count: one species, true positives clustered at a modal 514, so 98 positives is not
98 independent tests.

At 72% BGC-positive the screen saves only 11% here (277.4 -> 247.3 CPU-min for a screen
costing 3.6), against 6.5x on Erwiniaceae and 7.3x on the actinomycetes. **The saving is
proportional to how dilute the taxon is**, and this clade was chosen to test sensitivity,
not savings. Data: `docs/comparisons/pepm_prescreen/heldout_bacteroides/`.

**Expanding the reference set made it worse.** Mining MIBiG by HMM added eight
unique pepMs (15 total); sensitivity stayed at 298/298 while false positives rose
from 8 to 67. Those extras are pepMs from fosfomycin and dehydrophos clusters,
divergent enough to attract spurious matches without catching anything new. Note
the curated seven include a *Pantoea* pepM (HvrA, Pantaphos), which favours this
test set — the expansion is unproven rather than useless, and worth revisiting for
a taxonomically distant clade.

`prescreen_results/<taxon>/prescreen_*.tsv` records every genome with its mode,
CDS density, best bitscore and verdict, so what was skipped is auditable rather
than silent.

### Validation Matrix (2026-09-09)

Five configurations, all on Erwiniaceae, all reproducing the same GCF network.
**Comparisons are by nucleotide sequence and co-membership, not by counts** —
matching totals can hide a substitution.

| Run | Genomes to antiSMASH | BGCs | Families | vs baseline |
|---|---:|---:|---:|---|
| screen off, unpartitioned *(baseline)* | 2,771 | 333 | 19 | — |
| screen off, partitioned | 2,771 | 333 | 19 | ARI 1.0000 |
| screen on, unpartitioned | **306** | 333 | 19 | ARI 1.0000 |
| screen on, partitioned | **306** | 333 | 19 | ARI 1.0000 |
| screen on + antiSMASH reuse (*P. ananatis*) | **0** | 225 / 225 | 6 | 0 lost, 0 extra |

Every ARI comparison is over 23,995 co-membership pairs with **0 split and 0
merged**. The screened runs recovered all 333 BGCs with **333/333 identical
nucleotide sequences**, 5,845,237 bases either way.

**antiSMASH output is not byte-reproducible.** BiG-SCAPE's `gbk.hash` differed on
all 333 BGCs between two runs that were otherwise identical, because antiSMASH
stamps `Run date` into every region GenBank. Compare `nt_seq`, never file hashes.

**Reuse and the screen compound.** The *P. ananatis* run downloaded 343 genomes,
screened 192 through, and ran antiSMASH **zero** times — every screened genome
already had a result under Erwiniaceae. All 225 of its BGCs matched the
Erwiniaceae subset base-for-base (3,526,555 bases). The screen filters
`renamed_genomes` before the reuse branch, so both paths consume the narrowed set
and there is no second code path to keep in step.

The partitioned runs split 236/88/4/2/2/1 and `BIGSCAPE_CENTERS` measured
**171 of 171** centre pairs, so the family-centre tree has a fully measured
backbone in every configuration.

### Wall Time at Scale: antiSMASH Batching and GTDB-Tk Sharding

Elapsed time for a large run is set by two stages; everything else is under a day
combined.

| Stage | before | after | mechanism |
|---|---:|---:|---|
| antiSMASH | 34.7 d | 0.7 d submit / 2.4 d compute | 50-genome batches |
| GTDB-Tk | up to 23 d | ~2.3 d at 10 shards | 5,000-genome shards |

**antiSMASH was submission-bound, not compute-bound.** At one task per genome and 20
submissions/min, a million genomes spends 34.7 days being *submitted* against 2.4 days
computing. `params.antismash_batch_size` (default 50) puts submission at 0.7 days,
comfortably under compute; larger batches buy nothing and only coarsen retry granularity.

The batch loop **continues past a failed genome** rather than exiting. That mattered less
when a failure cost one genome; batched, an aborting task would forfeit 50. The task exits
non-zero only when *every* genome in the batch failed, which signals a broken environment
(the conda startup race, a missing database) rather than bad input — and that is what the
retry in `conf/labels.config` is for. `time` is raised to 8h since 50 genomes run
sequentially.

Results are written under `as_out/` so the output glob cannot match the staged database
directory, and `saveAs` strips the prefix so the published layout is unchanged.

**GTDB-Tk is sharded, and its tree is gone.** `classify_wf` still builds a tree internally
— pplacer placement is how it classifies — but a per-shard tree spans a disjoint genome
set, and N of them cannot be concatenated into one phylogeny. `MERGE_GTDBTK` concatenates
the summaries instead, which is lossless (one independent row per genome) and is what
every consumer actually reads. It fails if any genome appears in two shards, since that
would silently inflate every per-clade count.

Removing the tree touched more than the tree:

- **The GCF x genus heatmap lost its phylogenetic column ordering, and it has been
  restored from a better source.** The order used to come from the run's own pplacer
  tree; it now comes from GTDB's *reference* phylogeny, shipped in the GTDB-Tk data
  package at `pplacer/gtdb_r<rel>_bac120.refpkg/gtdb_r<rel>_bac120_decorated_unrooted.tree`.
  That is curated, identical between runs, and independent of which genomes happened to
  be sequenced — where a pplacer order could shift with the query set. See
  `genus_tree_from_gtdb_reference()`.
- **GTDB-Tk reuse would have broken silently.** `CHECK_GTDBTK_REUSE` required a tree file
  to exist before returning REUSE; with no run producing one, every reuse would have
  fallen back to a full re-run. It now checks the summary alone.
- `prune_tree` in `filter_gtdbtk_results.py`, `--outgroup` / `params.gtdbtk_outgroup`, and
  `prepare_phylo_tree_for_js` are all gone or orphaned — nothing produces a whole-set tree
  to root, prune, or render. `viz/tree_viz.py` itself is retained but no longer imported
  by `viz/__init__`.

**Batching antiSMASH invalidates every cached antiSMASH task**: the process source and its
input cardinality both change, so the first run after this costs a full re-analysis.

### Report JavaScript Is Not Covered by the Python Checks

`tests/check_undefined.py` parses Python with `ast`, but the report's JavaScript is
Python *string data* — `REPORT_JS` in `viz/report_assets.py` plus inline fragments in
`viz/clustering.py` and `visualize_results.py` — so `ast` sees opaque text. That blind
spot shipped a Genomes-tab search box wired to `filterGenomes()`, a function defined
nowhere: every keystroke threw a `ReferenceError` and filtered nothing, silently, for
the life of the feature.

`utils/report_lint.py` closes it. `check_report(html)` cross-references inline
`on*="name(...)"` handlers against `function name(` definitions and returns readable
problems. `visualize_results.py` calls it **before writing** the file and exits 1
rather than emitting a report with dead handlers — verified: breaking a function name
gives exit 1 and leaves any existing report untouched.

The linter runs against the assembled HTML because that is the only point where all
the JS fragments exist together; checking the Python sources individually would report
false positives, since a handler defined in one fragment is called from another.

`tests/check_report_js.py` self-tests the linter (7 cases) and optionally checks a
report passed as an argument; it runs in `run_tests.sh`. Pointed at the pre-fix
published report it correctly reports `filterGenomes()`.

Only the undefined direction is checked. "Defined but never called" was tried and
dropped as too noisy — `searchNorm` and `searchMatches` are invoked from other JS
rather than from markup, and `filterKCBHits` is legitimately uncalled when the KCB tab
has no hits to render a search box for.

### Resuming a Run After Editing Scripts

Two traps, both hit on 2026-08-25.

**Editing a file under `scripts/` now invalidates the cache correctly** (fixed
2026-08-25; it did not before). Modules invoke scripts as
`python ${projectDir}/scripts/foo.py` — an interpolated path, not a declared `path`
input — so Nextflow's task hash never saw them and `-resume` happily reused output
built from code that had since changed. A *silent* wrong answer, the worst kind.

Each script-running process now embeds a digest of the scripts it depends on:

```groovy
# scripts-version: ${Utils.scriptsHash(projectDir, ['visualize_results.py', 'utils', 'viz'])}
```

The script block's text is part of the task hash, so a changed digest re-runs the
task. Dependencies are listed **per process**, not hashed as one tree: a whole-tree
digest would make `RENAME_GENOMES` depend on plotting code, and since it feeds
antiSMASH, editing a chart would invalidate 1,735 antiSMASH tasks. Measured: editing
`viz/rarefaction.py` changes the `VISUALIZE_RESULTS` digest and leaves
`RENAME_GENOMES` and `GCF_BIOSYNTHETIC_TREE` untouched.

Note `path` inputs were tried first and rejected. A directory `path` input does **not**
hash its contents — a process staging `scripts/` served stale output while reporting
`cached=1` — and staging individual files breaks the scripts'
`sys.path.insert(0, Path(__file__).parent)` imports.

**When adding a process that runs a script**, add the marker and list its
dependencies. `tests/check_script_deps.py` (in `run_tests.sh`) fails if a declared
list stops covering a script's real imports.

**`-resume <run-name>` can silently fall back to the wrong session.** The task hash
begins with the session UUID, so resuming the wrong session misses every entry and
the pipeline starts from scratch — including the NCBI download. A name that fails to
resolve does not error; it quietly resumes the most recent session, which is easily a
3-second `-preview`. Resume by **UUID**, taken from `.nextflow/history` (column 6):

```bash
awk -F'\t' '{print $3, $6}' .nextflow/history   # run name -> session UUID
nextflow run main.nf -resume <uuid> --taxon "Pantoea"
```

Confirm it bound before letting it run: `-dump-hashes` prints the session UUID as the
first hash entry, and the summary line should report a large `cached=` count. If you
see `cached=0` and `NCBI_DATASETS_DOWNLOAD` starting, kill it — the resume missed.

### BiG-SCAPE Partitioning (`--bigscape_partition`)

**Off by default.** Validated at 185-518 BGCs, not at the scale that needs it, and never
run on SLURM. Turn it on deliberately.

```
PARTITION_BGCS -> BIGSCAPE_PARTITION (one per partition) -> MERGE_BIGSCAPE -> CLUSTERING_STATS
```

| Parameter | Default | Notes |
|-----------|---------|-------|
| `bigscape_partition` | `false` | Enables the partitioned path |
| `bigscape_partition_identity` | `0.60` | 0.60-0.80 is the verified safe window; 0.90 splits real families |
| `bigscape_partition_threshold` | `10000` | **Total BGCs in the run** — not base pairs, not a partition size — below which everything goes in one partition |
| `bigscape_partition_max_size` | `0` | Largest partition, in BGCs. `0` derives it from the task's memory allocation. A **safety valve, not a tuning knob** |
| `bigscape_partition_memory_margin` | `0.85` | Fraction of the allocation to budget when deriving the cap |

**The threshold is a dataset-level switch, the cap is per-partition.** They are not a
matched pair, which an earlier `_min`/`_max` naming wrongly implied.

**Partitioning is a net loss at small scale**, which is what the threshold is for: every
partition re-pays BiG-SCAPE's fixed Pfam-load cost. Measured on identical inputs, 185 BGCs
took **30 s as one job against 112 s across 19 partitions** (3.7x slower), and 518 BGCs
took 93 s against 191 s (2.1x). Modelling the fixed cost against the quadratic term puts
the crossover near **10,000-12,000 BGCs**, hence the default. Below it, partitioning costs
time and buys nothing, because memory is not yet a constraint either.

**The cap should follow the memory allocation, not a guess**, so it is derived by
inverting the measured fit `GB = 1.14 + 1.29e-7*n^2` against `task.memory`:

| RAM | largest partition |
|----:|------------------:|
| 16 GB | 10,700 BGCs |
| 32 GB | 15,500 |
| 48 GB | 19,100 |
| 64 GB | 22,100 |
| 128 GB | 31,400 |
| 256 GB | 44,400 |

Those are at 85% of the allocation, not 100%, because **the fit is being extrapolated well
past its data**: it was measured to 10,000 BGCs, and 128 GB implies ~31,000 — a 3.1x reach,
4.4x at 256 GB. The margin costs ~8% of the cap and buys ~19 GB of headroom at 128 GB,
against an OOM kill that discards hours of clustering.

**To run smaller partitions, lower the memory allocation**, which moves the cap and the
request together. `bigscape_partition_max_size` is a safety valve for when the derivation
is wrong for your data — BGCs with unusually rich domain content could need more memory
per BGC than the fit predicts, and that is the one case the allocation lever cannot
express, since you want a smaller cap at the *same* memory.

A component above the cap is chunked, which **can** split a real family, so the partitioner
warns when it fires.

**`PARTITION_BGCS` runs before BiG-SCAPE**, so it cannot read pepM from a clustering
database — it extracts CDS translations from the region GenBanks and finds pepM by
`hmmsearch` against PF13714. Cross-checked against the database route: both give 6
partitions, largest 236, work 57% on Erwiniaceae.

**Single linkage, deliberately.** Two BGCs must share a partition whenever they could
possibly cluster. Single linkage merges on *any* qualifying link, so it never separates a
pair a stricter criterion would have joined.

**`MERGE_BIGSCAPE` remaps every id.** Partition databases carry colliding autoincrements,
so each table is copied with its primary key offset and every foreign key rewritten;
`SPEC` in the script is a topological sort. `run` and `edge_params` are deduplicated
rather than offset — every partition ran identical parameters, and offsetting would invent
parameter sets that never existed.

**The merged `distance` table is deliberately incomplete.** It holds only within-partition
comparisons — that is the point. Readers must treat a missing pair as *not compared*, not
as distance zero.

**Two guards, both earned in testing:**

- *Completeness.* A missing partition database merges cleanly, passes referential
  integrity, and is simply short some BGCs. Caught when 14 of 19 partition databases gave
  179 regions against the reference's 185. `--partitions` now derives the expected count
  and fails on a mismatch.
- *Filename uniqueness.* Partitions stage region GenBanks flat, so two regions sharing a
  filename would silently overwrite. Region files are named by contig accession and are
  unique in practice; the partitioner fails loudly rather than lose BGCs if that stops
  holding.

**`CLUSTERING_STATS` replaced `EXTRACT_CLUSTERING_STATS`.** The old process parsed
BiG-SCAPE's `output_files/*_clustering_c*.tsv`, which ties it to one output *directory*;
partitioned runs have one per partition. `stats_from_db.py` reads the database instead and
was verified to emit **byte-identical** JSON on the unpartitioned Erwiniaceae run, so both
paths use it and the TSV-parsing module is gone.

`PEPM_ALL_BY_ALL` output now reaches the report: the pepM-vs-similarity figure and
correlation table land in **GCF Analysis** (biological evidence, beside the clustering it
justifies), the partitioning feasibility table in **Pipeline Info** (operational). Both
sections return '' when the analysis did not run, so a report without it still renders.

**Verified end to end (2026-09-04).** A partitioned Erwiniaceae run reproduces the
unpartitioned clustering exactly: 333 BGCs, 19 families, 8 singletons, largest family 215,
and **23,995 co-membership pairs with 0 split and 0 merged**. Six partitions, largest
236; 31,566 within-partition comparisons against 55,278 all-pairs.

Four bugs surfaced only under Nextflow, all invisible to standalone testing:

1. **`Path.rglob` does not descend into symlinked directories.** Nextflow stages antiSMASH
   results as symlinks; rglob found 0 GBKs where `glob.glob(recursive=True)` found 333.
   Python 3.13 added `recurse_symlinks` but defaults it to False. BiG-SCAPE's own loader
   carries the same workaround.
2. **The partition manifest held paths relative to the partitioner's work directory**,
   which does not exist downstream — Nextflow resolved them against the launch dir into
   broken symlinks. Paths are now `Path(g).resolve()`.
3. **A stale cached `PARTITION_BGCS` task served the old manifest** after the fix. Removing
   the task directory forced the re-run.
4. **`BIGSCAPE_PARTITION` had no conda environment.** `withName` matches exactly, so the
   existing `'BIGSCAPE'` selector never covered it, and unlike a *stale* selector Nextflow
   does not warn about a missing one — it fails at runtime with `command not found`.

**Trees under partitioning: one global centre tree.**
A partitioned run's merged `distance` table holds only within-partition comparisons, so a
global all-BGCs tree substitutes a constant for every cross-partition pair — 23,712 of
55,278 cells on Erwiniaceae. The clustering is unaffected (BiG-SCAPE never compared those
pairs either) but a tree built on a uniform constant has an arbitrary backbone.

- **`BIGSCAPE_CENTERS`** re-runs BiG-SCAPE over one representative GBK per family, so
  every centre pair is *measured*. On Erwiniaceae that turned **92 of 171 substituted
  centre pairs into 0**, in 16 seconds over 19 centres; a later full run measured
  **171 of 171**. `bgc_gcf_tree.py --centers_db` consumes it. This scales because centre
  count tracks diversity rather than BGC count (19 Erwiniaceae, 81 Streptomyces, 100
  combined).

**`PARTITION_TREES` was built and then removed.** It drew one all-BGCs tree per
partition, on the theory that the centre tree gives the global view and these give the
detail inside it. The composition table killed that argument: a partition is a pepM
identity component sized to bound BiG-SCAPE's memory, not a biological unit. On
Erwiniaceae partition 0 held 236 BGCs across **2 families** — 215 in one — so 91% of the
figure was within-family variation at a leaf count nobody can read, while partition 2 was
4 BGCs in 1 family. The trees also appeared only on partitioned runs, reintroducing the
mode-dependence that removing the all-BGCs tree had just eliminated. If per-BGC detail is
wanted, the unit to draw is a **family**, not a partition.

`BIGSCAPE_CENTERS` runs only when `--bigscape_partition` is on; the unpartitioned path is
untouched. Partitions with fewer than three BGCs are skipped rather than failed, since singleton
partitions are normal.

The all-BGCs tree does not scale regardless of partitioning: 121,000 leaves is 7.3e9 pairs
and 1,890 GB, and is not a readable figure at any resolution. The centre tree is the
global view that survives; at extreme diversity even it needs its own partitioning
(52,000 centres would want 358 GB), but that ceiling is 5-270x further out.

`GCF_BIOSYNTHETIC_TREE` failed to exit in the first partitioned run. Both tree scripts were
subsequently run by hand against the merged database and completed normally (exit 0), so
that was environmental rather than a hang in the tree code.

Verified earlier: the merge preserves clustering exactly (0 split, 0 merged co-membership
against the reference on the 181 regions compared), BiG-SCAPE runs on a single-BGC
partition, and the DAG resolves. **Not yet verified: a full pipeline run on the
partitioned path**, which is the remaining gap before trusting it.

### Region size: `antismash_phosphonate_neighbourhood`

antiSMASH sizes a region as the rule core plus a fixed neighbourhood. The strict
`phosphonate` rule declares `NEIGHBOURHOOD 5`, and that is too small for this chemistry:
on *P. ananatis* LMG 5342 the region ends at 810,246 while the HiVir cluster runs to
814,435, so **the MFS transporter, the hypothetical, the FMN reductase and the second
ATP-grasp are missing from all 215 regions of that family** — every gene-content metric
in the report sees 8 of 12 genes. Measured on that genome:

| NEIGHBOURHOOD | region | size | of the 12,526 bp cluster | flank added |
|---:|---|---:|---|---|
| 5 (antiSMASH's own) | 796,909-810,246 | 13,338 bp | 8,337 bp | 5.0 kb up, 0 down |
| **10 (default here)** | 791,909-815,246 | 23,338 bp | **all of it** | 10.0 kb up, 0.8 kb down |
| 20 | 781,909-825,246 | 43,338 bp | all of it | 20.0 kb up, 10.8 kb down |

**`--hmmdetection-strictness relaxed` is not an alternative.** `phosphonate-like`, the
relaxed rule carrying `NEIGHBOURHOOD 20`, never fires on this cluster — verified on the
same genome, with no trace of it anywhere in the output. Relaxed produced a byte-identical
region to strict. The strict rule matches, so its neighbourhood is what applies.

antiSMASH exposes no option for a bacterial neighbourhood (only fungal multipliers) and
none for pointing at a different rule file, so `scripts/genome/patch_antismash_neighbourhood.py`
edits the installed `strict.txt` before each antiSMASH invocation. It is **idempotent**
(a file already at the value is untouched) and **atomic** (temp file plus `os.replace`),
because the conda environment is shared by every concurrent ANTISMASH task and is
recreated whenever the spec changes or the cache is cleared — which is why it runs per
task rather than once.

**Caveats worth holding onto.** The neighbourhood is symmetric, so 10 also pulls in 10 kb
upstream that the cluster does not contain; 10 is tuned to HiVir, not a general truth.
The value is part of `Utils.antismashParamsHash` **only when set**, so a `null` run
reproduces the hash results were already produced under (verified: `55968d66...`, which
matches the `.antismash_meta` of the Erwiniaceae run) and keeps them reusable, while any
other value invalidates them by design — the regions genuinely differ. Changing it also
moves every BiG-SCAPE distance, since complete regions are compared end to end.

### `BIGSCAPE_REFERENCES` — distance to the characterised clusters

Runs BiG-SCAPE a second time over a **copy** of the finished clustering database, with
`--reference-dir` pointing at `assets/phosphonate_reference_bgcs`. BiG-SCAPE recognises
the work already in the database and computes only the pairs involving a reference; the
published clustering is never touched.

**Why not simply pass `--reference-dir` to the main run**, which is what the code did
until this process existed: a reference inside the GCF cutoff *joins a family*, and the
~15 downstream scripts read the database with no way to tell a reference from the
dataset. Measured on a 14-genome subset with one reference loaded:

| | without | with |
|---|---:|---:|
| `total_bgcs` reported by `stats_from_db.py` | 18 | **19** |
| members in the pantaphos family | 10 | **11** |
| genomes in `genome_gcf_mapping` | 14 | **15** — the extra one named `phosphonate_reference_bgcs` |

Query family membership itself was identical in every comparison, so the damage is to
the counts and labels rather than to the clustering.

**Cost** (8 cores, same flags as the main run):

| BGCs | main clustering | reference pass | share |
|---:|---:|---:|---:|
| 334 (Erwiniaceae) | 56.7 s | **9.4 s** | 17% |
| 2,000 (replicated) | 388.9 s | 40.3 s | 10% |
| 4,000 (replicated) | 1,049.7 s | 134.4 s | 13% |

`--db-only-output` is what makes it cheap: without it the pass spends 42 s of 78 s at
2,000 BGCs regenerating HTML and trees nothing here reads. Verified to leave all 10,010
and 20,010 reference distances identical. Copying the database costs 0.82 s at 856 MB.

**Not wired on the partitioned path.** A merged partition database holds only
within-partition distances, so this pass would compute every cross-partition pair that
partitioning exists to avoid. The subworkflow warns and skips; measuring references there
means running the pass per partition.

**References are matched by content hash, not path.** BiG-SCAPE deduplicates input on the
sha256 of the file — read as *text*, so a CRLF file does not hash as its raw bytes — and
when a reference is byte-identical to a query BGC it drops the reference and keeps the
query, logged at INFO. `reference_distances.py` hashes the same way, so that case is
reported as an exact match instead of vanishing. It also warns when a reference never
loaded at all, which is this directory's characteristic failure: a file without an
antiSMASH region feature is ignored silently.

**References carry `contig_edge=True`, deliberately.** A curated reference is the cluster;
an antiSMASH query region is a rule core plus a neighbourhood, so the two are delimited
differently and comparing them end to end scores the difference in *boundaries* as a
difference in biology. The curated pantaphos sat 0.364 from LMG 5342's **own** region for
exactly that reason. Under `auto`, BiG-SCAPE compares a pair by its shared part whenever
either record is on a contig edge, and the edge parameter is stored per run — so setting
the flag on the reference costs nothing and recomputes nothing. Measured on a 20-BGC set:
BGCs within 0.30 of pantaphos went from 3 to 11, LMG 5342's own region from 0.364 to
0.000, and genuinely different clusters did not move (0.946 and 0.953 either way). The
flag is not a claim that the record runs off a contig, and the `note` qualifier
`make_reference_bgc.py` writes into every reference says so. See
`assets/phosphonate_reference_bgcs/README.md`.

### BiG-SCAPE distances are not reproducible without `PYTHONHASHSEED=0`

BiG-SCAPE 2.0.1 returns different distances for identical input. Measured 2026-09-16 on
Erwiniaceae subsets: **8-14% of pairs differed between repeat runs**, by up to 0.18 (0.27
for a reference pair).

`file_input/load_files.py` dedupes GBKs through a Python `set` keyed on the sha256
*string*, so load order follows the per-process hash salt. That order sets record ids,
which set which record is A in each pair, and `comparison/workflow.py` extends pairs
asymmetrically. In all 226 differing comparisons the A/B orientation had flipped; no
distance ever differed without a flip.

`export PYTHONHASHSEED=0` — now set in `BIGSCAPE`, `BIGSCAPE_PARTITION`,
`BIGSCAPE_CENTERS` and `BIGSCAPE_REFERENCES` — made two runs identical across all 190
pairs and all family centres.

What it does and does not reach: differences appeared only at distances >= 0.66, and
across 6,319 pair comparisons nothing below 0.5 ever moved, so **GCF membership at the
0.30 cutoff was identical every time**. Family *centres* did move, 2 of 5 between
identical runs, so GCF representatives were unstable before this. Even pinned, changing
the genome set permutes load order again — like GCF ids, far distances are comparable
only within one run.

### `NOVELTY_SCORE` — ranking families for laboratory follow-up

`priority = distance x evidence`. They multiply because both are necessary: a maximally
divergent single truncated region is not a lead.

**The distance axis was reference-biased, and the bias was the whole signal.** The first
version measured identity to the 7 characterised references. Six of those are
*Streptomyces*; the one Enterobacterial reference (HvrC) is the pantaphos synthase.
Measured over all 333 Erwiniaceae regions the identities are bimodal with an empty gap:

```
22-45%   every Decarboxylase, Reductase and Transaminase family   (n=364)
         nothing at all between 45.4% and 93.8%
94-100%  every Synthase family                                    (n=944)
```

That is a binary readout of *"does a same-taxon reference exist"*. Consequences:

- Ranks 1, 3, 4, 6 and 14 were **all Reductase** — top-ranked only because VlpB
  (*S. durhamensis*) is the most distant reference in the set.
- Within a class the value was effectively constant: every Reductase member scored
  0.75-0.78 regardless of its own sequence. Ranks 5-13 were separated by noise.
- Two families could not be ranked **at all**, having no coupling class to score.

**Isolation replaces it.** BiG-SCAPE already writes the complete all-pairs distance
matrix — 55,278 rows for 333 regions, exactly `n(n-1)/2`, no reference set involved.
Isolation is the **median** over members of that member's smallest distance to any
region *outside* its family. Median, so one atypical member cannot make a family look
either connected or isolated.

| check | result |
|---|---|
| Spearman rho(isolation, reference divergence) | **-0.17** — carries independent information |
| Spearman rho(isolation, family size) | -0.45, but singletons mean 0.490 vs multi-member 0.505 and singletons span 0.256-0.861 — **not a size artefact** |
| families rankable | **19 / 19** (was 17; GCF-16, the single most isolated family in the run at 0.861, was one of the two that could not be ranked) |

**Reference identity is kept, demoted to two jobs.** It is published as context —
`reference_status`, `reference_pct_id`, `reference_organism` — because seeing
"23.9% to *Streptomyces durhamensis*" is what tells a reader the number describes the
reference set rather than the family. And it gates: a family at >= `CHARACTERISED_PCT`
(60%) has its distance zeroed, because it is a solved cluster whatever its isolation.
**The threshold sits inside the empty 45.4-93.8% gap, so any value in that range gives
identical results** — robust, not tuned.

**Effect on the Erwiniaceae ranking:**

| GCF | old rank | new rank | note |
|---|---:|---:|---|
| 4 | 7 | **1** | isolation 0.814 |
| 11 | 12 | **3** | *Winslowiella iniecta* B149 — the cluster the lab independently chose to characterise |
| 16 | unranked | **7** | most isolated family in the run |
| 1 | **1** | 9 | Reductase; was top on VlpB distance alone |
| 9 | 3 | 11 | Reductase |
| 2 | 17 | 18 | pantaphos — `characterised`, distance zeroed, correctly at the bottom |

GCF-11 rising from 12th to 3rd is worth noting: the lab picked that cluster for
characterisation on independent grounds, and the unbiased axis agrees with them where
the reference-based one did not.

### `utils/domain_functions.py` — Pfam accession → biosynthetic role

**Why accessions and not product text.** A keyword metric over NCBI product names was
tried first and inverted on the two lab-confirmed clusters. The cause was not subtle:
`serine hydroxymethyltransferase` matched `methyltransferase` and
`aspartate-semialdehyde dehydrogenase` matched `dehydrogenase`, so two core
amino-acid metabolism genes counted as tailoring chemistry. Pfam separates them by
construction — SHMT is **PF00464**, a methyltransferase is **PF13649**.

Coverage is the second reason. antiSMASH scans every region with clusterhmmer whether
or not NCBI annotated the assembly, so domains reach **77.6% of CDS uniformly**, where
product text reached 40.2% and tracked annotation quality.

| category | what it means |
|---|---|
| `core` | the phosphonate pathway itself — pepM and the coupling enzymes |
| `tailoring` | chemistry past the coupling step |
| `lipid` | lipid handling, kept separate given the open lipid/small-molecule question |
| `transport` | moves the product |
| `regulation` | controls expression |
| `mobile` | how the cluster arrived, never counted as chemistry |
| `primary` | core metabolism swept in at region edges — **explicitly excluded from tailoring** |

`PRIMARY` exists because antiSMASH region boundaries catch chromosomal neighbours.
Peptidyl-tRNA hydrolase and DnaB are not tailoring enzymes however many sit beside a
BGC, and counting them is exactly what broke the first attempt.

`elaboration()` counts **distinct** tailoring/lipid domains, not occurrences: three
copies of one methyltransferase domain is one kind of chemistry, and counting
occurrences would let a tandem duplication look like elaboration.

**The map is curated and incomplete** — ~100 domains covering 93% of observed hits.
Everything else returns `other`, which means *not classified here*, never *not a
biosynthetic gene*. Do not read an absence as evidence.

**Effect on the three characterised clusters** (elaboration, domain-based):

| family | truth | keyword version | domain version |
|---|---|---:|---:|
| GCF-18 | confirmed **phosphonolipid** | 3.00 | **2.00** |
| GCF-11 | confirmed **not** a lipid | 1.00 | **4.00** |
| GCF-2 | pantaphos, small molecule | — | **5.33** |

The ordering is now lipid < non-lipid < small molecule, where the keyword version had
it inverted. GCF-18 also shows `primary` = 7.00, the highest of any family, confirming
that region is mostly swept-in housekeeping around a minimal cluster.

**This is not validation of the lipid hypothesis.** There is still exactly one
confirmed lipid, two metrics have now given two different answers, and GCF-14 and
GCF-16 tie GCF-18 at 2.00. What changed is that the metric is no longer confounded by
annotation quality or by keyword collisions — it is worth computing, not yet worth
predicting from.

### Consensus gene content in the report

`build_consensus_clusters_section()` renders one consensus cluster per family from
`gcf_consensus_clusters.tsv`, replacing the single-representative view. A representative
shows one genome's annotation, which on this data is usually a bad draw.

**Prevalence is the column that matters.** A gene at 1.00 is in every member and helps
define the family; one at 0.24 is accessory and may be a neighbour the region boundary
caught. The GCF-2 consensus shows the pantaphos cluster as eight genes at prevalence
0.98-1.00 — including a GNAT acetyltransferase and an ATP-grasp protein that
*P. ananatis* LMG 5342's own annotation calls "hypothetical" — with the SpoIIE
phosphatase correctly at 0.69 and tagged `primary`.

The Naming column carries provenance: a name backed by one genome is flagged as one
genome's opinion, not consensus.

### `GCF_ANNOTATION_TRANSFER` — complete a BGC's annotation from its relatives

**The problem it solves.** Only **40.2%** of CDS in the Erwiniaceae run carry an
informative product, and **170 of 333 regions carry none at all**. Those assemblies
are GenBank-only, with no functional annotation — the clusters are not bare, they
are unreadable. Any gene-content metric built on product text is therefore comparing
NCBI annotation pipelines rather than biology, and it fails *silently*: an
unannotated cluster scores zero on every functional category and looks minimal.

This was found the hard way. An "elaboration" metric (tailoring + transport genes
per cluster) appeared to separate the confirmed phosphonolipid from the confirmed
non-lipid, until GCF-8 scored a perfect zero — not a bare cluster, a genome with
0/22 products annotated.

**Why transfer works here.** GCF members are homologous by construction, and
annotation quality across them is extremely uneven: the typical family has a
*median* member at 0% and a *best* member at 88-100%. One RefSeq-quality genome
carries the whole family.

| | before | after |
|---|---:|---:|
| CDS with an informative product | **40.2%** | **80.2%** |

Per family, the ones that were unusable become usable: GCF-2 29.7 → 80.6, GCF-5
29.9 → 84.7, GCF-6 29.9 → 81.4, GCF-9 35.1 → 87.7, GCF-7 14.3 → 61.9.

**How.** `diamond blastp` all-vs-all *within* each family (never across), orthologue
groups by single linkage over edges passing `--min_identity` (50%) and
**mutual** `--min_coverage` (0.70, both directions — one-sided coverage would accept
a short fragment aligning inside a long multidomain protein). Consensus product per
group by majority vote, ties broken toward the longer (more specific) name.

**Validation.**

| control | result |
|---|---|
| pepM/Ppd groups get the correct consensus | every family with an annotated one |
| GCF-2 resolves into coherent groups | 8 groups at prevalence 0.98-1.00 across 215 members |
| transferred calls resting on a **single** source genome | 127 / 2,074 (**6.1%**) |
| groups where annotated members disagree | 73 / 313 (23.3%), mostly RefSeq-vs-GenBank synonyms |

The consensus GCF-2 cluster recovers genes LMG 5342's own annotation calls
"hypothetical" — a GNAT N-acetyltransferase and an ATP-grasp protein among them.

**What it is not: an observation.** A transferred product is an inference from a
homologue. Every row in `gcf_annotation_transfer.tsv` carries `origin`
(observed / transferred / none) plus the source genomes, `n_sources`, `n_agree` and
`n_disagree`. The failure mode is error propagation — one mis-annotated RefSeq gene
becomes N of them, and the agreement count looks reassuring *because they share a
single origin*. **`n_sources` is the column that exposes that**; treat a call with
`n_sources=1` as one genome's opinion, not as consensus.

**Three families gain nothing** (GCF-8, 13, 16 on Erwiniaceae): singletons with no
annotated relative. They are emitted with `origin=none` so downstream code can
exclude them rather than read absent annotation as absent genes.

**It does not resolve phosphonolipid vs small molecule.** The two lab-confirmed
clusters are the worst possible pair for it: GCF-18 is a well-annotated singleton
(nothing to transfer from) and GCF-11's two members are uniformly half-annotated
(both missing the same genes — its `NTP_transf_3` CDS is labelled "hypothetical
protein"). Neither improves. See "Predicting phosphonolipid vs. small molecule is
unsolved" above.

### `PEPM_ALL_BY_ALL` — pepM identity vs gene-neighbourhood similarity

Reproduces Yu et al. (PNAS 2013;110(51):20759) Fig. 2B on the run's own data, and answers
whether pepM identity could partition BiG-SCAPE's all-pairs problem.

This is **not** the pepM-vs-references divergence plot removed on 2026-08-27 (below).
That compared each pepM against a reference set, which the paper does not license. This
is the paper's actual analysis: all pepMs compared pairwise against each other.

**It is mostly a join.** BiG-SCAPE's `distance` table already holds the y-axis for every
pair — `jaccard` is shared domain content, the analogue of the paper's "fraction of
homologous genes shared". Only the pepM axis is new.

Method follows the paper: identity from one alignment with **pairwise deletion of missing
sites**, not BLAST. `hmmalign` against PF13714 gives that and is linear in sequence count
where all-by-all alignment is quadratic. Only match columns count, so fusion proteins are
not penalised for residues nobody was aligned against.

**Measured on Erwiniaceae (333 BGCs, 55,278 pairs, 100% pepM coverage):**

| | r | r² | slope |
|---|---|---|---|
| vs `jaccard` | +0.598 | 0.358 | +2.49 |
| vs BiG-SCAPE similarity | +0.641 | 0.411 | +2.71 |

The correlation above 60% identity is confirmed, but in this data the relationship is a
**step, not a line**: median neighbourhood similarity is ~0.03 below 0.88 identity and
jumps to 0.998 at ≥0.98. The paper's dataset spanned all known producers and had a
populated middle; one family does not.

**Every same-GCF pair has pepM identity ≥ 0.901** (median 1.000), so a cut anywhere from
0.50 to 0.90 is lossless here.

**But it does not partition well enough on its own.** Single-linkage at any threshold
leaves one component holding 71% of BGCs, because the dominant GCF is genuinely one
cluster of near-identical pepMs. Work falls to ~52% of a single job — a 2x saving, not
the 10-100x needed:

| cut | components | largest | same-GCF lost | work vs one job |
|----:|-----------:|--------:|--------------:|----------------:|
| 0.60 | 6 | 236 (71%) | 0 | 57% |
| 0.90 | 13 | 236 (71%) | 0 | 52% |

**Streptomyces settles it: partitioning works, and Erwiniaceae could not show it.**
1,573 complete genomes, 1,094 organism names, 185 BGCs in 81 families with a largest
family of 15 — against Erwiniaceae's 333 BGCs in 19 families with a largest of 215.

| | Erwiniaceae | Streptomyces |
|---|---:|---:|
| r vs BiG-SCAPE similarity | +0.641 | **+0.884** |
| r² | 0.411 | **0.781** |
| components @0.60 | 6 | **19** |
| largest component @0.60 | 236 (71%) | **38 (21%)** |
| same-GCF pairs lost @0.60 | 0 | **0** |
| work vs one job | 57% | **11%** |

The correlation is far stronger because Streptomyces has the populated middle range the
paper's kingdom-wide dataset had and one family does not. Cuts stay lossless to 0.80
(39 components, largest 25, work 6%); 0.90 separates same-GCF pairs, so **0.60-0.80 is
the safe window**.

**Verified by re-running, not just predicted** (`bench_bigscape_partitioned.py`). At 0.60
both taxa rebuild the reference clustering exactly — 81/81 families and 466/466
co-membership pairs for Streptomyces, 19/19 and 23,995/23,995 for Erwiniaceae, zero split,
zero merged, **ARI 1.0000**. Erwiniaceae is the stronger evidence despite partitioning
badly: its largest partition holds 71% of BGCs, so the agreement is not an artefact of
small partitions.

**A mixed set behaves as the middle case, and dilution is the mechanism.** Combining both
taxa (518 BGCs) gives 100 families — exactly 19 + 81, so *no cross-taxon families form*:
Erwiniaceae and Streptomyces phosphonate BGCs are disjoint lineages, not variations on
shared clusters.

| | Erwiniaceae | Combined | Streptomyces |
|---|---:|---:|---:|
| BGCs | 333 | 518 | 185 |
| families | 19 | 100 | 81 |
| r vs BiG-SCAPE similarity | +0.641 | +0.683 | +0.884 |
| components @0.60 | 6 | 24 | 19 |
| largest @0.60 | 236 (71%) | **236 (46%)** | 38 (21%) |
| work vs one job | 57% | **26%** | 11% |
| ARI vs unpartitioned | 1.0000 | **1.0000** | 1.0000 |

The *Pantoea* component stays **236 BGCs, unchanged** — only its share falls, because the
denominator grew. Adding diversity does not break dense clusters apart, it dilutes them.

**Measured, not inferred** (`bench_component_rarefaction.py --add-from`). Holding
Erwiniaceae fixed and adding Streptomyces in nine steps, the largest component is 236 at
every single step:

| Streptomyces added | total | largest | share | components |
|---:|---:|---:|---:|---:|
| 0 | 333 | 236 | 70.9% | 6.0 |
| 92 | 425 | **236** | 55.5% | 21.7 |
| 185 | 518 | **236** | 45.6% | 24.0 |

Growth across the whole addition: **+0 BGCs**. New diversity adds components beside the
dense cluster; it never enlarges it.

So **peak memory is set by the most deeply sampled single clade, not by dataset size.**
A million genomes spread across many clades does not enlarge any one component. A million
concentrated on one over-sequenced clade would — and NCBI *is* skewed, so a heavily
sequenced phosphonate-carrying clade is the case to watch.

Note the plain random-subsample rarefaction in the same script is **uninformative for
this question and says so**: drawing from a pool where 45.6% of BGCs sit in one component
returns ~45.6% at every depth (slope ratio 0.95, share pinned 44-47% over a 12x range).
That is arithmetic, not biology. `--add-from` is the mode that bears on scaling.

At 0.90 it breaks as predicted — 65 partitions, 81 -> 83 families, **23 same-family pairs
split**, ARI 0.974. That confirms the cheap check is directionally sound, but note it
predicted **7**, not 23: `lost` uses pairwise similarity >= 0.70 as a proxy for family
membership, while BiG-SCAPE families are transitively closed and 9% of same-family pairs
sit below that cutoff, joined through a third BGC. **Treat a non-zero `lost` as
disqualifying rather than as a damage budget**; zero remains a reliable all-clear.

Applying the measured memory fit `GB = 1.14 + 1.29e-7*n^2` at 121,000 BGCs:

| | peak RAM |
|---|---:|
| one job, no partition | 1,890 GB |
| Erwiniaceae-like split (71%) | 953 GB |
| **Streptomyces-like split (21%)** | **84 GB** |
| Streptomyces at 0.80 (14%) | 38 GB |

So a million genomes goes from *needing a 2 TB node that may not exist* to fitting the
existing 128 GB `process_high_memory` allocation. **Judge partitioning on diversity, not
BGC count** — the earlier "2x is not enough" verdict was measured on the one clade where
it could not work.

**Pantoea vs Erwiniaceae confirms the scope dependence, by being nearly degenerate.**
`--organism Pantoea` slices the same run without re-clustering (the `distance` values are
pairwise and scope-independent; only family labels come from the parent run). Pantoea is
**319 of Erwiniaceae's 333 BGCs, 95.8%** — so the two are almost the same data, and every
correlation moves by ~0.01:

| | Pantoea (319) | Erwiniaceae (333) |
|---|---:|---:|
| r vs `jaccard` | +0.588 | +0.598 |
| r vs BiG-SCAPE similarity | +0.631 | +0.641 |
| largest component @0.60 | 235 (74%) | 236 (71%) |
| work vs one job @0.60 | 61% | 57% |
| **components @0.60** | **3** | **6** |

The last row is the signal. Adding 14 BGCs (+4%) from four other genera — *Erwinia* 9,
*Mixta* 2, *Winslowiella* 2 — **doubled the component count**. Taxonomic breadth, not BGC
count, is what partitions pepM space. Judge partitioning on a genuinely diverse taxon;
neither of these is one.

`PF13714` is resolved to its versioned accession (`PF13714.13`) by scanning the HMM file:
hmmfetch's index keys on the exact string, and the version moves when `pfam_release` is
bumped.

### Why There Is No pepM Tree or pepM-Divergence Analysis

An earlier version built pepM and per-class coupling enzyme trees (`COUPLING_ENZYME_TREE`)
and a pepM-vs-coupling divergence plot. Both were removed on 2026-08-27.

Yu et al. (PNAS 2013;110(51):20759, doi:10.1073/pnas.1315107110) established that pepM
identity predicts gene-neighbourhood similarity above ~60% identity — but that is a
statement about pepM sequences compared **pairwise against each other**, not about a
pepM compared against a reference set. The divergence plot applied it the second way,
which the paper does not license.

More importantly, the correlation exists because pepM is a *proxy* for the neighbourhood.
BiG-SCAPE measures the neighbourhood directly, over full domain content and adjacency,
and the `distance` table already holds every pairwise comparison. Reproducing the
pairwise pepM analysis would rebuild a proxy for something already measured directly.

Classification now rests on two non-overlapping signals:

| Signal | Answers | Source |
|--------|---------|--------|
| GCF membership | which BGCs group together | BiG-SCAPE, whole neighbourhood |
| Coupling class + support | what chemistry, and how well evidenced | SMCOG markers + reference identity |

`scripts/bgc_coupling_tree.py` is retained as a standalone post-pipeline tool. It is no
longer a pipeline stage, but phylogenetic placement remains the right way to adjudicate
an individual ambiguous enzyme (the PalB/SMCOG1019 case), which a scalar identity cannot
do. Running it needs `hmmer` and `fasttree`, which are no longer in any pipeline conda
environment — install them separately or pass paths with `--hmmbuild` etc.


### Genome Name Conventions

Two different spellings of a genome name coexist, and they do **not** compare equal:

| Source | Form | Example |
|--------|------|---------|
| `region_counts.tsv` `record` column | **with** file extension | `Pantoea_ananatis_01-1.gbff` |
| `renamed_genomes/` filenames | **with** extension | `Pantoea_ananatis_01-1.gbff` |
| antiSMASH result directory | **without** | `Pantoea_ananatis_01-1` |
| BiG-SCAPE `gbk.path` component | **without** | `Pantoea_ananatis_01-1` |

`PHYLOGENY` joins `row.record` against `genome.name` — both extension-bearing, so that
join is correct. But anything correlating the BiG-SCAPE DB with `region_counts.tsv` must
normalise first: `viz/rarefaction.py` uses `strip_genome_suffix()`, which removes only
known genbank suffixes (`.gbff`, `.gbk`, `.gb`, `.genbank`). Do **not** use `Path.stem` —
it eats the version off accession-style names (`..._GCA_963520565.1` → `..._GCA_963520565`).

An unnormalised comparison does not fail loudly; it silently yields an empty intersection
and, if you are merging sets, doubles your genome count. `generate_rarefaction_curve`
guards against this by refusing to pad when the two name sets share nothing.

### Data Key Conventions

BGC regions are uniquely identified using `region_name` (e.g., "40.1" = record_index 40, region 1):

- **KCB lookup**: `(genome, region_name)` → KnownClusterBlast hit info
- **BGC-to-GCF mapping**: `(genome, region_name)` → GCF family assignment
- **JSON serialization**: `"genome|region_name"` format (e.g., `"Streptomyces_coelicolor_A32|40.1"`)

This avoids key collisions since `region` numbers are only unique within a record/contig, not within a genome. The `region_name` matches antiSMASH's naming convention directly.

Key files:
- `scripts/clustering/extract_gcf_representatives.py`: `load_kcb_lookup()`, `build_record_index_map()`, `extract_genome_gcf_mapping()`
- `scripts/analysis/tabulate_regions.py`: Creates `region_name` column in tabulation

### antiSMASH Link Paths

Links from `bgc_report.html` to antiSMASH results use paths relative to `main_analysis_results/{taxon}/`:
```
../../antismash_results/{taxon}/{genome}/index.html#r{record_index}c{region_number}
```
This is set in `scripts/clustering/extract_gcf_representatives.py` (`antismash_link`). If the report location changes, this depth must be updated accordingly.

### Known Issues

- **Duplicate gene names**: Some NCBI genomes have duplicate CDS feature names (e.g., `sapC`), causing antiSMASH to fail with "multiple CDS features have the same name"
- **DIAMOND memory errors**: `malloc(): corrupted top size` errors during ClusterBlast indicate memory issues; try increasing memory allocation or reducing concurrent jobs
- **NCBI dehydrated download corruption**: The `--dehydrated` download mode can produce null-filled files due to network timeouts. The module includes validation with `sync` + retry logic, but if corruption persists, delete the cached work directory and re-run
- **pyhmmer/BiG-SCAPE compatibility**: pyhmmer 0.12+ changed `profile.accession` from bytes to str, breaking BiG-SCAPE. The module pins `pyhmmer<0.11`
- **GTDB-Tk duplicate taxon labels**: GTDB-Tk normalizes genome names case-insensitively. If two genomes have names differing only in case (e.g., `MDCuke` vs `MDcuke`), GTDB-Tk will fail with `NewickReaderDuplicateTaxonError`. The `create_name_map.py` script now handles this by tracking names case-insensitively and adding numeric suffixes to duplicates
- **Run reports are scoped per taxon (fixed 2026-08-25)**: `trace`, `report` and
  `timeline` all set `overwrite = true`, so with one shared path any later run — a
  `-preview` included, which writes them while executing nothing — destroyed the previous
  run's benchmark data. They now live under `params.pipeline_info_dir`
  (`results/pipeline_info/<taxon>/`), defined once in `nextflow.config` and referenced by
  both the report scopes and `subworkflows/bgc_analysis.nf` so the two cannot drift.
  Analysing a new taxon no longer clobbers the last one; a repeat run of the *same* taxon
  still overwrites, which is intended. A clobbered trace can still be rebuilt from
  Nextflow's cache: `nextflow log <run_name> -f task_id,hash,name,status,exit,submit,start,complete,realtime,cpus,memory,peak_rss,peak_vmem,pcpu,pmem,rchar,wchar`
  (note `pcpu`/`pmem` — the template parser rejects `%cpu`/`%mem`). Note
  `Utils.sanitizeTaxon` is not callable from config scope, so the taxon transform is
  inlined there; it mirrors `lib/Utils.groovy` exactly and must stay in step.
- **Report output is reproducible as of 2026-08-25** (was not before): three sources of run-to-run drift were pinned. `generate_rarefaction_curve` now resamples from a local `random.Random(seed)` (`seed=0` default, `--seed` on `visualize_results.py`) instead of the global RNG; taxonomy `node_<id>` values use an md5 digest rather than the `PYTHONHASHSEED`-salted builtin `hash()`; and every SVG writer imports `utils/plotting.py`, which pins `svg.hashsalt`, plus passes `metadata=SVG_METADATA` to drop matplotlib's `<dc:date>` stamp. Two runs over the same data now produce byte-identical `rarefaction_curve.svg` and `bgc_report.html`. If you add a figure, import `utils.plotting` and pass `metadata=SVG_METADATA` to any SVG `savefig` or you reintroduce the drift
- **matplotlib clip-path ids are not covered by `svg.hashsalt`**: the salt pins marker and
  hatch ids, but clip paths are keyed on `(id(clippath), str(clippath_trans))` in
  `backend_svg.py` — `id()` being a CPython memory address. Addresses often repeat when the
  allocation sequence is identical, so simple figures compare equal *by luck*; anything that
  perturbs allocation order (a bigger distance matrix, different dict/set traffic) shifts
  them. `all_bgcs_biosynthetic_tree_circular.svg` diverged this way while
  `gcf_biosynthetic_tree.svg` did not. `utils/plotting.canonicalise_svg(path)` rewrites the
  ids in order of first appearance and is called after every SVG `savefig`; rendering is
  unaffected because definitions and references are rewritten together. **Testing note:**
  two calls in one Python process can reuse the same address and falsely pass — always
  compare across *separate* processes
- **Rarefaction query filters (added 2026-08-25)**: `viz/rarefaction.py` now joins
  `family` and filters `record_type = 'region' AND f.cutoff = ?`, matching the six other
  query sites. On real data the filters are a **no-op** — BiG-SCAPE assigns families only
  to region records (of 320 each of region/protocluster/proto_core/cand_cluster, only the
  regions appear in `bgc_record_family`), and both runs used a single cutoff. They matter
  only for a comma-separated `bigscape_cutoffs`, which would otherwise merge family ids
  across cutoffs into one set (1.4x inflation in a two-cutoff test).
- **Region membership is segment-based (fixed 2026-08-25)**: an origin-spanning region
  is a genuine join of two disjoint intervals — `join{[114343:121989](+), [0:32816](+)}`.
  Collapsing it to one `(start, end)` is lossy either way: the first pair drops the second
  segment, and `(min_start, max_end)` invents a span covering the gap, which on a wrapped
  region is most of the replicon. Both coupling scripts now use
  `parse_location_segments` + `cds_in_segments` and test overlap against every segment.
  **Measured on Pantoea:** exactly 5 of 320 regions have compound locations, and those
  were precisely the 5 classified `Unknown` — the coupling enzyme (SMCOG1271, Synthase)
  sits in the discarded second segment. The fix reclassifies all 5 to Synthase and changes
  nothing else, leaving zero Unknowns. `parse_location_bounds` is kept for callers that
  genuinely want a single extent, but its docstring now warns against using it for
  membership.

## Troubleshooting

```bash
nextflow clean -f -k              # Clear cache if -resume fails
rm -rf work/                      # Remove intermediate files (keeps databases)
du -sh work/                      # Check work dir size
```

**Preview a run without destroying `pipeline_info/`** — put this in a scratch config and
pass it with `-c`:

```groovy
trace    { enabled = false }
report   { enabled = false }
timeline { enabled = false }
```

```bash
nextflow run main.nf -preview -c no_reports.config --taxon "Pantoea"
```

- **GTDB-Tk OOM**: Requires 56-64GB RAM; keep `--pplacer_cpus 1`
- **Conda env conflicts**: Each tool has its own environment; don't mix in COLLECT_VERSIONS
- **NCBI download corruption**: If GENBANK_TO_FASTA fails with "No sequences found", check for null-filled files:
  ```bash
  # Find corrupted files (first bytes are null)
  for f in work/*/ncbi_dataset/data/*/genomic.gbff; do
    [ -z "$(head -c 10 "$f" | tr -d '\0')" ] && echo "Corrupted: $f"
  done
  # Fix: delete cached NCBI download and re-run
  grep "NCBI_DATASETS_DOWNLOAD" .nextflow.log | grep "workDir" | tail -1  # Find work dir
  rm -rf work/XX/XXXXXX  # Delete the cached directory
  nextflow run main.nf -resume
  ```
