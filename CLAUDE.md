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

**Note:** Detection is hardcoded to phosphonate rule only (`--hmmdetection-limit-to-rule-names phosphonate`). `--cb-knownclusters`, `--clusterhmmer`, and `--tigrfam` are always enabled.

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
│   └── gcf_representatives.json # GCF data with KCB hits and gene diagrams
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
            ├── Decarboxylase/   # Decarboxylase + Decarboxylase-Nucleotidyltransferase share one tree
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

### Measured: Erwiniaceae (2,548 genomes, 2026-06)

Local run, 5 h 01 m wall / 11.5 h task time over 7,951 tasks. Reference figures for
sizing resources — re-measure after any change to the antiSMASH parameter set.

| Process | Tasks | Mean | Max | %cpu (p50 / p99) | Peak RSS (p90 / p99) |
|---------|-------|------|-----|------------------|----------------------|
| `ANTISMASH` | 950 | 30 s | 9 m | 117% / 145% | 1.2 GB / 4.7 GB |
| `GTDBTK_CLASSIFY` | 1 | 1 h 31 m | — | 164% | 103 GB |
| `NCBI_DATASETS_DOWNLOAD` | 1 | 1 h 08 m | — | 13% | 24 MB |
| `BIGSCAPE` | 1 | 68 s | — | 473% | 13.8 GB |

Notes: antiSMASH never used two full cores, hence the 2-CPU override. GTDB-Tk's 103 GB
peak exceeds the 48 GB `process_high_memory` default — it only survives locally because
cgroups aren't enforced; the SLURM profile's 128 GB is the real requirement. The NCBI
download is serial and network-bound, so its runtime is unaffected by CPU allocation.

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
- **GCF Analysis**: GCF biosynthetic NJ tree (embedded as base64), dynamic coupling enzyme class table, BiG-SCAPE clustering statistics and GCF visualization
- **Novel BGCs**: BGC regions without KnownClusterBlast matches
- **KCB Hits**: Known cluster matches grouped by MIBiG entry

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
Decarboxylase-Nucleotidyltransferase 27, Decarboxylase 22, Transaminase 6,
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
| TPP+NTP | TPP_enzyme_C + NTP_transf_3 rules | → phosphonolipid (CDP-pathway) | GCF-5 | 84 |
| Ppd | SMCOG1055 (ThDP-decarboxylase) | → 2-phosphonoacetaldehyde → 2-AEP | GCF-1/8 | 72 |
| PalB | SMCOG1019 (Aminotran_1_2/PF00155) | → phosphonoalanine | GCF-7 | 20 |
| Unknown | — | — | — | 4 |

**Region boundary reading:** CDS scanning is limited to the region feature's extent, read
with `parse_location_bounds(..., span=False)` — the first coordinate pair only. This is
deliberately narrower than what `bgc_coupling_tree.py` uses (see Known Issues); widening
it reclassifies BGCs whose region feature has a compound location.

**Key insights:**
- Classification maps almost perfectly onto BiG-SCAPE GCF families — coupling enzyme type is the primary determinant of GCF membership.
- The `Fe-ADH` rule-based marker (iron-containing alcohol dehydrogenase / 2-Hacid_dh_C) is antiSMASH's marker for the phosphonopyruvate reductase (→ phosphonolactate) pathway.
- GCF-5 in that run (TPP+NTP) confirmed as **phosphonolipid BGCs**: Ppd-type ThDP enzyme + two NTP_transf_3 cytidylyltransferases + CDP-alcohol phosphatidyltransferases + Asn_synthase (CDP-phosphonate pathway). Well-annotated NCBI genomes explicitly label the ThDP enzyme as "phosphonopyruvate decarboxylase".
- AEP-pathway BGCs (GCF-1/8 in that run) use Ppd as coupling enzyme regardless of tailoring enzymes downstream.

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
- This handles divergent sequences that escape SMCOG thresholds (e.g. Ppd-CDP ThDP decarboxylases, which carry `TPP_enzyme_C` but not `SMCOG1055`)

**CLASS_MARKERS — extraction markers per class:**

| Class | Marker type | Marker | Rationale |
|-------|-------------|--------|-----------|
| Synthase (FrbC-like) | smcog | SMCOG1271 | HMGL-like phosphonomethylmalate synthase |
| Decarboxylase / Decarboxylase-Nucleotidyltransferase | domain | TPP_enzyme_C | Both classes carry this; the nucleotidyltransferase variant lacks SMCOG1055 |
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
