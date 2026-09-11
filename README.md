# ClusterQuest

**Large-scale Organization of Secondary Metabolites**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

ClusterQuest is a Nextflow DSL2 pipeline for discovering phosphonate biosynthetic gene clusters (BGCs) across bacterial genomes. It retrieves genomes from NCBI, screens them, detects BGCs, clusters them into gene cluster families (GCFs), places them phylogenetically, and ranks the families most worth investigating in the lab — in one resumable workflow with an interactive HTML report.

![Pipeline Overview](docs/pipeline_overview.png)

## Features

- **Automated genome retrieval** from NCBI by taxon name (species, genus, family, order, or higher)
- **pepM pre-screen** that skips antiSMASH on genomes with no PEP mutase — 11% of genomes retained on Erwiniaceae at 298/298 sensitivity
- **Phosphonate BGC detection** with antiSMASH, restricted to the phosphonate rule, with KnownClusterBlast against MIBiG
- **GCF clustering** via BiG-SCAPE, with optional pepM-identity partitioning for inputs too large to cluster in one job
- **Coupling enzyme classification** of each BGC, scored against characterised references
- **Novelty scoring** that ranks GCFs for lab follow-up and publishes its components, not just a composite
- **Phylogenetic placement** with GTDB-Tk, sharded for large genome sets
- **Cross-taxon result reuse** for antiSMASH and GTDB-Tk
- **Interactive HTML report** with sidebar navigation, searchable genomes, and linked BGC regions
- **Pinned reference databases** with provenance recorded in every report
- **HPC support** via a SLURM profile

## Installation

```bash
git clone https://github.com/alpole23/ClusterQuest.git
cd ClusterQuest

conda create -n nextflow -c conda-forge nextflow
conda activate nextflow
nextflow -version
```

Or download a ZIP from the green **Code** button on GitHub and extract it.

Per-process conda environments are created automatically on first run.

## Quick Start

```bash
# Analyse a taxon end to end
nextflow run main.nf --taxon "Pantoea"

# Resume after an interruption or a code change
nextflow run main.nf --taxon "Pantoea" -resume

# Send every genome to antiSMASH instead of pre-screening
nextflow run main.nf --taxon "Pantoea" --pepm_prescreen false
```

The report lands at `results/main_analysis_results/{taxon}/bgc_report.html`.

> **Open it from disk, not from a copy.** The report reads its sibling `genomes/`
> directory and links into `results/antismash_results/`. Moving the HTML file on its
> own breaks genome search and every region link.

### SLURM

```bash
# Direct
nextflow run main.nf -profile slurm --taxon "Streptomyces"

# Or via the submission script
sbatch submit_slurm.sh "Streptomyces"

# Follow the most recent job
tail -f $(ls -t clusterquest_*.out | head -1)
```

Before your first SLURM run, edit `submit_slurm.sh` for your cluster's partition
(`#SBATCH -p`), notification email, and module loads.

Allocations under `-profile slurm`:

| Process | CPUs | Memory | Time |
|---------|-----:|-------:|-----:|
| antiSMASH (batched) | 2 | 6 GB ×attempt | 8 h |
| BiG-SCAPE | — | 64 GB | 8 h |
| GTDB-Tk | 8 | 128 GB (highmem queue) | 24 h |
| Visualization | — | 16 GB | 4 h |

antiSMASH is allocated 2 CPUs deliberately: measured `%cpu` p99 was 145%, and the
old 4-CPU default throttled how many genomes could run at once.

## Requirements

- Nextflow ≥ 23.04
- Conda or Mamba
- 16+ GB RAM (48–128 GB for GTDB-Tk)
- ~200 GB disk for reference databases, plus ~16 MB per genome analysed

`outdir` and `workDir` should be on the same filesystem so published output can be
hard-linked (see `--publish_mode`).

## Scale

Three separate limits govern how large a taxon can be analysed, and each has its own
control.

**antiSMASH is the compute cost.** The pepM pre-screen runs between genome download
and antiSMASH, using DIAMOND to find PEP mutase — `blastp` against annotated proteins,
`blastx` against DNA for genomes below `--pepm_prescreen_min_density` proteins/Mb.
Genomes with no hit cannot carry a phosphonate BGC and are skipped. On Erwiniaceae this
took 2,771 genomes down to 306 while recovering all 333 BGCs. Every verdict is written
to `prescreen_results/{taxon}/`, so what was skipped stays auditable.

**BiG-SCAPE memory is the hard limit.** Peak memory is quadratic above ~4,000 BGCs,
measured as `GB = 1.14 + 1.29e-7·n²` — about 1.9 TB at the 121,000 BGCs a million
genomes would produce. Setting `--bigscape_partition true` splits the input by pepM
identity first, which rebuilds the same GCF network (ARI 1.0000 against unpartitioned
runs) with the largest job near 84 GB. It is **off by default**: validated at 185–518
BGCs, not yet at the scale that needs it, and never run on SLURM. Below
`--bigscape_partition_threshold` BGCs it is also a net loss, because every partition
re-pays BiG-SCAPE's fixed Pfam-load cost.

**Storage is mostly duplication.** `--publish_mode link` (the default) hard-links
published output so `results/` and `work/` share inodes; published files survive
`nextflow clean -f`. With raw downloads no longer published, peak storage is ~16 MB per
genome rather than ~52 MB. Set `copy` if `outdir` and `workDir` are on different
filesystems, or if anything edits published files in place.

Short per-genome steps are batched (`--task_batch_size`, `--antismash_batch_size`) since
they run in 0.2–1.5 s and were dominated by scheduler overhead; GTDB-Tk is sharded at
`--gtdbtk_shard_size` genomes.

Projected to 1M genomes at UIUC internal rates: antiSMASH ~25,000 CPU-h, BiG-SCAPE
~1,090 CPU-h, GTDB-Tk 563 CPU-h, ~16 TB storage. Compute is not the constraint at
roughly $1,320 — BiG-SCAPE's memory requirement is, and it forces partitioning.

See [`docs/benchmark_erwiniaceae.html`](docs/benchmark_erwiniaceae.html) for the full
cost and timing analysis.

## Parameters

### Global

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--taxon` | `Erwiniaceae` | NCBI taxon name |
| `--workflow` | `full` | `download`, `bgc_analysis`, or `full` |
| `--outdir` | `results` | Output directory |
| `--publish_mode` | `link` | `link` (hard-link) or `copy` |
| `--assembly_level` | — | Restrict download, e.g. `complete` |

### pepM pre-screen

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pepm_prescreen` | `true` | Skip antiSMASH on genomes with no PEP mutase |
| `--pepm_prescreen_bitscore` | `100` | DIAMOND cutoff; true positives score 154–552 |
| `--pepm_prescreen_min_density` | `500` | Proteins/Mb below which a genome is screened by DNA |
| `--pepm_prescreen_batch_size` | `200` | Genomes per screening task |

### antiSMASH

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--antismash_minimal` | `false` | Skip domain analysis for faster runs |
| `--antismash_batch_size` | `50` | Genomes per antiSMASH task |
| `--antismash_summary_gbk` | `false` | Write the whole-genome GenBank (~11 MB each) |
| `--reuse_antismash_from` | — | Reuse results from a previously analysed taxon |

> Detection is hardcoded to phosphonate BGCs. KnownClusterBlast, `--clusterhmmer` and
> `--tigrfam` are always on; `--no-zip-output` suppresses the redundant per-genome zip.

### Clustering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--clustering` | `bigscape` | `none` or `bigscape` |
| `--bigscape_cutoffs` | `0.30` | GCF distance threshold |
| `--bigscape_mibig_version` | — | Include MIBiG, e.g. `3.1` |
| `--bigscape_partition` | `false` | Partition by pepM identity before clustering |
| `--bigscape_partition_identity` | `0.60` | Partition threshold; 0.60–0.80 is the verified window |
| `--bigscape_partition_threshold` | `10000` | Total BGCs below which partitioning is skipped |
| `--bigscape_partition_max_size` | `0` | Largest partition; 0 derives it from the memory allocation |

### Phylogeny

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_gtdbtk` | `true` | Enable GTDB-Tk |
| `--gtdbtk_bgc_genomes_only` | `true` | Only place genomes that have BGCs |
| `--gtdbtk_shard_size` | `5000` | Genomes per shard |
| `--gtdbtk_outgroup` | `g__Escherichia` | Outgroup pattern for rooting |
| `--reuse_gtdbtk_from` | — | Reuse results from a previously analysed taxon |

### Reference databases

Versions are pinned rather than tracking "latest". Every database uses `storeDir`, which
skips the download whenever the output path exists — so an unpinned URL doesn't track
anything, it freezes whatever was current on the day of the first run with no record of
what that was. Changing a pin changes the output path, which is what triggers a fresh
download. **Two runs on different database versions are not directly comparable.**

| Parameter | Default |
|-----------|---------|
| `--gtdb_release` | `226` |
| `--pfam_release` | `38.2` |
| `--taxdump_date` | `2026-08-01` |
| `--check_db_updates` | `true` — warn (never fail) when a pin is behind upstream |

Full list: [Parameters Reference](docs/parameters.md) or [`nextflow.config`](nextflow.config).

## Report

`bgc_report.html` has a left sidebar with eight sections:

| Section | Contents |
|---------|----------|
| **Overview Stats** | Genome and BGC counts, KnownClusterBlast mapping, rarefaction/Chao2 coverage |
| **BGC Novelty** | Novelty-ranked families with their component scores, linked to representatives |
| **Gene Cluster Families → Analysis** | GCF × genus heatmap ordered by GTDB phylogeny, coupling enzyme support, pepM all-by-all |
| **Gene Cluster Families → Trees** | GCF biosynthetic phylogeny, family-centre tree |
| **Context → Phylogeny** | GTDB-Tk placement, taxonomy tree |
| **Context → Genomes Search** | Searchable genome table with per-genome pages |
| **Reference → KnownClusterBlast Hits** | MIBiG matches per region |
| **Reference → Pipeline Info** | Versions, database pins, parameters, partitioning summary |

Reports are deterministic: the same inputs and database pins produce the same HTML.

### Novelty score

`NOVELTY_SCORE` ranks GCFs as `priority = distance × evidence`, where distance is
0.7 × coupling-enzyme divergence + 0.3 × pepM divergence, and evidence weights genome
count, genus spread, region intactness, and class. The components are published
alongside the composite in `novelty_ranking.tsv` and shown in the report, so a ranking
can be argued with rather than taken on faith. Families with no coupling data go in a
separate unranked bucket instead of defaulting to zero distance.

> **GCF numbers are per-run.** `family.id` is a BiG-SCAPE `AUTOINCREMENT` column — it
> records write order, not a stable biological identity, and it shifts whenever the
> genome set changes. Identify a family by its member genomes, not its number.

## Modules

| Subworkflow | Description |
|-------------|-------------|
| `DOWNLOAD_GENOMES` | Retrieve and standardise genomes from NCBI |
| `ANTISMASH_ANALYSIS` | pepM pre-screen, batched BGC detection, cross-taxon reuse |
| `CLUSTERING` | BiG-SCAPE GCF clustering, with optional partitioning and merge |
| `PHYLOGENY` | Sharded GTDB-Tk placement, with cross-taxon reuse |
| `BGC_ANALYSIS` | Region counting, coupling classification, novelty scoring, visualization |

See the [Module Reference](docs/MODULES.md) for inputs, outputs, and configuration.

## Output Structure

```
results/
├── databases/                          # Cached, version-pinned reference databases
├── ncbi_genomes/{taxon}/
│   ├── renamed_genomes/                # Standardised genome files
│   └── name_map.json                   # Assembly ID → genome name
├── prescreen_results/{taxon}/          # pepM verdict per genome (auditable)
├── antismash_results/{taxon}/          # Per-genome antiSMASH output
├── bigscape_results/{taxon}/
│   ├── {taxon}.db                      # SQLite: BGCs, families, distances
│   ├── gcf_representatives.json        # GCF data with gene diagrams
│   └── bigscape_statistics.json
├── gtdbtk_results/{taxon}/
├── pipeline_info/{taxon}/              # trace, timeline, execution report
└── main_analysis_results/{taxon}/
    ├── bgc_report.html                 # ← the report
    ├── genomes/                        # Per-genome pages (required by the report)
    ├── novelty_ranking.tsv             # GCF priority + components
    ├── region_counts.tsv               # BGC counts per genome
    ├── region_tabulation.tsv           # Per-region detail
    ├── rarefaction_curve.png
    ├── pruned_phylo_tree.nwk
    ├── taxonomy_tree.json
    ├── pepm_all_by_all/                # pepM identity vs neighbourhood similarity
    └── gcf_heatmap/
        ├── gcf_species_heatmap.{png,svg}
        ├── gcf_biosynthetic_tree.{png,svg}
        ├── phosphonate_itol_coupling.txt
        └── phosphonate_coupling_support.tsv
```

## Cross-Taxon Result Reuse

```bash
# First: a broad taxon
nextflow run main.nf --taxon "Erwiniaceae"

# Later: a subset, reusing what was already computed
nextflow run main.nf --taxon "Pantoea" \
    --reuse_antismash_from "Erwiniaceae" \
    --reuse_gtdbtk_from "Erwiniaceae"
```

**antiSMASH reuse is per-genome and works in both directions.** Each genome is checked
individually; those found in the source are copied, the rest run fresh. Partial reuse
therefore works — reusing "Pantoea" results while running "Erwiniaceae" skips antiSMASH
for the Pantoea genomes and runs only the other genera.

**GTDB-Tk reuse requires superset → subset.** Tree placement needs all genomes analysed
together, so reuse applies only when every current genome exists in the source results.
Otherwise the pipeline falls back to a full GTDB-Tk run.

Reuse compounds with the pre-screen: a *P. ananatis* run downloaded 343 genomes, screened
192 through, and ran antiSMASH zero times, recovering 225/225 BGCs.

> **antiSMASH output is not byte-reproducible.** It stamps `Run date` into every region
> GenBank, so file hashes differ between otherwise identical runs. Compare BiG-SCAPE's
> `nt_seq`, never `gbk.hash`.

## Tests

```bash
bash tests/run_tests.sh
```

Twelve checks: Nextflow-level tests for `Utils` helpers and the batched per-genome
processes (including a corrupt genome and reuse-copy fidelity), plus Python-level
compile, undefined-call, script-dependency, and report-JS-handler checks.

`check_script_deps.py` is worth knowing about. Modules invoke their Python as
`python ${projectDir}/scripts/foo.py` — an interpolated path, not a declared input — so
Nextflow's task hash never sees the script, and `-resume` would happily reuse output from
code that has since changed. Every process therefore declares
`Utils.scriptsHash(projectDir, [...])` as a **`val` input**; Nextflow hashes the
unevaluated script source plus input values, so a `# scripts-version:` comment inside the
script block invalidates nothing. This test verifies each declared list still covers what
the process imports, and rejects the comment spelling outright so it cannot come back.

## Citation

- **antiSMASH**: Blin et al. (2023) *Nucleic Acids Research*
- **BiG-SCAPE**: Navarro-Muñoz et al. (2020) *Nature Chemical Biology*
- **GTDB-Tk**: Chaumeil et al. (2022) *Bioinformatics*
- **DIAMOND**: Buchfink et al. (2021) *Nature Methods*

## Acknowledgments

Developed with assistance from [Claude](https://claude.ai), using [Claude Code](https://docs.anthropic.com/en/docs/claude-code).

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.
