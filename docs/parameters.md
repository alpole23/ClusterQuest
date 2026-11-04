# ClusterQuest parameters

Every parameter, with the measurement or reasoning behind its default.

`nextflow.config` carries a one-line comment per parameter, because that line
is what `scripts/configure_run.py` shows as the prompt and what `--help` prints.
The notes that used to sit above each value live here instead: they were being
truncated mid-word at 240 characters in the prompt, which helped nobody.

**Do not edit a default here.** The value shown is read from `nextflow.config`,
and `tests/check_parameter_docs.py` fails if the two disagree or if a parameter
is missing a section below.

Set a parameter in a params file and pass `-params-file`, never as `--flag true`:
on the CLI a boolean becomes the **string** `"true"`, and `"false"` is truthy in
Groovy. `configure_run.py` writes real JSON booleans for this reason.

## Contents

- [Global Pipeline Options](#global-pipeline-options) — 3 parameters
- [Reference database versions](#reference-database-versions) — 4 parameters
- [Task batching](#task-batching) — 2 parameters
- [pepM pre-screen](#pepm-pre-screen) — 4 parameters
- [Download Genome Options](#download-genome-options) — 7 parameters
- [antiSMASH Analysis Options](#antismash-analysis-options) — 9 parameters
- [Region Analysis Options](#region-analysis-options) — 4 parameters
- [BiG-SCAPE Clustering Options](#big-scape-clustering-options) — 9 parameters
- [BiG-SCAPE partitioning](#big-scape-partitioning) — 5 parameters
- [Missing-gene recovery](#missing-gene-recovery) — 1 parameter
- [Within-GCF annotation transfer](#within-gcf-annotation-transfer) — 3 parameters
- [GTDB-Tk Phylogeny Options](#gtdb-tk-phylogeny-options) — 9 parameters

## Global Pipeline Options

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`workflow`](#workflow) | `"full"` | no |
| [`outdir`](#outdir) | `"results"` | yes |
| [`publish_mode`](#publish_mode) | `'link'` | `--advanced` |

### `workflow`

Pipeline mode: "download", "bgc_analysis", or "full"

**Default** `"full"`  ·  **Type** string

### `outdir`

Output directory for all results

**Default** `"results"`  ·  **Type** string  ·  **Modes** bgc_analysis, download, full

### `publish_mode`

How output reaches outdir: 'link' hard-links (no extra disk), 'copy' for a separate filesystem.

**Default** `'link'`  ·  **Type** string  ·  **Modes** bgc_analysis, download, full

How publishDir moves task output into outdir. 'link' makes hard links, so a
published file and its work/ counterpart share one inode: results/ costs no
extra disk, and the published copy survives `nextflow clean -f` (verified —
the inode's link count simply drops to 1). Measured on Erwiniaceae, copying
cost 26.25 MB/genome of pure duplication, ~25 TB at a million genomes.
Requires outdir and workDir on the same filesystem; set 'copy' if they are
not, or if anything edits published files in place, since a write through a
hard link also rewrites the cached task output and corrupts -resume.

## Reference database versions

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`gtdb_release`](#gtdb_release) | `'226'` | no |
| [`pfam_release`](#pfam_release) | `'38.2'` | no |
| [`taxdump_date`](#taxdump_date) | `'2026-08-01'` | no |
| [`check_db_updates`](#check_db_updates) | `true` | `--advanced` |

### `gtdb_release`

GTDB-Tk refuses reference data older than its own MIN_REF_DATA_VERSION.

**Default** `'226'`  ·  **Type** string

Pinned deliberately rather than tracking "latest". Every database uses
storeDir, which skips the download whenever its output path already
exists — so an unpinned URL does not track anything: it freezes whatever
was current on the day of the very first run and never re-checks, with no
record of what that was. This pipeline ran on a January taxdump for seven
months that way.

Pinning makes the version a decision you can cite in a paper, and changing
a pin changes the storeDir output path, which is what triggers a fresh
download. Do not "upgrade" these casually: two runs on different database
versions are not directly comparable.

GTDB-Tk refuses data older than its own MIN_REF_DATA_VERSION; 2.6.1 wants r226.
Bump the gtdbtk conda pin before bumping this.

### `pfam_release`

current as of 2026-08

**Default** `'38.2'`  ·  **Type** string

### `taxdump_date`

NCBI monthly archive; the live file changes daily

**Default** `'2026-08-01'`  ·  **Type** string

### `check_db_updates`

Warn when a pinned database is behind upstream. Set false for offline runs.

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, download, full

Warn (never fail) when a pinned database is behind what upstream offers.
One HEAD request per database at startup; set false for offline runs.

## Task batching

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`task_batch_size`](#task_batch_size) | `100` | `--advanced` |
| [`antismash_batch_size`](#antismash_batch_size) | `50` | `--advanced` |

### `task_batch_size`

Genomes per task for the short per-genome steps, to cut scheduler overhead.

**Default** `100`  ·  **Type** integer  ·  **Modes** bgc_analysis, download, full

Genomes per task for the short per-genome steps (rename, GenBank→FASTA,
antiSMASH result copy). These take 0.2–1.5 s each, so one task per genome is
dominated by scheduler overhead; batching cuts job count by this factor.

### `antismash_batch_size`

see antismashBatchSize() in subworkflows/helpers.nf

**Default** `50`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

## pepM pre-screen

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`pepm_prescreen`](#pepm_prescreen) | `true` | yes |
| [`pepm_prescreen_bitscore`](#pepm_prescreen_bitscore) | `100` | `--advanced` |
| [`pepm_prescreen_min_density`](#pepm_prescreen_min_density) | `500` | `--advanced` |
| [`pepm_prescreen_batch_size`](#pepm_prescreen_batch_size) | `200` | `--advanced` |

### `pepm_prescreen`

Skip antiSMASH on genomes with no PEP mutase, which cannot carry a phosphonate BGC.

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, download, full

pepM pre-screen: skip antiSMASH on genomes with no PEP mutase, which cannot
carry a phosphonate BGC. On by default: validated at 298/298 sensitivity over
2,771 genomes, and it is what makes order-scale runs tractable. Every verdict
is published to prescreen_results/<taxon>/, so what was skipped is auditable.
Set false to send every genome to antiSMASH. See CLAUDE.md.

### `pepm_prescreen_bitscore`

298/298 sensitivity; TPs score 154-552, background <= 51

**Default** `100`  ·  **Type** integer  ·  **Modes** bgc_analysis, download, full

### `pepm_prescreen_min_density`

proteins/Mb below which a genome is screened by DNA

**Default** `500`  ·  **Type** integer  ·  **Modes** bgc_analysis, download, full

### `pepm_prescreen_batch_size`

see pepmBatchSize() in subworkflows/helpers.nf

**Default** `200`  ·  **Type** integer  ·  **Modes** bgc_analysis, download, full

## Download Genome Options

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`taxon`](#taxon) | `"Erwiniaceae"` | yes |
| [`assembly_level`](#assembly_level) | `null` | `--advanced` |
| [`assembly_source`](#assembly_source) | `'GenBank'` | yes |
| [`assembly_reference_only`](#assembly_reference_only) | `false` | yes |
| [`metadata_shard_threshold`](#metadata_shard_threshold) | `100000` | `--advanced` |
| [`download_batch_size`](#download_batch_size) | `25` | `--advanced` |
| [`download_max_forks`](#download_max_forks) | `4` | `--advanced` |

### `taxon`

NCBI taxon (species, genus, family, order, etc.)

**Default** `"Erwiniaceae"`  ·  **Type** string  ·  **Modes** bgc_analysis, download, full

### `assembly_level`

Restrict the NCBI download by assembly level, e.g. 'complete'. Null takes everything.

**Default** `null`  ·  **Type** ['string', 'null']  ·  **Modes** download, full

Restrict the NCBI download by assembly level, e.g. 'complete' or
'chromosome,complete'. Null takes everything. Useful when taxonomic
breadth matters more than exhaustive coverage: Streptomyces has 16,030
GenBank genomes but 1,573 complete ones spanning 1,094 organism names,
which is a better diversity sample at a tenth of the compute.

### `assembly_source`

'GenBank' (GCA_) or 'RefSeq' (GCF_), case-insensitive.

**Default** `'GenBank'`  ·  **Type** string  ·  **Modes** download, full

'GenBank' (GCA_) or 'RefSeq' (GCF_), case-insensitive. GenBank was hardcoded
until 2026-09-23 and remains the default so existing runs reproduce.

The `datasets` CLI also accepts 'all'; this pipeline REJECTS it (2026-11-02).
'all' returns both sources, which is most genomes twice -- a GCA_ and a GCF_
accession for one assembly. Measured on Brenneria: GenBank 60, RefSeq 58,
all 119, of which 59 are the same assembly under both accessions. Both twins
survive renaming, because the accession is appended when organism names
collide, so both are screened, both reach antiSMASH and both enter BiG-SCAPE
-- inflating family member counts in a way that reads as biological rather
than as double counting.

**RefSeq is the one to use at order scale.** Measured 2026-09-23:
Enterobacterales has 1,257,000 GenBank assemblies against 150,690 RefSeq
ones — 8.3x fewer. What RefSeq drops is overwhelmingly redundant clinical
isolates differing by a few SNPs, which add compute and BiG-SCAPE memory
(quadratic) without adding chemistry. That is the same argument the
component rarefaction made: new diversity adds components, depth in one
over-sequenced clade does not.

### `assembly_reference_only`

Restrict to NCBI's designated REFERENCE genome per species.

**Default** `false`  ·  **Type** boolean  ·  **Modes** download, full

Restrict to NCBI's designated REFERENCE genome per species. Enterobacterales
has 150,690 RefSeq assemblies and 640 reference ones, so this is a ~235x
sample that still spans the order -- the same argument assembly_level makes
above, taken further. Intended for pilots: measuring a clade's actual BGC
prevalence and screen retention before committing to the full set, rather
than borrowing another clade's rates.

A SELECTION filter, not a property filter, so unlike assembly_source it is
applied only when resolving the taxon to accessions. FETCH_RENAME_SCREEN
fetches those accessions by name and must not re-filter them.

REQUIRES assembly_source = 'RefSeq'. NCBI designates reference status on GCF_
accessions only, so `--reference --assembly-source GenBank` matches nothing and
the run proceeds with ZERO genomes rather than failing. Measured 2026-10-17:

  Enterobacterales --reference --assembly-source GenBank ->      0
  Enterobacterales --reference --assembly-source RefSeq  ->    643
  Bacteria         --reference                           -> 22,871

Since assembly_source defaults to GenBank, turning this on and changing nothing
else is the silent-empty-run case. Set both.

### `metadata_shard_threshold`

Above this many records, the metadata fetch is split by release date. 0 disables.

**Default** `100000`  ·  **Type** integer  ·  **Modes** download, full

Above this many records, NCBI_FETCH_METADATA splits its fetch by release
date. NCBI cannot hold one stream open long enough for a whole domain:
3,350,527 records failed four times at 0.29-1.31 GB of ~1.72 GB, while
149,539 fetch cleanly in ~4.5 min. 0 disables sharding entirely.
Below the threshold the fetch is one call, exactly as before -- an extra
`datasets` invocation costs ~0.9 s fixed, so sharding a genus would be
~4x slower for no benefit.

### `download_batch_size`

Accessions fetched per download task. Larger archives fail more often at NCBI.

**Default** `25`  ·  **Type** integer  ·  **Modes** download, full

Accessions fetched per FETCH_RENAME_SCREEN task. Sets transient disk
(batch x download_max_forks x ~8.7 MB, about 440 MB at the defaults) and,
less obviously, the FAILURE RATE.

NCBI returns `Internal error (invalid zip archive)` for a share of requests,
and the share grows with archive size: 20 accessions (~40 MB) failed 7.5% of
the time across 188 measured fetches, while a 51-accession request (~130 MB)
failed all six attempts of a real run and took the run down with it. This was
200 at first, which would mean ~500 MB archives — the worst end of that curve.
25 keeps archives near 60 MB. Raising it trades reliability for task count.

### `download_max_forks`

Concurrent genome-fetch tasks. Link-limited, and batches start failing above 4.

**Default** `4`  ·  **Type** integer  ·  **Modes** download, full

Concurrent fetch tasks. MEASURED on the COMBINED fetch+screen stage, which is
what this process actually does — a download-only probe overstates it, because
the pepM screen runs inline before the next batch starts.

Re-probed wired on 2026-10-08 at 2/4/8/12/16, 150 genomes per level on disjoint
cold GenBank slices:

    forks   MB/s   failed batches
      2     22.4        0
      4     23.0        0          <- best clean throughput
      8     24.1        1
     12     24.7        1
     16      6.6        5          <- collapses

Throughput is LINK-limited, not concurrency-limited: 6x the forks buys 10% of
the bandwidth, and batches start failing at 8. Read MB/s, not genomes/min —
each level drew a different slice and genome size varies 2.4-16 MB, so
genomes/min mostly reports slice composition.

The earlier note claimed 33.5 MB/s peaking at 2 on a "saturated link". The
shape of that conclusion was right and still is; the rate was not — it was
download-only and over WiFi. Wired, the same forks=2 setting sustains ~5x what
the RefSeq run actually achieved (70.4 genomes/min over 43 h).

4 -- the level with the best CLEAN throughput above, and the reconciliation
of three numbers that have all appeared here is the point:

  * "4 rather than 2" came from a probe on a throttled ~264 Mbit WiFi link
    that no longer exists. Right in shape, wrong in degree, wrong hardware.
  * 24 was set when forks were not the binding limit at all. Raising them
    4 -> 12 -> 24 moved the rate 134 -> 194 genomes/min and then stopped,
    because FETCH_RENAME_SCREEN reserved 8 GB and 56 GB / 8 GB held the
    executor at SEVEN concurrent tasks whatever this said. Two tuning rounds
    went into the wrong parameter.
  * 12 was set on 2026-11-02 from conf/labels.config's fork probe, which
    measured aggregate DOWNLOAD streams (11 -> 56 MB/s from 4 to 8) and is
    exactly the download-only probe the top of this note warns overstates
    the stage. REVERTED the same day: against the combined measurement it
    buys 7% more MB/s and costs the zero-failure property.

Prefer 0 failed batches to 7% more bandwidth. A bad batch costs three inline
retries and then a fall back to one accession at a time, so on a multi-day
run its wall-clock cost dwarfs the margin -- and at 16 forks the stage does
not degrade, it collapses to 6.6 MB/s with 5 of its batches failing.

No need to scale this to the machine. It is a property of the LINK, and the
executor caps it anyway: fetch tasks take 1 cpu and `executor.cpus` is
2 x cores, so a 4-core box admits 8 of these regardless of what is set here.

## antiSMASH Analysis Options

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`input_genomes`](#input_genomes) | `null` | yes |
| [`reuse_antismash_from`](#reuse_antismash_from) | `null` | `--advanced` |
| [`antismash_minimal`](#antismash_minimal) | `false` | yes |
| [`antismash_cb_general`](#antismash_cb_general) | `false` | `--advanced` |
| [`antismash_cc_mibig`](#antismash_cc_mibig) | `false` | `--advanced` |
| [`antismash_smcog_trees`](#antismash_smcog_trees) | `false` | `--advanced` |
| [`antismash_summary_gbk`](#antismash_summary_gbk) | `false` | `--advanced` |
| [`antismash_html`](#antismash_html) | `true` | `--advanced` |
| [`antismash_phosphonate_neighbourhood`](#antismash_phosphonate_neighbourhood) | `10` | `--advanced` |

### `input_genomes`

only for --workflow bgc_analysis; Path to pre-downloaded genomes directory

**Default** `null`  ·  **Type** ['string', 'null']  ·  **Modes** bgc_analysis

### `reuse_antismash_from`

Taxon name to reuse antiSMASH results from, e.g. "Pantoea"

**Default** `null`  ·  **Type** ['string', 'null']  ·  **Modes** bgc_analysis, full

### `antismash_minimal`

antiSMASH minimal mode (faster, skips domain analysis)

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `antismash_cb_general`

only when antismash_minimal = false; ClusterBlast: Compare vs antiSMASH DB

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `antismash_cc_mibig`

only when antismash_minimal = false; ClusterCompare: Advanced MIBiG scoring

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `antismash_smcog_trees`

only when antismash_minimal = false; Phylogenetic trees for BGC genes

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `antismash_summary_gbk`

Write {genome}.gbk, the whole annotated genome (~11 MB each).

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

Write {genome}.gbk, the whole annotated genome (~11 MB each). No pipeline step reads it; useful for eyeballing small sets,
costly at scale — 1,736 Pantoea genomes would add ~20 GB.

### `antismash_html`

Keep antiSMASH's per-genome HTML viewer. ~1.13 MB per genome; set false at scale.

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

Keep antiSMASH's per-genome HTML viewer (index.html, regions.js, js/, css/, images/). ~1.13 MB per genome, of which 708 kB is
byte-identical boilerplate repeated in every genome directory.
Set false for order-scale runs: it frees ~19 GB per 16,500
analysed genomes. The report's per-genome antiSMASH links stop
resolving, so leave it on whenever anyone will read the report.
Packaging only, so like the two flags above it stays OUT of
Utils.antismashParamsHash and does not invalidate reuse.

### `antismash_phosphonate_neighbourhood`

Flank kept around a phosphonate rule core, in kb. Null uses antiSMASH's own 5.

**Default** `10`  ·  **Type** ['integer', 'null']  ·  **Modes** bgc_analysis, full

Size of the flank antiSMASH keeps around a phosphonate rule core, in kb. Its own
rule file says 5, which is too small for this chemistry: on LMG 5342 the region
stops 4.2 kb short of the HiVir cluster, so its last four genes are missing from
all 215 members of that family. 10 covers the cluster with 0.8 kb to spare; 20
triples the region and adds ~31 kb of unrelated DNA. `--hmmdetection-strictness
relaxed` is NOT an alternative — `phosphonate-like`, the rule that carries a 20 kb
neighbourhood, never fires on this cluster.

The neighbourhood is symmetric, so this also pulls in the same amount of upstream
DNA, and every region-based metric downstream sees it. Set null for antiSMASH's
own value, which also restores the parameter hash results were produced under.
Changing it invalidates antiSMASH results by design: the regions differ.

## Region Analysis Options

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`run_analysis`](#run_analysis) | `true` | `--advanced` |
| [`count_per_contig`](#count_per_contig) | `false` | `--advanced` |
| [`split_hybrids`](#split_hybrids) | `false` | `--advanced` |
| [`skip_tree`](#skip_tree) | `false` | `--advanced` |

### `run_analysis`

Enable region analysis and visualization

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `count_per_contig`

Count per contig (true) or per genome (false)

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `split_hybrids`

Split hybrid types (T1PKS-NRPS → T1PKS + NRPS)

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `skip_tree`

Skip phylogenetic tree in HTML report visualization

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

## BiG-SCAPE Clustering Options

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`clustering`](#clustering) | `"bigscape"` | yes |
| [`bigscape_cutoffs`](#bigscape_cutoffs) | `"0.30"` | yes |
| [`bigscape_min_bgc_length`](#bigscape_min_bgc_length) | `2000` | `--advanced` |
| [`bigscape_alignment_mode`](#bigscape_alignment_mode) | `"auto"` | `--advanced` |
| [`bigscape_mibig_version`](#bigscape_mibig_version) | `""` | `--advanced` |
| [`bigscape_reference_dir`](#bigscape_reference_dir) | `"${projectDir}/assets/phosphonate_reference_bgcs"` | no |
| [`bigscape_classify`](#bigscape_classify) | `"category"` | `--advanced` |
| [`bigscape_include_singletons`](#bigscape_include_singletons) | `true` | `--advanced` |
| [`bigscape_mix`](#bigscape_mix) | `false` | `--advanced` |

### `clustering`

"none" or "bigscape"

**Default** `"bigscape"`  ·  **Type** string  ·  **Modes** bgc_analysis, full

### `bigscape_cutoffs`

GCF distance threshold(s), comma-separated

**Default** `"0.30"`  ·  **Type** string  ·  **Modes** bgc_analysis, full

### `bigscape_min_bgc_length`

Exclude region records shorter than this (bp) from clustering.

**Default** `2000`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

Exclude region records shorter than this (bp) from clustering.

2000 bp is about two bacterial genes, which cannot be a functional BGC by any
reading -- while partial clusters are worth keeping, since a phosphonate BGC
can be five genes. It is set low deliberately for that reason.

Why it exists: a contig-edge fragment carrying only the hallmark gene is within
the family cutoff of EVERY cluster that carries it, so affinity propagation
makes it a family centre and piles unrelated regions onto it. On GenBank
Enterobacterales a 1,001 bp single-CDS record holding nothing but pepM became
the centre of a 66-region "family" spanning 18 genera.

2000 removes exactly 2 of 1,476 regions there -- one of 1 CDS and one of 3 --
and that dissolves the artefact into 18 real families, largest 24. Raising it
buys nothing: 3000, 5000 and 10000 all give the same 18 families while
discarding 6, 15 and 52 regions respectively, some with 5 CDS.
See docs/comparisons/bigscape_cutoff_sweep/.

### `bigscape_alignment_mode`

"auto", "global", or "glocal"

**Default** `"auto"`  ·  **Type** string  ·  **Modes** bgc_analysis, full

### `bigscape_mibig_version`

MIBiG version (e.g. "4.0") to add to the distance matrix, or "" to exclude.

**Default** `""`  ·  **Type** string  ·  **Modes** bgc_analysis, full

MIBiG version (e.g. "4.0") or "" to exclude. Setting it puts all ~2,432 MIBiG clusters in the
distance matrix, including the 13 that carry a pepM.
Cost on a 334-BGC run: 2,766 BGCs, ~2.1 GB peak
(32 GB allocated), ~62 CPU-min. Off by default
because MIBiG entries join the same matrix and can
shift family composition.

### `bigscape_reference_dir`

Characterised phosphonate clusters MIBiG does not carry. Set "" to disable.

**Default** `"${projectDir}/assets/phosphonate_reference_bgcs"`  ·  **Type** string

Characterised phosphonate clusters MIBiG does NOT carry, measured against every
run so "distance to a known phosphonate BGC" becomes a real column. antiSMASH's
KnownClusterBlast cannot supply this: across 334 Erwiniaceae regions its best-hit
similarity never exceeded 15 (the pipeline's own floor), and exactly ONE region hit
any of MIBiG's 13 pepM-bearing clusters. See the directory's README.

BIGSCAPE_REFERENCES measures these on a COPY of the finished database rather than
adding them to the clustering itself: a reference inside the GCF cutoff joins a
family, and the downstream scripts read the database without knowing which records
are the dataset. Cost is 9.4 s on Erwiniaceae's 334 BGCs (10-13% of the clustering
step at 2,000-4,000), and the published clustering is unchanged.
Files must be antiSMASH-processed GenBank; set "" to disable.

### `bigscape_classify`

"", "category", "class", or "legacy"

**Default** `"category"`  ·  **Type** string  ·  **Modes** bgc_analysis, full

### `bigscape_include_singletons`

Include unclustered BGCs

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `bigscape_mix`

Allow mixing BGC classes in same GCF

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

## BiG-SCAPE partitioning

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`bigscape_partition`](#bigscape_partition) | `false` | yes |
| [`bigscape_partition_identity`](#bigscape_partition_identity) | `0.60` | `--advanced` |
| [`bigscape_partition_threshold`](#bigscape_partition_threshold) | `10000` | `--advanced` |
| [`bigscape_partition_max_size`](#bigscape_partition_max_size) | `0` | `--advanced` |
| [`bigscape_partition_memory_margin`](#bigscape_partition_memory_margin) | `0.85` | `--advanced` |

### `bigscape_partition`

Split BiG-SCAPE input by pepM identity first, for runs too large to cluster in one pass.

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

BiG-SCAPE compares every BGC against every other and its peak memory is
quadratic above ~4,000 BGCs (GB = 1.14 + 1.29e-7*n^2), reaching ~1.9 TB at
the 121,000 BGCs a million genomes would yield. Splitting the input by
pepM sequence identity first rebuilds the same GCF network: ARI 1.0000
against unpartitioned runs on Streptomyces, Erwiniaceae and a combined
set (185-518 BGCs), and at 1,302 BGCs 3 of 1,302 members move — the same
3 that move between two UNPARTITIONED runs of that input, which disagree
with each other slightly more than the partitioned run does. Largest job
around 84 GB. See docs/comparisons/bigscape_partition_equivalence/.

That holds only while no component is force-chunked. Chunking a real
942-BGC component into 400/400/142 cost 4,605 co-membership pairs and
invented 7 families (Jaccard 0.852) — so treat PARTITION_BGCS's
"exceeded --max_partition_size" warning as invalidating the run.

Off by default: validated to 1,302 BGCs, not at the scale that needs it,
and it has never run on SLURM. Turn it on deliberately.

### `bigscape_partition_identity`

pepM identity threshold for partitioning. 0.60-0.80 is the verified safe window.

**Default** `0.60`  ·  **Type** number  ·  **Modes** bgc_analysis, full

0.60-0.80 is the verified safe window. Every same-GCF pair measured had
pepM identity >= 0.901, so 0.60 leaves wide margin; 0.90 demonstrably
splits real families (23 pairs, ARI 0.974). Failure is silent — a split
family looks like two families — so err low.

### `bigscape_partition_threshold`

Total BGCs at which partitioning starts paying, not a per-partition size.

**Default** `10000`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

Total BGCs in the run (not base pairs, and not a per-partition size) at
which partitioning starts paying. Below it everything goes in one
partition and the run behaves as if partitioning were off.

Partitioning is a net LOSS at small scale: every partition re-pays
BiG-SCAPE's fixed Pfam-load cost. Measured on identical inputs,
185 BGCs took 30 s as one job and 112 s across 19 (3.7x slower);
518 BGCs took 93 s against 191 s (2.1x). The quadratic term only
overtakes that fixed cost around 10,000-12,000 BGCs.

### `bigscape_partition_max_size`

Largest partition, in BGCs. 0 derives it from the task's memory allocation.

**Default** `0`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

Largest partition allowed, as a BGC count. Leave at 0 to derive it from the
memory BIGSCAPE_PARTITION is allocated — the cap exists to avoid an OOM
kill, so it should follow the allocation rather than be guessed.

This is a safety valve, NOT a tuning knob. To run smaller partitions,
lower the memory allocation: that moves the cap and the request together,
which is the honest lever.

BUT the derivation understates, because of the INSTRUMENT rather than the
data. `/usr/bin/time -f %M` is ru_maxrss, which for a parent plus
waited-for children reports the max of any ONE process -- and BiG-SCAPE
runs 11. Tree-wide PSS measured at five sizes gives a stable ratio of
1.58x (1.51-1.65), so the honest curve is the published fit x 1.58:

    GB = 1.80 + 2.04e-7*n^2

    n = 5,000  ->   7 GB     32 GB  -> safe cap 11,163 BGCs
    n = 8,000  ->  15 GB     56 GB  -> safe cap 14,990
    n = 11,707 ->  30 GB     64 GB  -> safe cap 16,064
                              128 GB -> safe cap 22,912

So the derivation returns ~19,000 at 56 GB where ~15,000 is honest -- 27%
high, worth setting explicitly for a run near the cap but not dangerous.
An ~11,700-BGC GenBank run needs ~30 GB and fits this box WITHOUT
partitioning at all.

RETRACTED: an earlier version of this note claimed ~131 GB at 11,707 and a
safe cap of ~6,900, from a quadratic fitted over 150-1,300 BGCs where the
quadratic term is only 36% of the signal. Extrapolating that 93x in n^2
against a benchmark reaching 10,000 -- 1.62 GB at 333 rising to 14.21 GB,
the quadratic in full view -- was wrong. It also blamed BGC replication,
which the data does not support: at n=333 the pool is 1.0x, no replication
at all, and distinct BGCs still read 1.16x, the same as at n=600 with a
1.8x pool. That ratio is dataset variance, not replication.
See docs/comparisons/bigscape_partition_equivalence/.

### `bigscape_partition_memory_margin`

Fraction of the allocation to budget when deriving the partition cap.

**Default** `0.85`  ·  **Type** number  ·  **Modes** bgc_analysis, full

Fraction of the allocation to budget when deriving the cap. Deliberately
under 1.0 because of that extrapolation: costs ~8% of the cap, buys ~19 GB
of headroom at 128 GB, against an OOM that discards hours of clustering.

## Missing-gene recovery

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`recover_orfs`](#recover_orfs) | `true` | yes |

### `recover_orfs`

Call genes antiSMASH would miss: it runs its own finder only on records with zero CDS.

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

Genomes are often only partly annotated, and antiSMASH runs its own gene
finder ONLY on records with zero CDS features -- so a genome that shipped
with some genes keeps only those. Measured on Erwiniaceae: 327 of 333 BGC
regions contain a 300+ bp unannotated gap, and the lab-confirmed
phosphonolipid cluster in P. ananatis LMG 5342 was missing 9 of its 15
genes, including the transaminase that sets the product's headgroup.

Recovery runs after the pre-screen and feeds antiSMASH via
--genefinding-gff3, which merges extra calls into an annotated record.
prodigal finds ORFs de novo; a blastx against the other screened genomes'
proteins catches what prodigal truncates and supplies a name.

## Within-GCF annotation transfer

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`annotation_transfer`](#annotation_transfer) | `true` | `--advanced` |
| [`transfer_min_identity`](#transfer_min_identity) | `50.0` | `--advanced` |
| [`transfer_min_coverage`](#transfer_min_coverage) | `0.70` | `--advanced` |

### `annotation_transfer`

Propagate a product found on any GCF member to its orthologues, with provenance.

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

Half the CDS in a typical run carry no functional annotation, because many
NCBI assemblies are GenBank-only. GCF members are homologous and annotation
quality across them is very uneven — the typical family has a median member
at 0% and a best member near 100% — so a product found on any member is
propagated to its orthologues. Measured on Erwiniaceae: 41.5% -> 83.3%.

Every transferred call carries its provenance (source genome, agreement
count, number of independent sources), because the failure mode is error
propagation: one mis-annotated RefSeq gene becomes N of them, and the
agreement count looks reassuring because they share a single origin.

### `transfer_min_identity`

percent identity floor

**Default** `50.0`  ·  **Type** number  ·  **Modes** bgc_analysis, full

Orthology cuts. 50/0.70 is deliberately strict: a false merge propagates a
wrong name across a family, which is worse than leaving a gene unnamed.

### `transfer_min_coverage`

mutual alignment coverage, both directions

**Default** `0.70`  ·  **Type** number  ·  **Modes** bgc_analysis, full

## GTDB-Tk Phylogeny Options

| Parameter | Default | Asked by the wizard |
|---|---|---|
| [`run_gtdbtk`](#run_gtdbtk) | `true` | yes |
| [`reuse_gtdbtk_from`](#reuse_gtdbtk_from) | `null` | `--advanced` |
| [`gtdbtk_bgc_genomes_only`](#gtdbtk_bgc_genomes_only) | `true` | `--advanced` |
| [`gtdbtk_cpus`](#gtdbtk_cpus) | `8` | `--advanced` |
| [`gtdbtk_pplacer_cpus`](#gtdbtk_pplacer_cpus) | `1` | `--advanced` |
| [`gtdbtk_shard_size`](#gtdbtk_shard_size) | `5000` | `--advanced` |
| [`gtdbtk_min_perc_aa`](#gtdbtk_min_perc_aa) | `10` | `--advanced` |
| [`gtdbtk_skip_ani_screen`](#gtdbtk_skip_ani_screen) | `false` | `--advanced` |
| [`gtdbtk_outgroup`](#gtdbtk_outgroup) | `"g__Escherichia"` | no |

### `run_gtdbtk`

Enable phylogenetic analysis

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `reuse_gtdbtk_from`

Taxon name to reuse GTDB-Tk results from

**Default** `null`  ·  **Type** ['string', 'null']  ·  **Modes** bgc_analysis, full

### `gtdbtk_bgc_genomes_only`

Only analyze genomes with detected BGCs

**Default** `true`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

### `gtdbtk_cpus`

CPUs for identify/align steps

**Default** `8`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

### `gtdbtk_pplacer_cpus`

CPUs for pplacer (keep low, memory-bound)

**Default** `1`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

### `gtdbtk_shard_size`

see gtdbtkShardSize() in subworkflows/helpers.nf

**Default** `5000`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

### `gtdbtk_min_perc_aa`

Minimum % amino acids in MSA

**Default** `10`  ·  **Type** integer  ·  **Modes** bgc_analysis, full

### `gtdbtk_skip_ani_screen`

Send every genome through the full GTDB-Tk path, skipping its ANI pre-screen.

**Default** `false`  ·  **Type** boolean  ·  **Modes** bgc_analysis, full

Skip GTDB-Tk's ANI pre-screen, sending every genome through the full
markers -> align -> pplacer path. OFF, because the screen is a fast lane
rather than a substitute: a genome without a close GTDB reference still
takes the long path, so a miss costs wall-clock, not accuracy.

Measured twice, on the easy input and the hard one:

                         fell through   screen on   screen off
  2,530 RefSeq genomes      9 (0.4%)       50 min        --
    450 GenBank drafts     24 (5.3%)       44 min      143 min

GenBank misses the screen 13x more often -- 11.3% at contig level against
4.0% at complete -- and it still changes no call: all 426 ANI-decided
genomes matched their tree-placed classification on the full taxonomy
string, at every rank including species. The only output difference is the
`classification_method` column, which nothing in this pipeline reads.

Kept as a toggle rather than deleted so the old behaviour stays reachable
if a future GTDB release changes what the screen does.
See docs/comparisons/gtdbtk_ani_screen_genbank/.

### `gtdbtk_outgroup`

Outgroup pattern for tree rooting (e.g., "g__Escherichia")

**Default** `"g__Escherichia"`  ·  **Type** string
