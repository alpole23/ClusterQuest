# ClusterQuest demonstration runbook

Written 2026-11-02 for a same-day demo. Two runs are set up: an **offline** one that
needs no network, and a **live** one that fetches from NCBI. Run the offline one if
you only want to show one thing.

---

## What was broken before this was set up

Two things would have stopped a demo cold, both now handled:

1. **The reference databases were gone from this box.** `results/databases/` was
   empty and `work/` was 4 KB, so every database symlink inside
   `results_entero_refseq/` and `results_entero_norecover/` was dangling. Nothing
   that calls antiSMASH or BiG-SCAPE could have run.

2. **One stale conda environment on this box cannot verify TLS.** The cached
   `bioconda::antismash` environment has a CA path baked into its OpenSSL pointing
   at `work/conda/`, where the conda cache lived before it moved to
   `~/.clusterquest/conda`. `work/` has since been cleaned, so that path resolves to
   nothing and every HTTPS download from inside that environment fails with
   `CERTIFICATE_VERIFY_FAILED`. `demo/run.sh` exports
   `SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt`, which sidesteps it.

**Neither point is a pipeline bug.** An earlier draft of this file claimed both
were, and that was wrong on both counts:

- **`hmmpress` is not missing.** It is in that environment along with 23 other
  HMMER binaries. antiSMASH locates its tools by probing `PATH`, and the failure
  came from invoking the downloader by absolute path without activating the
  environment. Under Nextflow, which activates it, the downloader runs clean.
  Nothing needs adding to `conf/conda.config`.
- **A fresh machine does not hit the TLS problem.** Verified by creating a
  throwaway environment at the current cache directory: its `openssl_cafile`
  resolves correctly and HTTPS verification succeeds. The failure needs an
  environment that was built at one prefix and is being served from another.

The permanent fix for this box, when there is no meeting pending, is to delete the
stale environment from `~/.clusterquest/conda` and let Nextflow re-solve it. **Do
not do that today** — a re-solve can land a different antiSMASH or BiG-SCAPE
version, which is not a thing to discover an hour before presenting.

---

## Before the meeting

Everything below is already done unless a step says otherwise. Verify with the
commands given; each should print what it says it will.

### 1. Confirm the databases are in place

```bash
bash /home/alexp/pipeline/demo/check.sh
```

Prints one line per database with OK or MISSING. All four must read OK:
antiSMASH, Pfam 38.2, taxonkit, and the demo genome set. If Pfam reads MISSING,
the download did not finish — see *If Pfam is missing* below.

### 2. Do one warm-up run

**Do this before the meeting, not during it.** The first run of any Nextflow
pipeline solves conda environments and stages files; the second run of the same
thing finishes in seconds because `-resume` reuses the cache.

```bash
bash /home/alexp/pipeline/demo/run.sh offline
```

**Already done for you** — it completed in 5 min 24 s on 2026-11-02 and the results
are in `results_demo/`. Run it again anyway if you want to see it happen, or to
confirm nothing has drifted; `-resume` makes the repeat replay from cache in
seconds rather than redoing the 5 minutes. Measured: the second run took **11
seconds**, with 30 of 31 processes cached (only the report regenerates).

### 3. Open the report once so you know where it is

```bash
ls /home/alexp/pipeline/results_demo/main_analysis_results/Demo/bgc_report.html
```

Open it in a browser and leave the tab up. If the live run breaks in front of
people, this is what you switch to.

---

## Running it

```bash
bash /home/alexp/pipeline/demo/run.sh offline
```

```bash
bash /home/alexp/pipeline/demo/run.sh live
```

Both write to `results_demo/`. Both pass `-resume`, so re-running after an
interruption picks up where it stopped instead of starting over.

To stop a run: `Ctrl-C` once and wait a few seconds for Nextflow to shut its
tasks down. Pressing it repeatedly leaves orphaned processes.

---

## What each run shows, and the numbers to expect

### offline — 24 genomes, no network

`demo/genomes_offline/` holds 24 genomes already on disk, chosen to span six
gene-cluster families plus four singletons from the Enterobacterales RefSeq run.
Six are *Pantoea ananatis*, which carry the pantaphos cluster — the largest
family in that run at 186 members.

What it exercises: pepM pre-screen → ORF recovery → antiSMASH → BiG-SCAPE → report.

**These are measured, from a real run on this box on 2026-11-02, not estimates:**

| Step | Result |
|---|---|
| wall time | **5 min 24 s** (cold, first run) |
| peak memory | 1.8 GB |
| Nextflow processes | 31 completed, 0 failed |
| genomes in | 24 |
| pass the pepM screen | 24 — **all of them** |
| BGC regions found | **24, all phosphonate** |
| gene-cluster families at c0.30 | **10: sizes 8, 5, 3, 2, and 6 singletons** |

Note the family count came out as 4 multi-member families plus 6 singletons, not
the 6-plus-4 those same genomes showed inside the full Enterobacterales run. That is
the caveat below doing exactly what it warns about — BiG-SCAPE families depend on
the whole input set, and 24 regions on their own cluster differently from 24 regions
among 1,302. Nothing is wrong; just quote the number above, not the RefSeq one.

**Say this out loud:** all 24 pass the screen because they were *selected* from
genomes that already passed it. This run demonstrates detection and clustering,
not the screen's selectivity. Use the live run for that.

### live — 58 *Brenneria* genomes from NCBI

| Step | Expect |
|---|---|
| RefSeq genomes for *Brenneria* | 58 |
| pass the pepM screen | **25 (43%)** |
| BGC regions found | ~25 |

The 58 and the 25 come from this pipeline's own committed screen verdicts for the
Enterobacterales RefSeq run, so the screen numbers should reproduce closely.

**Be careful with the family count.** In the RefSeq run those 25 *Brenneria*
regions fell into 5 families (16, 4, 2, 2, 1), but they were clustered alongside
1,277 other regions. BiG-SCAPE families are a property of the whole input set, so
clustering 25 regions on their own can give a different answer. Do not promise a
number here. Separately, `bigscape_min_bgc_length = 2000` was added after that
run, so region counts can differ by a region or two.

---

## What is deliberately switched off

| Setting | Why |
|---|---|
| `run_gtdbtk: false` | The GTDB-Tk reference database is ~110 GB and is not on this box. |
| `skip_tree: true` | The report's phylogenetic tree needs GTDB-Tk output. |
| `antismash_minimal: true` | Skips `--clusterhmmer` and `--tigrfam`. Roughly 3× faster per genome, and it drops a 400 MB Pfam 35.0 download. Trade-off: no Pfam domain annotation inside regions. |
| `bigscape_partition: false` | This is the default. Partitioning does nothing below 10,000 BGCs, and this demo has ~25. |
| `check_db_updates: false` | Stops the run reaching out to check database versions while an audience waits. |

If someone asks to see domain annotation, `antismash_minimal: false` turns it on,
but it needs Pfam 35.0 downloaded first and makes antiSMASH ~3× slower. Not a
during-the-meeting change.

---

## If something goes wrong

**`Cannot find nextflow` / Java errors.** `demo/run.sh` sets `JAVA_HOME` and `PATH`
itself. If you are running `nextflow` by hand instead, you need:

```bash
export JAVA_HOME="$HOME/miniconda3/envs/nextflow"; export PATH="$PATH:$HOME/miniconda3/envs/nextflow/bin"
```

**`CERTIFICATE_VERIFY_FAILED`.** A database download is being attempted from the
stale antiSMASH environment and `SSL_CERT_FILE` is not set. `demo/run.sh` sets it;
running `nextflow` by hand does not. It only ever bites on a database download, so
if the databases are cached you will never see it.

**antiSMASH tasks fail immediately.** Its database directory is incomplete. Re-run
the downloader — note that `PATH` must include the environment's own `bin`, or
antiSMASH will not find `hmmpress` and will fail at the final pre-build step:

```bash
AS=/home/alexp/.clusterquest/conda/env-b7173d34e69e6dfd50c3a7b5ba6c44f6; SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt PATH="$AS/bin:$PATH" "$AS/bin/download-antismash-databases" --database-dir /home/alexp/pipeline/results_demo/databases/antismash
```

**Careful:** this downloader fetches Pfam 35.0 **second**, before any of the real
databases, and it does so with a single connection that EBI throttles to ~20 KB/s
— about 4 hours for 279 MB. If `<db>/antismash/pfam/35.0/Pfam-A.hmm.gz` is already
present and matches its checksum, the downloader skips it and goes straight to the
fast host. That file is pre-placed, so the command above should run through. If it
ever starts crawling on `Downloading Pfam-A.hmm.gz: 0.0x%`, kill it, fetch that one
file with `aria2c -x16 -s16` into that exact path, and re-run.

**If Pfam is missing.** BiG-SCAPE cannot run without it. Fetch and compile it:

```bash
aria2c -x16 -s16 -k1M -d /home/alexp/pipeline/.demo_dbstage/pfam38.2 -o Pfam-A.hmm.gz --continue=true "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz"
```

Then run `bash demo/finish_pfam.sh` to unpack, `hmmpress` it, and move it into
`results_demo/databases/pfam38.2/`. Budget ~12 minutes for the download plus a few
for `hmmpress`. **Do not** let the pipeline's own `DOWNLOAD_PFAM` do this — it uses a
single `wget` connection, which EBI throttles to ~20 KB/s, so it would take ~6 hours.

**A run hangs with no output.** Nextflow is probably solving a conda environment.
Watch `.nextflow.log` in another terminal:

```bash
tail -f /home/alexp/pipeline/.nextflow.log
```

---

## One honest caveat about the demo as a whole

These runs show the pipeline working end to end on 24 and 58 genomes. They do not
demonstrate the scalability claims — the pre-screen's value shows up at 150,000 to
1,350,000 genomes, and BiG-SCAPE partitioning does nothing at all below 10,000
BGCs. For those, point at the committed run records under `docs/comparisons/` and
the figures in `docs/figures/`, not at a live run.
