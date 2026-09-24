# Before / after, per change

Each change to how this pipeline detects and characterises phosphonate BGCs is supposed
to be justified by data, not by recollection. The runs that produce that data live under
`results/`, which is gitignored and gets overwritten by the next run — so the small,
comparable core of each comparison is extracted here and committed.

`scripts/compare_runs.py` builds the per-BGC form of this from two result directories:

```bash
python scripts/compare_runs.py --results results \
    --before Erwiniaceae_pre_recovery --after Erwiniaceae \
    --outdir docs/comparisons/orf_recovery
```

It matches BGCs on genome plus region filename — never on file hash, since antiSMASH
stamps a run date into every region GenBank and no two runs agree byte for byte.

## Link the databases before running an arm into a new `--outdir`

```bash
mkdir -p results_new/databases
ln -s "$(readlink -f results/databases)"/* results_new/databases/
```

**Do this first, every time.** Every database process uses
`storeDir "${params.outdir}/databases"`, and `storeDir` skips a download only when
that exact path already exists. A second `--outdir` is therefore a second copy of
everything — **153 GB** (GTDB-Tk 139, antiSMASH 9.4, Pfam 4.5, TaxonKit 0.5), about
77 minutes of transfer before any analysis starts, and another chance for a download
to stall. One comparison run died exactly that way: the antiSMASH fetch hung and
Nextflow killed the run with `process hasn't exited`.

The databases are **version-pinned** (`gtdb_release`, `pfam_release`, `taxdump_date`)
and byte-identical across runs by construction, so sharing them is always correct.
Nothing about them depends on the taxon or the output directory.

This never affects a normal single-`outdir` user, which is why it stays a documented
step rather than a parameter: results are already namespaced by taxon *inside* one
outdir (`${outdir}/antismash_results/<taxon>/`, and `Utils.buildReusePath` resolves
within the same outdir), so analysing ten taxa downloads the databases once. It is
specifically the A/B pattern below — a second outdir per comparison — that pays.

A `params.database_dir` was considered and rejected: it would add API surface for a
problem only this project's own methodology creates, and the default would have to
stay `${outdir}/databases` anyway, since changing it would make every existing
install silently re-download 153 GB.

## What is here

| change | question | data | status |
|---|---|---|---|
| pepM pre-screen | what does it cost, what does it save, and does it generalise? | `pepm_prescreen/` | **controlled** (P. ananatis), **second clade** (actinomycetes), **genuinely held-out clade** (*B. fragilis*), + an older confounded pair |
| ORF recovery ("gene refactor") | what does it change about BGC gene content and classification? | `orf_recovery/` | measured, same taxon both sides |
| KCB vs BiG-SCAPE | can KnownClusterBlast measure distance to known clusters? | `kcb_vs_bigscape/` | measured |
| reference clusters | do they belong in the clustering, what does the pass cost, does `contig_edge` fix the boundary artefact? | `bigscape_references/` | measured |
| antiSMASH neighbourhood | does a wider flank capture the whole cluster, at what cost? | `antismash_neighbourhood/` | one genome for boundaries, family scale in **two clades** — which disagree |
| pepM screen discriminators | can more references or a profile HMM make the screen more specific? | `pepm_screen_discriminators/` | **no** — tested three alternatives, all worse |
| BiG-SCAPE determinism | are distances reproducible run to run? | `bigscape_determinism/` | measured |
| pantaphos family split | is the 186/29 split biosynthetic, and does antiSMASH's CUTOFF-chaining manufacture families? | `pantaphos_family_split/` | **yes** for the ATP-grasp (ablation merges 18 of 29); **no** for the chaining — negative result, nothing merged across 334 regions |

### pepM pre-screen — `pepm_prescreen/`

antiSMASH CPU-minutes over *Erwiniaceae*: **3,680 with the screen off, 569 with it on**
(6.5x), against a screen costing **91 CPU-minutes**. Whole-run total falls from 4,508 to
1,497 CPU-minutes. The screen sends 307 of 2,771 genomes to antiSMASH.

**Controlled, 2026-09-17 — `pepm_prescreen/controlled_pantoea_ananatis/`.** 344 *P. ananatis*
genomes, one code version, screen off vs on:

| | screen off | screen on |
|---|---:|---:|
| genomes reaching antiSMASH | 344 | 193 |
| ANTISMASH | 528.7 CPU-min | 311.3 |
| RECOVER_ORFS | 73.1 | 23.7 |
| PEPM_PRESCREEN | — | **10.5** |
| total | 603.0 | **346.0** |

**257 CPU-minutes saved for 10.5 spent**, and the screen loses nothing: 226 BGCs either way,
identical region set, all 226 base-for-base identical. The saving is smaller than the 6.5x
below because this taxon is phosphonate-rich — 56% of its genomes carry a pepM against 11%
of Erwiniaceae — so the screen pays in proportion to how dilute the taxon is.

It is not perfectly neutral, though: **4 of 226 BGCs each lost one recovered gene**
(3,227 → 3,223 CDS). `BUILD_PROTEIN_POOL` draws homology evidence from the screened genomes
only, 193 instead of 344, so ORF recovery is slightly weaker behind the screen. The BGC
sequences are identical; the annotation on four of them is not.

**The older pair below is confounded, and the confound is not small.** The two sides are different code: the
screen-off side is the 2026-08-29 benchmark run, the screen-on side is 2026-09-15, which
also batches antiSMASH and adds ORF recovery. antiSMASH CPU-minutes are comparable
between them; task counts are not. The per-task traces of the 2026-09-09 validation runs,
which *were* one code version with the screen off and on, are gone — their Nextflow cache
index files no longer exist, so they cannot be rebuilt. It is kept because it is the only
measurement at family scale — 2,771 genomes rather than 344 — and the controlled run above
is what the claim now rests on.

The validation matrix in `CLAUDE.md` records 333 BGCs and ARI 1.0000 with the screen on or
off; that matrix's own per-pair data was never stored.

### Second clade for the pre-screen — `pepm_prescreen/controlled_actinomycetes/`

323 actinomycete genomes across three documented phosphonate lineages. Ground truth from
an unscreened arm: 31 BGC-positive genomes, 32 regions.

**This is not a held-out clade, and an earlier version of this file said it was.** Six of
the seven curated pepM references are actinomycete; only HvrA is *Pantoea*. One reference
comes from *Glycomyces* sp. NRRL B-16210, which is **in this test set** and scores 828 on
its own sequence. Sensitivity excluding that self-match is 30/30. *Erwiniaceae*, with one
reference of seven, was the cross-clade test.

**Sensitivity 31/31 — nothing missed**, with 37 points of headroom between the cut at 100
and the weakest true positive at 137. antiSMASH fell 7.3x (944 → 129 CPU-min) for a screen
costing 12.7, because only 9.6% of this clade is BGC-positive against 56% of *P. ananatis*.

False positives are **25x more common** than on Erwiniaceae (23 of 292 negatives, 7.9%,
against 0.32%), and one true positive falls inside the negative range: 30 of 31 sit above
the top negative of 157, but the weakest sits at 137. A cut at 160 would give 30/31 with
zero false positives.

The false positives are full-length homologues at 31.9-37.7% identity — the isocitrate
lyase / PEP mutase superfamily, which actinomycetes carry in quantity. Adding four more
characterised actinomycete pepMs does **not** help: the weak true positive rises to 147,
still below 157, and one false positive is added. Neither identity nor coverage separates
the classes, so the discriminator is wrong rather than the reference set.

### Genuinely held-out clade — `pepm_prescreen/heldout_bacteroides/`

136 *Bacteroides fragilis* genomes. Phylum **Bacteroidota contributes none of the seven
pepM references**, and the clade is held out in BGC space too: all 5 characterised
reference clusters sit 0.82-0.93 from their nearest *B. fragilis* BGC, none inside the
0.30 GCF cutoff. Ground truth from an unscreened arm: **98 of 136 genomes BGC-positive,
143 regions, 19 GCFs** — 72% positive, against 11% for Erwiniaceae and 9.6% for the
actinomycete set, so this is a **sensitivity** test rather than a savings test.

**It found a real bug.** The screen missed 2 of 98. Both scored exactly **0.0** — no hit
at all — while antiSMASH called a phosphonate region in each. Both loci are annotated
`phosphoenolpyruvate mutase` *and* flagged `/pseudo`, and **NCBI withholds `/translation`
from a `/pseudo` CDS**, so the pepM protein never reached diamond. The screen was not
failing to recognise a pepM; it was never shown one. `--min_density` cannot catch this —
both genomes run 763 and 726 CDS/Mb against the 500 guard, because density is a
whole-genome proxy for a single-gene problem.

`parse_genome` now translates from the CDS's own coordinates when the translation is
absent, keeping internal stops as `X` rather than truncating.

| | before fix | after fix |
|---|---:|---:|
| sensitivity | 96/98 | **98/98** |
| true-positive bitscores | 0.0-567 | **342-567** |
| top negative | 255 | 330 |
| highest lossless cut | 0.0 (38 FPs) | **342, with zero FPs** |
| false positives | 1 of 38 | 2 of 38 |

So on this clade the classes separate **completely** — any cut in (330, 342] gives 98/98
with no false positives — where the actinomycete set had no such cut at all.
**Erwiniaceae is unchanged in every field** (299/299, 307 passed, 8 FPs of 2,485, true
positives 154-552), so the fix is strictly an improvement.

**The honest caveat is diversity, not count.** *B. fragilis* is one species and its true
positives cluster tightly (modal bitscore 514), so 98 positives is not 98 independent
tests. What the run establishes beyond doubt is the pseudogene blind spot, which both
earlier validations were blind to because neither clade contains a pseudogenised pepM.

Cost, for completeness: 277.4 -> 247.3 CPU-min, an **11%** saving for a screen costing
3.6 — against 6.5x on Erwiniaceae and 7.3x on the actinomycetes. The saving tracks how
dilute the taxon is, and at 72% positive there is almost nothing to remove.

### antiSMASH neighbourhood — `antismash_neighbourhood/`

On one genome, 5 kb truncates the HiVir cluster and 10 kb captures it (see `summary.json`).
At family scale the two clades disagree, which is the reason to have tested a second one.
On *P. ananatis* (`controlled_pantoea_ananatis/`, 344 genomes) the same 226 BGCs give 5
families either way with membership ARI **1.0000** — nothing split or merged. On the
actinomycete set (`controlled_actinomycetes/`, 323 genomes) the same 32 BGCs give 20
families either way but membership ARI **0.9524**: two *Glycomyces* clusters merge at 10 kb
and two *S. griseus* clusters split, four of 32 BGCs in all. So "the flank does not move GCF
structure" holds in Enterobacterales and **not** in actinomycetes; it is a small change, and
arguably a better one, but it is a change. What
changes is what the regions contain: median 16.7 → 29.2 kb, **3,223 → 5,253 CDS** (+63%),
223 of 226 BGCs gaining genes, for 9% more antiSMASH CPU. The flank is symmetric, so some of
those genes are upstream context the cluster does not contain.

### ORF recovery — `orf_recovery/`

Same taxon, before and after, 333 vs 334 BGCs:

- **117 BGCs gained genes; 268 genes in total** (5,186 → 5,463 CDS inside BGC regions)
- **19 families → 18**, ARI 0.9990 — the change is one family, not a reshuffle
- no product-class changed
- *P. ananatis* LMG 5342 region 2: **14 → 28 genes, family 18 → 4** — the false singleton
  the missing annotation had manufactured, dissolving into a real family
- *Winslowiella iniecta* B149 and B120: 12 → 23 genes, family 11 → 9

`per_bgc.tsv` carries every BGC's gene count and family on both sides; `per_gcf.tsv`
carries each run's novelty ranking.

## Gaps

Changes with no stored before/after data. Listed so the absence is visible rather than
assumed:

- **a comparison against alternative tools** (BiG-SLiCE, plain antiSMASH + BiG-SCAPE)
- **NOVELTY_SCORE**, and the later switch to within-run isolation
- **GCF_ANNOTATION_TRANSFER** — the current run publishes its output, but no "before"
- **Pfam gene categories and the consensus cluster**
- **the phosphonate rule set** — `bgc_rules.py --validate` exists; its output is not kept
- **headgroup prediction**

Also noted while assembling this: `region_tabulation.tsv` has **empty KCB columns for all
334 regions** in both the current and the pre-recovery run, although antiSMASH did write
`knownclusterblast/` output (248 regions carry a significant hit). The report's KCB fields
are therefore blank. Flagged here, not investigated.
