# Paper figures

Regenerate any of these with the script named in its section; all read committed
data under `docs/comparisons/` and `docs/benchmark_data/`, so a figure can always
be traced back to the run that produced its numbers. Output is deterministic —
`scripts/figures/figure_style.py` pins matplotlib's hash salt and drops the SVG
timestamp, so re-running gives byte-identical files.

```bash
python scripts/figures/fig1_orf_recovery.py --outdir docs/figures
python scripts/figures/fig2_prescreen.py    --outdir docs/figures
python scripts/figures/fig3_partitioning.py --outdir docs/figures
python scripts/figures/fig4_scalability.py  --outdir docs/figures
python scripts/figures/fig5_architecture.py --outdir docs/figures   # 2 figures
```

## Figure 1 — window x recovery (`fig1_window_x_recovery.py`)

*P. ananatis* LMG 5342's two phosphonate clusters, each at all four combinations
of neighbourhood (5 / 10 kb) and ORF recovery (off / on), drawn on genome
coordinates so the four captured intervals can be overlaid.

| | region 1 (pantaphos) | region 2 (phosphonolipid) |
|---|---:|---:|
| 5 kb, deposited only | 10 | 14 |
| 5 kb, + recovery | 14 | 28 |
| 10 kb, deposited only | 19 | 22 |
| **10 kb, + recovery** | **25** | **37** |

**A 2x2 rather than a sequence, because the two are not independent.** Recovered
genes can match the detection rule, which extends the rule core, which moves the
flank: region 2's 3' boundary shifts **1,566 bp from recovery alone**, at *both*
neighbourhood settings. Presenting one change as happening "before" the other
would assert an independence the data does not have.

**The two clusters in one genome need opposite fixes**, which is the argument
for making both changes. Pantaphos gains more from the window (+9) than from
recovery (+4) — 4.2 kb of the HiVir cluster lies outside the 5 kb flank
entirely. The phosphonolipid gains more from recovery (+14) than from the
window (+8). Neither change alone would have been enough for both.

Coordinates come from each run's `region_tabulation.tsv`; the region GenBanks
are individually re-based to 0 and cannot be overlaid without them.

## Figure 2 — ORF recovery across clades (`fig1_orf_recovery.py`)

antiSMASH runs gene finding only on records with **zero** CDS features, so a
GenBank deposit that annotates *some* of its genes is trusted for all of them.

**Panels A–C are three clades, and they do not agree — which is the point.**

| panel | clade | gain | what was recovered |
|---|---|---|---|
| A | *P. ananatis* LMG 5342 | 14 → 28 | the AEP transaminase, **both** CDP-alcohol phosphatidyltransferases, phosphocholine CT, an MFS transporter — the genes that say what the cluster makes |
| B | *S. griseus* | best of its regions | see the rendered panel |
| C | *B. fragilis* BFG-525 | 27 → 32 | five hypotheticals; the core cluster was already annotated |

Across all 143 *B. fragilis* regions, **82 of 99 recovered genes are
"hypothetical" and only 3 get any antiSMASH functional call**. So the claim is
not "recovery adds genes" but the sharper one: **recovery is targeted rather
than indiscriminate.** It transforms a 2012 deposit annotating 6 of 15 genes,
and is near-silent on modern complete assemblies.

**Panel D exists because panel A is an outlier.** LMG 5342 gained 14 genes and
is **rank 1 of 334**; the median gain among regions that gained anything is
**+2**, and 65% gained nothing. Showing +14 alone would invite the reader to
take it as typical. Ranks 2 and 3 are *W. iniecta* B149 and B120 at +11 — the
clusters the lab characterised were among the worst annotated in the run.

*Streptomyces* is **S. griseus, not S. hygroscopicus.** The bialaphos lineage
yields only 2 phosphonate regions across 39 genomes — which independently
reproduces its count of 2 in the actinomycete comparison — too few to choose an
example from. Finding that also surfaced a bug: `PEPM_ALL_BY_ALL` exited 1 below
three pepM sequences, killing a run after detection and clustering had finished.

Genes are coloured by antiSMASH's own `gene_kind`, not by product keyword. The
keyword version was tried first and inverted on the two lab-confirmed clusters,
because `serine hydroxymethyltransferase` matches `methyltransferase`. Region
GenBanks carry no Pfam accessions (checked: 0 of 28 CDS), so
`utils/domain_functions` has nothing to read at this level, and `unclassified`
covering 15 of 28 genes is part of the point — the CDP-alcohol
phosphatidyltransferases that decide the product are among the genes antiSMASH
does not classify.

Note `+13 recovered` against `14 → 28`: recovered genes lengthen the cluster, so
the region boundary moves and swept in one gene that was always annotated. The
figure reports both rather than letting one stand for the other.

Generating the "before" arm also surfaced a second bug: with
`recover_orfs = false` every genome in an antiSMASH batch was handed the same
`NO_RECOVERED_GFF` placeholder, and Nextflow refuses a batch whose input files
collide on name — so the ORF-recovery ablation could not be re-run at all.

## Figure 2 — pepM pre-screen

**Panel B is a negative result, deliberately.** The screen cannot *find* more
BGCs; it only decides which genomes antiSMASH sees. Its design goal is to be
detection-neutral while removing most of the compute, so the evidence that
matters is that the BGC count is **unchanged**.

Panel A is ordered by BGC prevalence rather than by clade, because the saving is
monotonic in it: 9.0× at 9.6% positive, 3.0× at 11%, 1.7× at 56%, 1.1× at 72%.
Bars are **total pipeline CPU-minutes**. Elsewhere in this repository the
*Erwiniaceae* saving is quoted as 5.1×, which is the antiSMASH-only figure
from the controlled 2026-11-01 pair on a `realtime × %cpu` basis (1,984 → 386
CPU-min). The superseded uncontrolled pair read 3,680 → 569 and 6.5×, which is a
`realtime × cpus` figure from arms that differed in code; the controlled pair on
that same basis is 5.49×. On whole-run totals it is 4.6×, and this figure uses
totals throughout so the four clades are comparable.

The *B. fragilis* bar is the one place neutrality failed: the screen returned
141 of 143 regions because a pseudogene-flagged pepM carries no `/translation`.
After the fix it returns 143. The corrected bar is **not drawn** — 143 against
141 is 1.4% of that axis and would be an invisible mark pretending to be
evidence — so it is stated in words instead.

## Figure 3 — BiG-SCAPE partitioning

The honest framing, and the one the data supports: partitioning costs
**1.7–3.7× more compute** at every scale measured, because each partition re-pays
BiG-SCAPE's fixed Pfam-load cost. What it buys is memory — **2,792 GB** in one job
at the 117,000 BGCs a million genomes yields, against **1,396 GB** for the largest
partition. A 2.0× saving.

**That was 22.4× until 2026-11-02, and the correction is worth understanding.**
The panel used `PARTITIONED_SHARE = 0.21` — the *Streptomyces* largest-partition
share, 38 of 185 BGCs — against a BGC count projected from the *Erwiniaceae* rate.
Share and population came from different taxa. Every share ever measured, all at
the 0.60 pepM identity the pipeline ships:

| run | largest partition | share |
|---|---:|---:|
| *Streptomyces* 185 BGCs | 38 | 20.5% |
| Combined (Erw + Strep) 518 | 236 | 45.6% |
| *Erwiniaceae* 5 kb, 333 | 236 | 70.9% |
| *Erwiniaceae* 10 kb, 334 | 236 | **70.7%** ← matches the projection's basis |
| Enterobacterales RefSeq, 1,302 | 942 | 72.4% |

Real taxa cluster at ~71%: one dominant pepM lineage holds most of the BGCs, and
that is what caps the saving. 20.5% was the most favourable of the five, measured
on the one test set with no dominant component.

**"1.7–3.7×" is compute, not wall time.** `bigscape_seconds_partitioned` is
accumulated in a serial loop over partitions (`total_s += secs` in
`bench_bigscape_partitioned.py`), so it is the sum of N subprocess elapsed times,
while the one-job arm is a single process's elapsed time. No partitioned *wall*
time is committed anywhere in this repository; under Nextflow the partitions run
concurrently, so the true elapsed penalty is smaller and unmeasured. Panel A said
"slower" until 2026-11-02.

**The memory figures carry a 1.58× correction, applied 2026-10-08.** The
benchmark recorded `/usr/bin/time -f %M`, which is `ru_maxrss`: for a parent plus
waited-for children it reports the max of any ONE process, and BiG-SCAPE runs
eleven. Summing PSS across the process tree is 1.58× higher, stable at
1.51–1.65 over five sizes, so both the plotted points and the fit are scaled by
it — scaling only the curve left it floating above its own data at a 63.5%
residual. See `docs/comparisons/bigscape_partition_equivalence/`.

**Panel C's equivalence is conditional**, and the caption says so. At the natural
pepM cut the family counts are identical and ARI is 1.0000; forcing a split of a
component larger than `bigscape_partition_max_size` breaks real families —
measured at 1,302 BGCs as 4,605 co-membership pairs broken and 7 families
invented. The cap is a safety valve, not a tuning knob.

Panel B draws the memory fit **solid over its measured range and dashed where it
is extrapolated**. The quadratic was fitted on 1,500–4,000 BGCs; below ~4,000
peak memory looks flat and the fit does not describe it, and 121,000 is a 12×
reach. The figure must not read as 12× more measurement than exists.

## Figure 4 — scalability on one workstation (`fig4_scalability.py`)

Why the two scalability mechanisms exist: this pipeline was developed and run on a
single 16-core / 56 GB box with no cluster, and the pepM pre-screen and BiG-SCAPE
partitioning are what make a large taxon fit on it.

**The three panels differ in evidential strength, and the figure says which is
which rather than letting a uniform style imply they are equal.**

| panel | claim | strength |
|---|---|---|
| A | pre-screen cuts antiSMASH CPU **5.5×** on *Erwiniaceae*, **443×** on Enterobacterales | *Erwiniaceae* both arms measured; the Enterobacterales without-screen bar is **hatched** — that arm was never run |
| B | partitioning saves **10 GB / 1.8×** at 10,000 BGCs | projected from the measured memory curve and the measured partition shares; **no partitioned run has ever had its memory measured** |
| C | and costs **1.7–3.7×** more BiG-SCAPE compute | measured, but cumulative compute, not elapsed time |

**Panel A's claim is not "the screen saves 5×".** The saving is `1/retention`, and
retention is a property of the taxon: *Erwiniaceae* keeps 10.76% of genomes and
saves 5.5×, Enterobacterales keeps 0.113% and saves ~443×. A 95× difference in
retention gives an 81× difference in saving, so a single number for "what the
screen saves" would be quoting *Erwiniaceae*.

**Panel A is `realtime × cpus` throughout, including *Erwiniaceae*, which is
reported on `realtime × %cpu` everywhere else.** The GenBank Enterobacterales run
has no committed trace and its per-process CPU figures are only recoverable as
`realtime × cpus` — BIGSCAPE reports 1.0 CPU-h over a 7 m 21 s realtime, which is
8.16 effective cores against its 8 allocated, a utilisation no real process
reaches. Putting the two taxa on different bases in one panel is precisely the
error this project has already published twice. On `× %cpu` the *Erwiniaceae*
ratio reads 5.14× rather than 5.49×.

**Panel B refuses to leave the measured range.** It interpolates between observed
points rather than using the quadratic fit, so every gigabyte it prints sits
between two sizes that were actually run. The shaded band spans all five measured
partition shares (20.5–72.4%) and the solid line is the 72.4% that the largest real
taxon measured — because the saving goes as `1/share²`, so the band is 5.6× at one
end and 1.8× at the other and a single multiple would be a choice, not a result.

**The honest limit, stated in the panel:** BiG-SCAPE's memory is flat below ~4,000
BGCs and quadratic above it, and every run ever partitioned is at or below 1,302
BGCs — inside the flat region. The measured memory saving from partitioning, on
every run on record, is between 0.00 and 0.43 GB. The benefit is real arithmetic
from a measured curve; it has never been demonstrated.

**Panels B and C are split rather than drawn on one pair of axes.** They are two
consequences of one intervention in different units, and a dual axis would invite a
reader to compare a gigabyte against a second.

## Figure 5 — the three mechanisms, and S1 — full architecture (`fig5_architecture.py`)

One generator, two figures, so they cannot drift:

| file | role | contents |
|---|---|---|
| `fig5_mechanisms` | main body | three panels, one per mechanism |
| `figS1_architecture` | supplemental | the stage spine, artefacts, per-stage measurements and three proportional strips, plus the same three mechanisms |

Every number is read at render time from the Enterobacterales **RefSeq** run's own
records — `summary.json` for the funnel, `pipeline_trace.tsv.gz` for every
per-stage task count, elapsed span and CPU figure, and `disk_samples.tsv` for every
peak. Nothing is typed in by hand, so the figures cannot drift from the run.

**Why RefSeq and not the 8.9× larger GenBank run.** These figures sat on GenBank
from 2026-10-28 until 2026-11-02. That run has no committed trace and no disk or
memory sampler, so it can support a funnel and a CPU share and nothing else — and
panel A's brief is real data for *every* step. Taking per-stage wall and memory
from RefSeq while taking the funnel from GenBank would put two runs in one spine,
which is the same error class as mixing two CPU bases.

**The funnel:** 150,690 genomes → 1,309 pass the pepM screen (0.87%) → 1,303
regions → 72 families, so 99.13% is discarded before anything expensive runs.

**Panel A draws three strips, because "which stage is expensive" has three
different right answers.** A single strip would have to pick one, and an earlier
version picked the one that supported its caption.

| stage | tasks | elapsed | CPU used (`× %cpu`) | CPU reserved (`× cpus`) |
|---|---:|---:|---:|---:|
| `DOWNLOAD_GENOMES` | 6,032 | 32.5 h (**75.4%**) | 42.1 (40.4%) | 252.0 (**68.2%**) |
| `ANTISMASH_ANALYSIS` | 35 | 5.5 h (12.8%) | 51.6 (**49.5%**) | 76.7 (20.8%) |
| `CLUSTERING ∥ PHYLOGENY` | 19 | 4.9 h (11.4%) | 10.2 (9.8%) | 39.5 (10.7%) |
| `BGC_ANALYSIS` + report | 12 | scattered | 0.3 (0.3%) | 1.0 (0.3%) |

Fetch-and-screen wins on elapsed time and on reserved cores but loses on CPU
actually used: it is network-bound on four allocated cores, so *reserving* them is
most of what it does. BiG-SCAPE — quadratic in both time and memory — is 0.7% of
CPU used, because the funnel hands it 1,303 regions out of 150,690 genomes.

**Peak memory and disk are from the 60 s sampler, never from Nextflow's
`peak_rss`.** That column sums RSS across the process tree: it reads 97 GB for
GTDB-Tk and 18.7 GB for a BiG-SCAPE job that `ru_maxrss` puts near 2 GB, on a
machine with 56 GB of RAM. Per stage the sampler gives 2 / 19 / 13 GB, and peak
work-directory occupancy of 21 / 34 / **150** GB — 95 GB of that last figure is
pplacer scratch inside the single GTDB-Tk task, set by GTDB's reference data rather
than by how many genomes reach it.

**The 72 families carries a caveat the figure prints.** It was measured at
`bigscape_min_bgc_length = 0`, the default until 2026-10-30. The shipped default is
2,000 bp, this run contains exactly one region below it (1,001 bp, one CDS), and it
was never re-clustered. On the GenBank run removing the equivalent record moved the
count 112 → 126, so this is not a rounding difference.

**The strips are proportional and the stage boxes are not.** A 0.3% stage cannot be
drawn as a box with a label in it, so correspondence is carried by colour — each
stage box takes a tinted band in its segment's hue. Segment labels are *measured*
for fit rather than gated on a width threshold: a fixed threshold dropped a 9.8%
segment as too narrow while letting a 10.7% one overflow, because the two differ by
a percent and their label lengths differ by three characters.

### Layout constraints worth knowing before editing

These are recorded as comments in the generator too, because each cost a round of
rework:

- `box()` draws with `boxstyle pad=0.012`, so the visible edge is outside the rect
  the caller passes. Arrows are `zorder=2` against boxes at `3`, so an arrow
  ending on the rect has its head drawn **under** the box. `arrow()` trims both
  ends in data units to clear it; `FancyArrowPatch`'s `shrinkA`/`shrinkB` are in
  points and these panels differ 3× in physical width, so a point value that
  works in one gouges a hole in the other.
- **Above `rad` ≈ 0.87, matplotlib draws an arc in two pieces** with a visible gap
  in the middle. Measured by rasterising the patch and counting runs of coloured
  columns: 1 run at 0.86, 2 runs with a 29 px hole at 0.88. The path stays inside
  the axes, so it is a curvature artefact rather than clipping.
- `fit_text()` wraps to the box and then shrinks until matplotlib's own
  measurement of the rendered extent fits. Hand-tuned font sizes do not survive
  the first edit to a string — the first draft had five labels overflowing into
  their neighbours.
- `aspect_rad()` corrects curvature for a panel's aspect ratio (these run 1.12 and
  2.59). Useful in general, but **not** used for the feedback arc: clearing a card
  is a data-space requirement that does not shrink when curvature is corrected for
  a wide panel, so correcting it drove the arc straight through the card.

## Figure formats

Every render writes both `.svg` and `.png`. The SVG is the one to submit — it is
true vector at any scale, and output is byte-reproducible.

**Text is outlined by default**, so a figure renders identically on a machine
without the font, but labels are not selectable or editable in Illustrator or
Inkscape. Journals usually want live text. Set `CQ_SVG_EDITABLE_TEXT=1` to emit
`<text>` elements instead:

```bash
CQ_SVG_EDITABLE_TEXT=1 python scripts/figures/fig2_prescreen.py --outdir docs/figures/editable
```

Both are genuine vector output; this changes only whether glyphs are paths or
characters. Editable-text files are roughly 40-60% smaller and live in
`docs/figures/editable/`.

`save()` also strips matplotlib's SVG 1.1 DOCTYPE. It is legacy boilerplate no
browser validates against, and an XML parser configured against external
entities rejects the whole document rather than ignoring it — which is exactly
what happened when these were first published.
