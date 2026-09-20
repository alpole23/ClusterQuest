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
*Erwiniaceae* saving is quoted as 6.5×, which is the antiSMASH-only figure
(3,680 → 569 CPU-min); on whole-run totals it is 3.0×, and this figure uses
totals throughout so the four clades are comparable.

The *B. fragilis* bar is the one place neutrality failed: the screen returned
141 of 143 regions because a pseudogene-flagged pepM carries no `/translation`.
After the fix it returns 143. The corrected bar is **not drawn** — 143 against
141 is 1.4% of that axis and would be an invisible mark pretending to be
evidence — so it is stated in words instead.

## Figure 3 — BiG-SCAPE partitioning

The honest framing, and the one the data supports: partitioning is **1.7–3.7×
slower** at every scale measured, because each partition re-pays BiG-SCAPE's
fixed Pfam-load cost. What it buys is memory — ~1,890 GB in one job at the
121,000 BGCs a million genomes yields, against ~84 GB for the largest partition
— and it buys that without changing a single cluster assignment (ARI 1.0000,
identical family counts, 0 split and 0 merged).

Panel B draws the memory fit **solid over its measured range and dashed where it
is extrapolated**. The quadratic was fitted on 1,500–4,000 BGCs; below ~4,000
peak RSS looks flat and the fit does not describe it, and 121,000 is a 12×
reach. The figure must not read as 12× more measurement than exists.

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
