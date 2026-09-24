# Committed inputs for Figure 1 and the window × recovery figure

The manuscript says every figure regenerates from committed data. That was not
true: `fig1_orf_recovery.py` and `fig1_window_x_recovery.py` read region
GenBanks out of six A/B run directories under `results_*/`, which is gitignored
and was 177 GB. Deleting those directories would have made both figures
unreproducible, silently — the scripts skipped missing clades with a printed
note and drew whatever was left.

This directory is what they need instead, at 2.5 MB.

| | what | why only this |
|---|---|---|
| `<clade>_before/`, `<clade>_after/` | the region GenBanks the panels **draw** | 14 files. Everything else in those runs was only ever counted. |
| `gene_gains.tsv` | per-region CDS gain for every region of every clade | the distribution panel needs 717 numbers, not 717 GenBanks |
| `gcf_annotation_transfer.tsv` | transferred product names | labels recovered genes that reach the GenBank as "hypothetical" |

## Why the gains are a table and not recomputed

`gain_distribution()` reads `gene_gains.tsv` and never globs. That is deliberate:
only the drawn regions are committed here, so a glob over this directory would
find a handful of GenBanks, return a **partial** distribution and draw a wrong
panel with no error. One source of truth removes that failure mode.

Both scripts still prefer a live run directory when one is present
(`src_dir()`), so re-running the pipeline and regenerating from real output
stays possible. To refresh the table after such a run:

```bash
python scripts/figures/fig1_orf_recovery.py --rebuild-gains
```

## Verified

With the six source directories moved aside, both figures regenerate
**byte-identically** to the versions produced with them present. That check is
the only reason the directories were deleted.

Separately, the committed SVGs were found to be stale — produced by earlier
versions of the scripts and never refreshed. Output is deterministic across
processes (checked), so the drift was staleness, not non-determinism. All four
figures were regenerated at the same commit that added this directory.
