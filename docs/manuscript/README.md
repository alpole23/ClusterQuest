# Manuscript source

`index.html` is the working draft, published privately as a Claude Artifact. This copy is
the version-controlled one: the artifact has its own history, but nothing there is
reviewable in a diff alongside the code and data that produced the numbers.

## The figures are NOT in this directory

They are `docs/figures/fig1_window_x_recovery.svg`, `fig1_orf_recovery.svg`,
`fig2_prescreen.svg` and `fig3_partitioning.svg`, referenced by the published artifact as
`fig1.svg` … `fig4.svg` in that order. Manuscript numbering and script numbering do not
line up — Figure 2 is `fig1_orf_recovery`, Figure 3 is `fig2_prescreen` — which is a trap
worth knowing before editing either.

**They are published separately from the page, and that has already gone wrong once.**
Republishing `index.html` alone leaves the artifact's figures untouched, so the text can
say 56% and 32× while the embedded figures show neither. It did, between versions 6 and 8.
When a figure changes, republish with an explicit file map:

```
files = {"fig1.svg": "docs/figures/fig1_window_x_recovery.svg",
         "fig2.svg": "docs/figures/fig1_orf_recovery.svg",
         "fig3.svg": "docs/figures/fig2_prescreen.svg",
         "fig4.svg": "docs/figures/fig3_partitioning.svg"}
```

A quick check: compare byte sizes of the artifact's files against `docs/figures/`. They
should match exactly, since SVG output is byte-deterministic.

## Regenerating everything

```bash
for f in fig1_window_x_recovery fig1_orf_recovery fig2_prescreen fig3_partitioning; do
    python scripts/figures/$f.py --outdir docs/figures
done
```

No run directory needed — inputs are committed under `docs/comparisons/figure_inputs/`.
See the "Paper figures" section of CLAUDE.md.

## Still to complete

The author list, contributions, funding and reference list are all marked `to be
completed` in the draft.
