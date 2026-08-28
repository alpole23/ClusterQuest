#!/usr/bin/env python3
"""Two-axis divergence view: how unusual is each BGC's scaffold vs its chemistry?

Every phosphonate BGC carries pepM (the hallmark gene) plus a coupling enzyme that sets
the downstream chemistry. Plotting divergence from characterised references on both axes
separates cases that a single tree conflates:

    typical pepM, divergent coupling   novel chemistry in a familiar scaffold
    divergent pepM, typical coupling   a distant relative doing known chemistry
    divergent in both                  unlike any characterised pathway
    typical in both                    well-characterised territory

The 60% pepM line is not arbitrary. Yu et al. (2013) correlated 342 pepM gene
neighbourhoods against PepM identity over 58,311 pairwise comparisons and found the
correlation holds only above ~60% identity, with essentially no neighbourhood
similarity below it. A BGC left of that line is not merely divergent — its pathway is
outside the range where pepM identity predicts anything about the surrounding genes.

Reads the coupling support TSV written by bgc_coupling_annotation.py --reference_faa
--reference_pepm.

    python scripts/bgc_divergence_plot.py \
        --support results/.../phosphonate_coupling_support.tsv \
        --outdir  results/.../divergence

Outputs bgc_divergence.png/.svg plus bgc_divergence_outliers.tsv.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.plotting import SVG_METADATA, canonicalise_svg          # noqa: E402
from utils.constants import COUPLING_COLORS, COUPLING_ORDER        # noqa: E402

import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.lines import Line2D                                 # noqa: E402

# Secondary encoding. The coupling palette fails a colourblind check on the
# Decarboxylase/Decarboxylase-Nucleotidyltransferase pair (deutan ΔE 3.5), and those two
# classes are chemically adjacent, so confusing them is not harmless. Shape carries the
# identity independently of hue; colour stays as-is for consistency with the GCF tree,
# heatmap and report badges, where the same class must look the same.
CLASS_MARKERS = {
    'Synthase':                             'o',
    'Reductase':                            's',
    'Decarboxylase':                        '^',
    'Decarboxylase-Nucleotidyltransferase': 'D',
    'Transaminase':                         'v',
    'Unknown':                              'X',
}

# pepM axis: Yu et al. (2013) compared 342 pepM gene neighbourhoods (6 genes either
# side of pepM, 13 total) against PepM amino-acid identity across 58,311 pairwise
# comparisons, and found a highly significant linear correlation ONLY above ~60%
# identity — below it, "essentially no similarity in the pepM gene neighborhood".
# So 60% is a floor rather than a gradient: below it, pepM identity tells you nothing
# about the surrounding pathway.
#   Yu X, Doroghazi JR, Janga SC, Zhang JK, Circello B, Griffin BM, Labeda DP,
#   Metcalf WW. Diversity and abundance of phosphonate biosynthetic genes in nature.
#   PNAS 2013;110(51):20759-20764. doi:10.1073/pnas.1315107110
PEPM_NEIGHBOURHOOD_PCT = 60.0
CITATION = ('pepM 60% line: Yu et al. 2013 PNAS 110(51):20759 '
            '(doi:10.1073/pnas.1315107110)')

# Coupling axis: no published threshold exists, so these are empirical poles measured
# on this pipeline's own Pantoea run — characterised orthologues scored 93.8-100%
# identity, unrelated members of the same superfamily 21.8-30.8%. Drawn as guides only.
BACKGROUND_PCT = 30.0
ORTHOLOGUE_PCT = 90.0

INK        = '#333333'
INK_MUTED  = '#777777'
GRID       = '#e6e6e6'
SURFACE    = '#ffffff'


def load_support(path, region_only=True):
    rows = []
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith('#')]
    for r in csv.DictReader(lines, delimiter='\t'):
        label = r['bgc']
        if region_only and '_' in label.rsplit('region', 1)[-1]:
            continue          # skip cand_cluster / protocluster sub-records
        try:
            coupling = float(r['assigned_pct_id'])
            pepm = float(r['pepm_pct_id'])
        except (KeyError, ValueError):
            continue
        if pepm <= 0:
            continue          # no pepM found — nothing to place on the x axis
        rows.append({'bgc': label, 'cls': r['assigned_class'], 'pepm': pepm,
                     'coupling': coupling, 'n_refs': r.get('assigned_n_refs', '?'),
                     'ref': r.get('assigned_ref', '-'), 'pepm_ref': r.get('pepm_ref', '-')})
    return rows


def plot(rows, outdir, taxon):
    fig, ax = plt.subplots(figsize=(9, 7.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)

    # Everything left of the 60% line is outside the range where pepM identity predicts
    # the gene neighbourhood at all — shaded because that is a statement about the
    # published correlation, not about these BGCs.
    ax.axvspan(0, PEPM_NEIGHBOURHOOD_PCT, color='#f4f4f2', zorder=0)
    ax.axvline(PEPM_NEIGHBOURHOOD_PCT, color='#b0b0b0', lw=1.3, ls='--', zorder=1)
    ax.text(PEPM_NEIGHBOURHOOD_PCT - 1.5, 50,
            'pepM 60% — below this, neighbourhood\nsimilarity is undetectable (Yu 2013)',
            rotation=90, ha='right', va='center', fontsize=8, color=INK_MUTED, style='italic')

    # Coupling-axis guides are empirical, so drawn fainter than the published line
    for v in (BACKGROUND_PCT, ORTHOLOGUE_PCT):
        ax.axhline(v, color=GRID, lw=1, zorder=0)

    # Aggregate coincident points. 236 Synthase BGCs sit at essentially one position;
    # drawn raw they overplot into a handful of dots and the cluster's weight is
    # invisible. Binning to 0.5% and sizing by sqrt(count) keeps every BGC represented
    # and makes the dominant cluster read as dominant.
    present = [c for c in COUPLING_ORDER if any(r['cls'] == c for r in rows)]
    for cls in present:
        pts = [r for r in rows if r['cls'] == cls]
        binned = {}
        for r in pts:
            key = (round(r['pepm'] * 2) / 2, round(r['coupling'] * 2) / 2)
            binned[key] = binned.get(key, 0) + 1
        xs = [k[0] for k in binned]
        ys = [k[1] for k in binned]
        sizes = [40 + 26 * (n ** 0.5) for n in binned.values()]
        ax.scatter(xs, ys, s=sizes, c=COUPLING_COLORS.get(cls, '#999999'),
                   marker=CLASS_MARKERS.get(cls, 'o'),
                   edgecolors=SURFACE, linewidths=1.2,   # 2px-equivalent surface ring
                   alpha=0.85, zorder=3, label=f'{cls} (n={len(pts)})')

    ax.set_xlabel('pepM identity to nearest characterised reference (%)',
                  fontsize=11, color=INK)
    ax.set_ylabel('Coupling enzyme identity to nearest reference of its class (%)',
                  fontsize=11, color=INK)
    ax.set_title(f'BGC divergence from characterised phosphonate chemistry — {taxon}',
                 fontsize=13, color=INK, pad=14)
    ax.set_xlim(0, 105); ax.set_ylim(0, 105)
    ax.tick_params(colors=INK_MUTED, labelsize=9)
    for sp in ax.spines.values():
        sp.set_color(GRID)

    # Quadrant captions in axes coordinates, pinned to the corners so they cannot
    # collide with marks wherever the data happens to fall.
    cap = dict(fontsize=8.5, color=INK_MUTED, style='italic',
               transform=ax.transAxes, zorder=2)
    ax.text(0.02, 0.97, 'scaffold beyond the reach of\nneighbourhood prediction',
            ha='left', va='top', **cap)
    ax.text(0.98, 0.03, 'familiar scaffold\nnovel chemistry', ha='right', va='bottom', **cap)
    ax.text(0.02, 0.03, 'unlike any characterised\npathway on either axis',
            ha='left', va='bottom', **cap)
    ax.text(0.5, -0.105, 'marker area \u221d number of BGCs at that position   \u00b7   ' + CITATION,
            ha='center', va='top', transform=ax.transAxes,
            fontsize=8, color=INK_MUTED, style='italic')

    leg = ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False,
                    fontsize=9, labelcolor=INK, title='Coupling enzyme class')
    leg.get_title().set_color(INK)
    leg.get_title().set_fontsize(9.5)

    fig.tight_layout()
    png = os.path.join(outdir, 'bgc_divergence.png')
    svg = os.path.join(outdir, 'bgc_divergence.svg')
    fig.savefig(png, dpi=200, bbox_inches='tight', facecolor=SURFACE)
    fig.savefig(svg, bbox_inches='tight', facecolor=SURFACE, metadata=SVG_METADATA)
    canonicalise_svg(svg)
    plt.close(fig)
    return png, svg


def write_outliers(rows, outdir):
    """Review shortlist, ordered by the published pepM criterion first.

    `below_pepm_60` is the stronger signal: Yu et al. (2013) found pepM identity
    predicts the gene neighbourhood only above ~60%, so a BGC beneath that line is
    outside the range where anything can be inferred about its pathway from pepM.
    `below_coupling_background` is this pipeline's own empirical guide and is weaker —
    several classes have a single reference, from a different genus.
    """
    out = os.path.join(outdir, 'bgc_divergence_outliers.tsv')
    flagged = [r for r in rows
               if r['pepm'] < PEPM_NEIGHBOURHOOD_PCT or r['coupling'] < BACKGROUND_PCT]
    flagged.sort(key=lambda r: (r['pepm'], r['coupling']))
    with open(out, 'w') as f:
        f.write('# Review shortlist for BGCs distant from characterised phosphonate chemistry.\n')
        f.write('# below_pepm_60: pepM identity < 60%% to the nearest characterised reference.\n')
        f.write('#   Yu et al. 2013 PNAS 110(51):20759 (doi:10.1073/pnas.1315107110) correlated\n')
        f.write('#   342 pepM gene neighbourhoods against PepM identity over 58,311 pairwise\n')
        f.write('#   comparisons: the correlation holds only above ~60%%, with essentially no\n')
        f.write('#   neighbourhood similarity below it. This is the stronger criterion.\n')
        f.write('# below_coupling_background: coupling enzyme < %.0f%% identity — an empirical\n' % BACKGROUND_PCT)
        f.write('#   guide from this pipeline, weakened by classes with a single reference.\n')
        f.write('# Ambiguous by construction: distant may mean wrong class, or novel variant.\n')
        f.write('# A shortlist for manual review, not a verdict.\n')
        f.write('bgc\tassigned_class\tpepm_pct_id\tpepm_ref\tbelow_pepm_60\t'
                'coupling_pct_id\tcoupling_ref\tn_refs\tbelow_coupling_background\n')
        for r in flagged:
            f.write(f"{r['bgc']}\t{r['cls']}\t{r['pepm']:.1f}\t{r['pepm_ref']}\t"
                    f"{'yes' if r['pepm'] < PEPM_NEIGHBOURHOOD_PCT else 'no'}\t"
                    f"{r['coupling']:.1f}\t{r['ref']}\t{r['n_refs']}\t"
                    f"{'yes' if r['coupling'] < BACKGROUND_PCT else 'no'}\n")
    return out, len(flagged)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--support', required=True, help='phosphonate_coupling_support.tsv')
    ap.add_argument('--outdir', default='.', help='Output directory')
    ap.add_argument('--taxon', default='', help='Taxon name for the title')
    ap.add_argument('--all_records', action='store_true',
                    help='Include cand_cluster/protocluster sub-records (default: regions only)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rows = load_support(args.support, region_only=not args.all_records)
    if not rows:
        print('No BGCs with both pepM and coupling support — nothing to plot.', file=sys.stderr)
        return 1
    print(f'  {len(rows)} BGCs with both axes')
    png, svg = plot(rows, args.outdir, args.taxon or 'phosphonate BGCs')
    print(f'  wrote {png}')
    print(f'  wrote {svg}')
    tsv, n = write_outliers(rows, args.outdir)
    print(f'  wrote {tsv}  ({n} BGCs below background)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
