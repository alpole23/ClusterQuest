#!/usr/bin/env python3
"""Figure 1 — what ORF recovery restores to a phosphonate BGC.

antiSMASH runs gene finding only on records with ZERO CDS features, so a
GenBank deposit that annotates *some* of its genes is trusted for all of them.
On *P. ananatis* LMG 5342 that meant 14 of 28 genes were invisible -- including
the AEP transaminase and both CDP-alcohol phosphatidyltransferases, which are
the genes that say what the cluster makes. Every gene-content metric built on
top was therefore measuring annotation quality rather than biology.

Genes are coloured by Pfam category via utils/domain_functions, NOT by product
keyword: a keyword metric was tried first and inverted on the two lab-confirmed
clusters, because "serine hydroxymethyltransferase" matches "methyltransferase"
and "aspartate-semialdehyde dehydrogenase" matches "dehydrogenase". Recovered
genes carry a heavy outline so the before/after difference is readable even in
greyscale.

Usage:
    python scripts/figures/fig1_orf_recovery.py --outdir docs/figures
"""
import argparse
import sys
from pathlib import Path

from figure_style import ACCENT, FAINT, GOOD, INK, panel_label, plt, save

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))

from utils import domain_functions as df  # noqa: E402

from Bio import SeqIO  # noqa: E402

CATEGORY_COLOR = {
    'core': '#b3392b',
    'tailoring': '#2f7d5d',
    'lipid': '#7a5ea8',
    'transport': '#1b5e7e',
    'regulation': '#c98a1b',
    'mobile': '#8a939f',
    'primary': '#cfd5dd',
    'other': '#e4e8ed',
}
CATEGORY_LABEL = {
    'core': 'phosphonate core',
    'tailoring': 'tailoring',
    'lipid': 'lipid handling',
    'transport': 'transport',
    'regulation': 'regulation',
    'mobile': 'mobile element',
    'primary': 'primary metabolism',
    'other': 'unclassified',
}

# Products worth naming on the figure. Matched as a lowercase substring against
# the product text, and only ever used to place a LABEL -- never to classify,
# which is what domain_functions is for.
LABEL_HINTS = [
    ('phosphoenolpyruvate phosphomutase', 'pepM'),
    ('phosphonopyruvate decarboxylase', 'Ppd'),
    ('2-aminoethylphosphonate', 'AEP transaminase'),
    ('cdp-alcohol', 'CDP-alcohol PT'),
    ('phosphocholine', 'phosphocholine CT'),
    ('mfs', 'MFS transporter'),
    ('ntp', 'NTP transferase'),
]


def is_recovered(feat):
    """Prodigal-recovered CDS are tagged recovered_<record>_<n> by recover_orfs.py."""
    return feat.qualifiers.get('locus_tag', [''])[0].startswith('recovered_')


def gene_category(feat):
    """antiSMASH's own functional call for a CDS.

    Region GenBanks carry NO Pfam accessions -- checked, 0 of 28 on this region
    -- so utils/domain_functions has nothing to read here and `gene_kind` is the
    authoritative signal antiSMASH does write. Product keywords are deliberately
    not used: that metric was tried and inverted on the two lab-confirmed
    clusters, because "serine hydroxymethyltransferase" contains
    "methyltransferase".

    `(none)` is the honest answer for 15 of 28 genes, and it is part of the
    point: the CDP-alcohol phosphatidyltransferases that decide what this
    cluster makes are among the genes antiSMASH does not classify.
    """
    kind = feat.qualifiers.get('gene_kind', [''])[0]
    return {
        'biosynthetic': 'core',
        'biosynthetic-additional': 'tailoring',
        'transport': 'transport',
        'regulatory': 'regulation',
    }.get(kind, 'other')


def label_for(feat):
    prod = feat.qualifiers.get('product', [''])[0].lower()
    for needle, short in LABEL_HINTS:
        if needle in prod:
            return short
    return None


def read_region(path):
    rec = next(SeqIO.parse(str(path), 'genbank'))
    genes = []
    for f in rec.features:
        if f.type != 'CDS':
            continue
        genes.append({
            'start': int(f.location.start),
            'end': int(f.location.end),
            'strand': 1 if f.location.strand in (None, 1) else -1,
            'recovered': is_recovered(f),
            'category': gene_category(f),
            'label': label_for(f),
        })
    return len(rec.seq), sorted(genes, key=lambda g: g['start'])


def draw_track(ax, genes, span, y, height=0.30, head=0.30):
    """One row of gene arrows on a shared coordinate axis."""
    ax.plot([0, span], [y, y], color='#c9ced6', lw=1.0, zorder=1)
    for g in genes:
        x0, x1 = g['start'], g['end']
        width = x1 - x0
        head_len = min(width * head, span * 0.010)
        body = width - head_len
        col = CATEGORY_COLOR[g['category']]
        edge = INK if g['recovered'] else '#ffffff'
        lw = 1.5 if g['recovered'] else 0.5
        if g['strand'] >= 0:
            verts = [(x0, y - height / 2), (x0 + body, y - height / 2),
                     (x1, y), (x0 + body, y + height / 2), (x0, y + height / 2)]
        else:
            verts = [(x1, y - height / 2), (x1 - body, y - height / 2),
                     (x0, y), (x1 - body, y + height / 2), (x1, y + height / 2)]
        ax.add_patch(plt.Polygon(verts, closed=True, facecolor=col,
                                 edgecolor=edge, lw=lw, zorder=3))


def annotate_labels(ax, genes, y, span, above=True):
    """Name the genes that carry the argument, staggered so they do not collide."""
    labelled = sorted((g for g in genes if g['label']),
                      key=lambda g: g['start'])
    # Three tiers, and each label nudged off its gene only as far as needed:
    # at two tiers "MFS transporter" and "NTP transferase" overlapped outright.
    for i, g in enumerate(labelled):
        mid = (g['start'] + g['end']) / 2
        tier = i % 3
        dy = (0.40 + tier * 0.30) * (1 if above else -1)
        ax.annotate(g['label'], xy=(mid, y + (0.16 if above else -0.16)),
                    xytext=(mid, y + dy),
                    ha='center', va='bottom' if above else 'top',
                    fontsize=7.0,
                    color=ACCENT if g['recovered'] else FAINT,
                    fontweight='bold' if g['recovered'] else 'normal',
                    arrowprops=dict(arrowstyle='-', lw=0.6,
                                    color=ACCENT if g['recovered'] else '#c9ced6'))


def draw_clade(ax, title, before_path, after_path):
    span_b, before = read_region(before_path)
    span_a, after = read_region(after_path)
    span = max(span_b, span_a)

    draw_track(ax, before, span_b, y=1.0)
    draw_track(ax, after, span_a, y=0.0)
    annotate_labels(ax, after, y=0.0, span=span, above=False)

    n_rec = sum(1 for g in after if g['recovered'])
    # len(after) - len(before) can exceed n_rec: recovered genes lengthen the
    # cluster, so the region boundary moves and can sweep in a gene that was
    # always annotated. Report both rather than let one stand for the other.
    swept = (len(after) - len(before)) - n_rec
    ax.text(-span * 0.015, 1.0, f'deposited\n{len(before)} genes',
            ha='right', va='center', fontsize=8, color=FAINT)
    ax.text(-span * 0.015, 0.0, f'after recovery\n{len(after)} genes',
            ha='right', va='center', fontsize=8, color=INK, fontweight='bold')
    tail = f'\n+{swept} by wider\nboundary' if swept else ''
    ax.text(span * 1.02, 0.5, f'+{n_rec} genes\nrecovered{tail}',
            ha='left', va='center', fontsize=8.5, color=GOOD, fontweight='bold')

    ax.set_xlim(-span * 0.20, span * 1.22)
    ax.set_ylim(-1.55, 1.75)
    ax.set_yticks([])
    ax.set_xticks([0, span / 2, span])
    ax.set_xticklabels(['0', f'{span / 2000:.0f} kb', f'{span / 1000:.0f} kb'],
                       fontsize=7.5)
    for side in ('left', 'top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_title(title, loc='left', fontsize=9.5, pad=6)
    return len(before), len(after)


def legend_handles(used):
    import matplotlib.patches as mpatches
    h = [mpatches.Patch(facecolor=CATEGORY_COLOR[c], edgecolor='white',
                        label=CATEGORY_LABEL[c]) for c in used]
    h.append(mpatches.Patch(facecolor='#ffffff', edgecolor=INK, lw=1.5,
                            label='recovered by gene calling'))
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='docs/figures')
    ap.add_argument('--clade', action='append', default=[], metavar='NAME:BEFORE:AFTER',
                    help='repeatable; defaults to the P. ananatis example alone')
    args = ap.parse_args()

    clades = []
    if args.clade:
        for spec in args.clade:
            name, before, after = spec.split(':', 2)
            clades.append((name, Path(before), Path(after)))
    else:
        base = ROOT / 'results' / 'antismash_results'
        clades = [(
            'Pantoea ananatis LMG 5342  ·  HE617160.1 region 2  '
            '(confirmed phosphonolipid)',
            base / 'Erwiniaceae_pre_recovery/Pantoea_ananatis_LMG_5342/HE617160.1.region002.gbk',
            base / 'Erwiniaceae/Pantoea_ananatis_LMG_5342/HE617160.1.region002.gbk',
        )]

    fig, axes = plt.subplots(len(clades), 1,
                             figsize=(11.6, 2.5 * len(clades) + 1.2),
                             squeeze=False)
    used = []
    for ax, (name, before, after) in zip(axes[:, 0], clades):
        nb, na = draw_clade(ax, name, before, after)
        print(f'{name}: {nb} -> {na} genes')
        for g in read_region(after)[1]:
            if g['category'] not in used:
                used.append(g['category'])
    order = [c for c in CATEGORY_COLOR if c in used]

    for ax, letter in zip(axes[:, 0], 'ABC'):
        panel_label(ax, letter, dx=-0.175, dy=1.22)

    axes[-1, 0].legend(handles=legend_handles(order), loc='upper center',
                       bbox_to_anchor=(0.5, -0.22), ncol=4, fontsize=8)
    fig.suptitle('Gene calling recovers the genes that say what a cluster makes',
                 fontsize=11.5, fontweight='bold', y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, args.outdir, 'fig1_orf_recovery')


if __name__ == '__main__':
    main()
