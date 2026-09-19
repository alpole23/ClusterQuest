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
    'transferred': '#8c6d4f',
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
    'transferred': 'named by GCF annotation transfer \u2020',
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
    ('leud', 'LeuD  (small subunit)'),
    ('isopropylmalate dehydratase', 'LeuC  (large subunit)'),
    ('phytanoyl', 'phytanoyl-CoA dioxygenase'),
    ('aspartate aminotransferase', 'Asp aminotransferase'),
    ('homoaconitate', 'homoaconitate synthase'),
    ('monooxygenase', 'monooxygenase'),
    ('atp-grasp', 'ATP-grasp'),
    ('fad/nad', 'FAD/NAD(P)-binding'),
    ('carbamoyl-phosphate synthase', 'carbamoyl-P synthase (ATP-grasp)'),
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


def load_transferred(path):
    """{locus_tag: transferred_product} from GCF_ANNOTATION_TRANSFER, if present.

    A recovered ORF usually reaches the region GenBank as "hypothetical" --
    prodigal finds the reading frame, it does not name the protein. The name
    arrives later, from annotated relatives in the same GCF. Reading only the
    GenBank therefore draws the most interesting recovered genes as
    unclassified: LMG 5342's recovered_HE617160.1_0066 is a GNAT family
    N-acetyltransferase on five source genomes, and the first version of this
    figure showed it as a grey unknown.
    """
    import csv as _csv
    out = {}
    if not path or not Path(path).exists():
        return out
    with open(path) as fh:
        for row in _csv.DictReader(fh, delimiter='\t'):
            prod = (row.get('transferred_product') or '').strip()
            if prod and row.get('origin') == 'transferred':
                out[row['locus_tag']] = prod
    return out


TRANSFER_TSV = (ROOT / 'results/main_analysis_results/Erwiniaceae/'
                'annotation_transfer/gcf_annotation_transfer.tsv')


def read_region(path, transferred=None):
    transferred = transferred if transferred is not None else {}
    rec = next(SeqIO.parse(str(path), 'genbank'))
    genes = []
    for f in rec.features:
        if f.type != 'CDS':
            continue
        tag = f.qualifiers.get('locus_tag', [''])[0]
        name = transferred.get(tag)
        label = label_for(f)
        cat = gene_category(f)
        if name:
            # A transferred name is an inference from homologues, not an
            # observation, so it is marked on the figure rather than passed off
            # as antiSMASH's own call.
            label = f'{name[:34]}  †'
            if cat == 'other':
                cat = 'transferred'
        genes.append({
            'start': int(f.location.start),
            'end': int(f.location.end),
            'strand': 1 if f.location.strand in (None, 1) else -1,
            'recovered': is_recovered(f),
            'category': cat,
            'label': label,
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


def draw_clade(ax, title, before_path, after_path, transferred=None):
    span_b, before = read_region(before_path, transferred)
    span_a, after = read_region(after_path, transferred)
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


def gain_distribution(before_dir, after_dir):
    """Per-region gene gain across a whole clade, so one example can be placed.

    Without this the figure invites the reader to take the illustrated cluster
    as typical. LMG 5342 gained 14 genes and is rank 1 of 334; the median gain
    among regions that gained anything at all is 2, and 65% gained nothing.
    """
    before_dir, after_dir = Path(before_dir), Path(after_dir)
    gains = {}
    for bg in before_dir.glob('*/*region*.gbk'):
        ag = after_dir / bg.parent.name / bg.name
        if not ag.exists():
            continue
        nb = sum(1 for f in next(SeqIO.parse(str(bg), 'genbank')).features
                 if f.type == 'CDS')
        na = sum(1 for f in next(SeqIO.parse(str(ag), 'genbank')).features
                 if f.type == 'CDS')
        gains[f'{bg.parent.name}/{bg.name}'] = na - nb
    return gains


def panel_distribution(ax, clade_gains, highlights):
    """Per-region gains, one row per SOURCE RUN, with every example marked.

    Rows are keyed on the run, not on the panel: three of the panels above come
    from the same Erwiniaceae before/after pair, and giving each its own row
    drew one distribution three times over and printed "117/333 gained" thrice
    as though they were independent measurements.
    """
    ytick, ylab = [], []
    for row, (clade, gains, marks) in enumerate(clade_gains):
        vals = sorted(gains.values())
        y = len(clade_gains) - row - 1
        ytick.append(y)
        n_gained = sum(1 for v in vals if v > 0)
        ylab.append(f'{clade}\n{n_gained}/{len(vals)} gained')
        # Deterministic beeswarm: spread each tied group about its own row, so
        # the height of a column reads as its count. Cycling i % 5 instead
        # stacked every value into the same five rows and looked like dashes.
        from collections import Counter
        seen, counts = Counter(), Counter(vals)
        ys = []
        for v in vals:
            k, n = seen[v], counts[v]
            seen[v] += 1
            offset = 0.0 if n == 1 else (k / (n - 1) - 0.5) * min(0.62, 0.05 * n)
            ys.append(y + offset)
        ax.scatter(vals, ys, s=11, color='#c3ccd6', edgecolor='none',
                   zorder=2, alpha=0.85)
        # Several panels can share one run, so a row carries several markers.
        for i, (letter, key) in enumerate(marks):
            if key not in gains:
                continue
            ax.scatter([gains[key]], [y], s=95, color=ACCENT, zorder=5,
                       edgecolor='white', lw=1.2)
            ax.annotate(f'{letter}  +{gains[key]}', xy=(gains[key], y),
                        xytext=(gains[key], y + (0.30 if i % 2 == 0 else -0.42)),
                        ha='center', fontsize=7.6,
                        color=ACCENT, fontweight='bold')
    ax.set_yticks(ytick)
    ax.set_yticklabels(ylab, fontsize=8)
    ax.set_xlabel('genes gained by ORF recovery, per BGC region')
    ax.set_ylim(-0.55, len(clade_gains) - 0.25)
    ax.axvline(0, color='#c9ced6', lw=0.9, zorder=1)
    ax.set_title('Where each example sits in its clade', loc='left',
                 fontsize=9.5, pad=6)
    for side in ('left', 'top', 'right'):
        ax.spines[side].set_visible(False)


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
    ap.add_argument('--no-distribution', action='store_true',
                    help='omit the panel placing each example within its clade')
    args = ap.parse_args()

    base = ROOT / 'results' / 'antismash_results'
    # (short clade name, panel title, before dir, after dir, region key)
    #
    # Ordered as a gradient in how much the deposit left out, because that is
    # the argument: recovery is targeted, not indiscriminate.
    CLADES = [
        ('Winslowiella',
         'Winslowiella iniecta B149  —  recovery restores CORE and '
         'TAILORING enzymes; GCF 11 → 9',
         base / 'Erwiniaceae_pre_recovery', base / 'Erwiniaceae',
         'Winslowiella_iniecta_B149/JRXF01000012.1.region001.gbk'),
        ('Pantoea',
         'Pantoea ananatis LMG 5342 region 2  —  confirmed '
         'phosphonolipid; most affected BGC of 334',
         base / 'Erwiniaceae_pre_recovery', base / 'Erwiniaceae',
         'Pantoea_ananatis_LMG_5342/HE617160.1.region002.gbk'),
        # Same genome, same deposit, same year as the panel above. One region
        # gains 14 genes and this one gains 4, so the difference is which genes
        # NCBI's pipeline happened to call rather than anything about the
        # assembly -- an internal control no separate strain can provide.
        # NOT a tailoring example: its monooxygenase, homoaconitate synthase,
        # isopropylmalate dehydratase, methyltransferase and SanS are all
        # DEPOSITED. All four recovered genes are unclassified by antiSMASH.
        ('pantaphos',
         'Pantoea ananatis LMG 5342 region 1  ·  pantaphos / HiVir  '
         '—  same genome, already well annotated; recovery completes '
         'the LeuC/LeuD dehydratase',
         base / 'Erwiniaceae_pre_recovery', base / 'Erwiniaceae',
         'Pantoea_ananatis_LMG_5342/HE617160.1.region001.gbk'),
        # S. griseus, not S. hygroscopicus: the bialaphos lineage yields only 2
        # phosphonate regions across 39 genomes (reproducing its count of 2 in
        # the actinomycete comparison), too few to choose an example from.
        ('Streptomyces',
         'Streptomyces griseus NRRL B-2929  \u2014  recovery restores two '
         'ATP-grasp tailoring enzymes',
         ROOT / 'results_fig1_gris_norecover/antismash_results/Streptomyces_griseus',
         ROOT / 'results_fig1_gris_recover/antismash_results/Streptomyces_griseus',
         None),
        ('Bacteroides',
         'Bacteroides fragilis BFG-525  ·  region 1  —  deposit '
         'already complete; recovery adds hypotheticals only',
         ROOT / 'results_fig1_bact_norecover/antismash_results/Bacteroides_fragilis',
         ROOT / 'results_heldout_off/antismash_results/Bacteroides_fragilis',
         'Bacteroides_fragilis_BFG-525/CP103089.1.region001.gbk'),
    ]

    # Rows of the distribution panel are keyed on the SOURCE RUN, and several
    # panels may share one; the Erwiniaceae pair supplies three of them.
    ROW_NAME = {'Winslowiella': 'Erwiniaceae', 'Pantoea': 'Erwiniaceae',
                'pantaphos': 'Erwiniaceae'}

    transferred = load_transferred(TRANSFER_TSV)
    print(f'{len(transferred)} transferred gene names available')
    clades, rows, cache = [], {}, {}
    for short, title, bdir, adir, key in CLADES:
        if not (bdir.exists() and adir.exists()):
            print(f'skipping {short}: {bdir if not bdir.exists() else adir} missing')
            continue
        run = (str(bdir), str(adir))
        if run not in cache:
            cache[run] = gain_distribution(bdir, adir)
        gains = cache[run]
        if not gains:
            print(f'skipping {short}: no paired regions')
            continue
        if key is None:                       # pick the clade's best example
            key = max(gains, key=lambda k: gains[k])
            title = f'{title}  —  best of {len(gains)} regions'
        letter = 'ABCDEF'[len(clades)]
        rows.setdefault(ROW_NAME.get(short, short), (gains, []))[1].append((letter, key))
        clades.append((title, bdir / key, adir / key))
    clade_gains = [(name, gains, marks) for name, (gains, marks) in rows.items()]

    nrows = len(clades) + (0 if args.no_distribution else 1)
    fig, axes = plt.subplots(nrows, 1, figsize=(11.6, 2.15 * len(clades) + 2.6),
                             squeeze=False,
                             gridspec_kw={'height_ratios':
                                          [2.15] * len(clades) +
                                          ([1.7] if not args.no_distribution else [])})
    used = []
    for ax, (name, before, after) in zip(axes[:, 0], clades):
        nb, na = draw_clade(ax, name, before, after, transferred)
        print(f'{name}: {nb} -> {na} genes')
        for g in read_region(after, transferred)[1]:
            if g['category'] not in used:
                used.append(g['category'])
    order = [c for c in CATEGORY_COLOR if c in used]

    if not args.no_distribution:
        panel_distribution(axes[len(clades), 0], clade_gains, {})

    for ax, letter in zip(axes[:, 0], 'ABCDEF'):
        panel_label(ax, letter, dx=-0.175, dy=1.22)

    axes[-1, 0].legend(handles=legend_handles(order), loc='upper center',
                       bbox_to_anchor=(0.5, -0.38), ncol=4, fontsize=8)
    fig.suptitle('Gene calling recovers the genes that say what a cluster makes '
                 '— where the deposit left them out',
                 fontsize=11.5, fontweight='bold', y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, args.outdir, 'fig1_orf_recovery')


if __name__ == '__main__':
    main()
