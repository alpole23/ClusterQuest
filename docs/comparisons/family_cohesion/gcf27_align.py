# Recorded experiment, not a pipeline stage. Run from the repo root.
"""GCF-27 member regions aligned on pepM, to show what contig truncation does.

Every member is anchored at its pepM start and flipped so pepM reads left-to-right,
because the question is what sits AROUND the anchor and strand is an artefact of
which way the contig was assembled.
"""
import csv, collections, sys
from pathlib import Path
from Bio import SeqIO
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch
sys.path.insert(0, 'scripts')
from utils.plotting import SVG_METADATA, canonicalise_svg

W = 'work/55/25510a8b8ea5d6670b01b9cc270f18/'
R = Path('/home/alexp/pipeline/work/7b/ce057e0769639e3546bb333b2cafd0/antismash_input')
pc = [r for r in csv.DictReader(open(W+'gcf_annotation_transfer.tsv'), delimiter='\t')
      if r['family'] == '27']
grp = {(r['region'], r['locus_tag']): r['group'] for r in pc}
cons = {r['group']: r for r in csv.DictReader(open(W+'gcf_consensus_clusters.tsv'), delimiter='\t')
        if r['family'] == '27'}
regions = sorted({(r['genome'], r['region']) for r in pc})

rows = []
for genome, region in regions:
    p = R/genome/f'{region}.gbk'
    if not p.exists(): continue
    recs = list(SeqIO.parse(str(p), 'genbank'))
    pf = {(f.qualifiers.get('locus_tag') or ['?'])[0]
          for rec in recs for f in rec.features if f.type == 'PFAM_domain'
          and any(x.split('.')[0] == 'PF01066' for x in f.qualifiers.get('db_xref', []))}
    edge = any('True' in str(f.qualifiers.get('contig_edge', ''))
               for rec in recs for f in rec.features if f.type == 'region')
    genes, anchor, L = [], None, 0
    for rec in recs:
        L = max(L, len(rec.seq))
        for f in rec.features:
            if f.type != 'CDS': continue
            tag = (f.qualifiers.get('locus_tag') or ['?'])[0]
            prod = (f.qualifiers.get('product') or [''])[0]
            s, e = int(f.location.start), int(f.location.end)
            st = 1 if f.location.strand in (None, 1) else -1
            g = {'s': s, 'e': e, 'st': st, 'tag': tag, 'prod': prod,
                 'grp': grp.get((region, tag)), 'pf': tag in pf}
            if 'phosphoenolpyruvate mutase' in prod.lower():
                anchor = g
            genes.append(g)
    if anchor is None: continue
    flip = anchor['st'] == -1
    off = anchor['s'] if not flip else anchor['e']
    for g in genes:
        if flip:
            g['x0'], g['x1'], g['st'] = off - g['e'], off - g['s'], -g['st']
        else:
            g['x0'], g['x1'] = g['s'] - off, g['e'] - off
    lo = (off - L) if flip else (0 - off)
    hi = off if flip else (L - off)
    rows.append({'genome': genome, 'genes': genes, 'edge': edge,
                 'lo': lo, 'hi': hi, 'n': len(genes), 'has_pf': bool(pf)})

rows.sort(key=lambda r: (not r['has_pf'], -(r['hi']-r['lo'])))

# colour the groups that are widespread enough to show synteny
freq = collections.Counter(g['grp'] for r in rows for g in r['genes'] if g['grp'])
PALETTE = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b',
           '#e377c2','#7f7f7f','#bcbd22','#17becf','#aec7e8','#ffbb78',
           '#98df8a','#ff9896','#c5b0d5']
# GCF-27 fragments into 505 groups across 38 members, so almost nothing clears a
# 50% threshold. Colour the most frequent groups instead -- the point is to show
# which genes are shared and where they sit, not to assert a core.
top = [g for g, c in freq.most_common() if c >= 4][:len(PALETTE)]
cmap = dict(zip(top, PALETTE))
PFC = '#111111'

fig, ax = plt.subplots(figsize=(15, 0.33*len(rows) + 2.1))
for i, r in enumerate(rows):
    y = len(rows) - i
    ax.plot([r['lo']/1000, r['hi']/1000], [y, y], color='#d8dce0', lw=1.1, zorder=1)
    for g in r['genes']:
        c = PFC if g['pf'] else cmap.get(g['grp'], '#e8eaec')
        w = (g['x1']-g['x0'])/1000
        ax.add_patch(FancyArrow(
            g['x0']/1000 if g['st'] == 1 else g['x1']/1000, y,
            w*g['st'], 0, width=0.44, length_includes_head=True,
            head_length=min(0.55, w*0.5), head_width=0.62,
            facecolor=c, linewidth=0.35, edgecolor='#4a4a4a', zorder=3))
    if r['edge']:
        for xend in (r['lo']/1000, r['hi']/1000):
            ax.plot([xend, xend], [y-0.46, y+0.46], color='#d62728', lw=2.4, zorder=5)
    ax.text(-0.5 + min(rr['lo'] for rr in rows)/1000, y,
            f"{r['genome'][:34]}  ({r['n']})", ha='right', va='center', fontsize=6.6,
            color='#c0392b' if r['edge'] else '#222')
ax.axvline(0, color='#2c5aa0', lw=1.0, ls='--', zorder=6)
ax.set_xlabel('kb relative to pepM start (every region flipped so pepM reads left-to-right)')
ax.set_yticks([]); ax.set_ylim(0.2, len(rows)+0.9)
ax.spines[['top','right','left']].set_visible(False)
ax.legend(handles=[Patch(color=PFC, label='CDP-alcohol phosphatidyltransferase (PF01066)'),
                   Patch(color='#e8eaec', label='gene in <4 members'),
                   plt.Line2D([],[], color='#d62728', lw=2.4, label='contig edge (region truncated)'),
                   plt.Line2D([],[], color='#2c5aa0', lw=1.0, ls='--', label='pepM anchor')],
          loc='upper left', bbox_to_anchor=(0, -0.055/ (0.33*len(rows)+2.1) * 12),
          ncol=4, frameon=False, fontsize=8)
fig.suptitle('GCF-27 — the family holding the confirmed phosphonolipid (P. ananatis LMG 5342 region 2)\n'
             f'{len(rows)} members aligned on pepM · {sum(r["has_pf"] for r in rows)} carry PF01066 (top), '
             f'{sum(not r["has_pf"] for r in rows)} do not (below) · '
             f'{sum(r["edge"] for r in rows)} truncated at a contig edge',
             fontsize=10.5, y=1 - 0.30/(0.33*len(rows)+2.1))
fig.subplots_adjust(top=1 - 0.78/(0.33*len(rows)+2.1), left=0.235, right=0.985,
                    bottom=1.05/(0.33*len(rows)+2.1))
out = '/tmp/claude-1000/-home-alexp/53e768e3-330d-4a71-9869-049e46d48788/scratchpad/gcf27_alignment'
fig.savefig(out+'.png', dpi=170)
fig.savefig(out+'.svg', metadata=SVG_METADATA); canonicalise_svg(out+'.svg')
print('rows', len(rows), '| coloured groups', len(top), '| edge', sum(r['edge'] for r in rows))
print('saved', out+'.png')
