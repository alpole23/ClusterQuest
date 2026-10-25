#!/usr/bin/env python3
"""HGT candidate display: arcs over a GTDB-ordered axis.

Why an arc diagram and not a scatter. Measured on this run, median BGC distance
is 0.93-0.97 at EVERY level of taxonomic separation -- same species through
different phylum. There is no trend to plot, which is the same reason this repo
removed the Yu et al. similarity-vs-taxonomy figure. The whole signal is in the
tail: pairs inside the 0.30 GCF cutoff collapse 846 (same genus) -> 149 (family)
-> 1 (order) -> 3 (class) -> 1 (different phylum).

So draw the tail, over a backbone that makes distance legible. Genomes sit on a
GTDB-ordered axis; every BGC pair within the cutoff gets an arc. Sharing is
overwhelmingly local, so short arcs pile into a band at the baseline and the few
long-range ones leap clear of it -- the exceptions are visible against the rule
rather than cherry-picked out of it.
"""
import json
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict

sys.path.insert(0, '/home/alexp/pipeline/scripts')
from utils import plotting  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Rectangle

D = '/mnt/bigdisk/clusterquest/work_pilot/4c/121b1cc8629a33e31769c842ab7b11'
OUT = '/tmp/claude-1000/-home-alexp/53e768e3-330d-4a71-9869-049e46d48788/scratchpad'
CUTOFF = 0.30
RANKS = [('species', 's'), ('genus', 'g'), ('family', 'f'),
         ('order', 'o'), ('class', 'c'), ('phylum', 'p')]
LV = [r[0] for r in RANKS]

tx = defaultdict(dict)
for line in open(f'{D}/gtdbtk.bac120.summary.tsv').readlines()[1:]:
    p = line.split('\t'); g = p[0].replace('usr_', '', 1)
    for n, c in RANKS:
        m = re.search(rf'{c}__([A-Za-z0-9_.\- ]+)', p[1])
        tx[g][n] = m.group(1).strip() if m and m.group(1).strip() else None

gen = {}
for rec in json.load(open(f'{D}/phosphonate_metadata.json')):
    pr = os.path.normpath(rec['gbk_path']).split(os.sep)
    if 'antismash_input' in pr:
        gen[rec['label']] = pr[pr.index('antismash_input') + 1]

con = sqlite3.connect(f'file:{D}/Bacteria.db?mode=ro', uri=True)
rid = {i: os.path.basename(p)[:-4] for i, p in con.execute(
    "select r.id,g.path from bgc_record r join gbk g on g.id=r.gbk_id "
    "where r.record_type='region'")}
fam = {os.path.basename(p)[:-4]: f for p, f in con.execute(
    "select g.path,f.family_id from bgc_record r join gbk g on g.id=r.gbk_id "
    "join bgc_record_family f on f.record_id=r.id join family y on y.id=f.family_id "
    "where r.record_type='region' and y.cutoff=0.3")}


def sep(a, b):
    for k, n in enumerate(LV):
        if tx[a][n] is not None and tx[a][n] == tx[b][n]:
            return k
    return 6


# axis: one slot per GENOME, ordered by full GTDB lineage so related sit adjacent
genomes = sorted({g for g in gen.values() if g in tx},
                 key=lambda g: tuple(tx[g][n] or 'zz' for n in reversed(LV)))
xi = {g: i for i, g in enumerate(genomes)}
N = len(genomes)

arcs = []
for a, b, d in con.execute('select record_a_id,record_b_id,distance from distance '
                           'where distance < ?', (CUTOFF,)):
    la, lb = rid.get(a), rid.get(b)
    ga, gb = gen.get(la), gen.get(lb)
    if ga in xi and gb in xi and ga != gb:
        arcs.append((xi[ga], xi[gb], sep(ga, gb), d, fam.get(la), ga, gb))

# sequential ramp for within-family separation; the accent is reserved for the
# beyond-family pairs, which are the figure's subject and carry direct labels
RAMP = {0: '#cfe0e8', 1: '#9ec4d4', 2: '#5f9bb4'}
ACCENT = '#8d3a2c'
LABELS = {0: 'same species', 1: 'same genus', 2: 'same family',
          3: 'same order', 4: 'same class', 5: 'same phylum', 6: 'different phylum'}

fig = plt.figure(figsize=(15, 6.4))
ax = fig.add_axes([0.045, 0.20, 0.80, 0.70])
far = [a for a in arcs if a[2] >= 3]
near = [a for a in arcs if a[2] < 3]

for x0, x1, k, d, f, ga, gb in near:
    w = abs(x1 - x0)
    ax.add_patch(Arc(((x0 + x1) / 2, 0), w, 2 * w ** 0.78, theta1=0, theta2=180,
                     lw=0.55, color=RAMP[k], alpha=0.5, zorder=1))
for x0, x1, k, d, f, ga, gb in far:
    w = abs(x1 - x0)
    ax.add_patch(Arc(((x0 + x1) / 2, 0), w, 2 * w ** 0.78, theta1=0, theta2=180,
                     lw=2.0, color=ACCENT, alpha=0.95, zorder=3))

ymax = max((abs(a[1] - a[0]) ** 0.78 for a in arcs), default=1) * 1.18
ax.set_xlim(-N * 0.012, N * 1.012)
ax.set_ylim(-ymax * 0.055, ymax)
ax.axis('off')

# phylum bands along the baseline
start, cur = 0, tx[genomes[0]]['phylum']
bands = []
for i, g in enumerate(genomes + [None]):
    ph = tx[g]['phylum'] if g else None
    if ph != cur:
        bands.append((start, i, cur)); start, cur = i, ph
pal = ['#1b5e7e', '#2f7d5d', '#6d7683', '#9aa5b1', '#c2a35b', '#7a5ea8']
for j, (s0, s1, ph) in enumerate(bands):
    ax.add_patch(Rectangle((s0, -ymax * 0.05), s1 - s0, ymax * 0.035,
                           facecolor=pal[j % len(pal)], edgecolor='none', zorder=4))
    if s1 - s0 > N * 0.035:
        ax.text((s0 + s1) / 2, -ymax * 0.062, (ph or '?')[:18], ha='center',
                va='top', fontsize=6.6, color='#15181d')

# direct labels on the beyond-family arcs
seen = set()
for x0, x1, k, d, f, ga, gb in sorted(far, key=lambda a: -abs(a[1] - a[0])):
    key = (f, k)
    if key in seen:
        continue
    seen.add(key)
    npair = sum(1 for z in far if (z[4], z[2]) == key)
    w = abs(x1 - x0); top = w ** 0.78
    ax.annotate(f'GCF-{f}  d={d:.2f}' + (f'  ({npair} pairs)' if npair > 1 else '') + f'\n{LABELS[k]}\n'
                f'{tx[ga]["genus"]} ↔ {tx[gb]["genus"]}',
                xy=((x0 + x1) / 2, top), xytext=((x0 + x1) / 2, top + ymax * 0.085),
                ha='center', fontsize=7.0, color=ACCENT, fontweight='bold',
                arrowprops=dict(arrowstyle='-', color=ACCENT, lw=0.8), zorder=5)

h = [plt.Line2D([], [], color=RAMP[k], lw=2, label=f'{LABELS[k]}  (n={sum(1 for a in arcs if a[2]==k)})')
     for k in (0, 1, 2)]
h.append(plt.Line2D([], [], color=ACCENT, lw=2.4,
                    label=f'beyond family  (n={len(far)})'))
ax.legend(handles=h, loc='upper left', bbox_to_anchor=(1.005, 1.0), frameon=False,
          fontsize=8, title='BGC pair within 0.30 cutoff', title_fontsize=8.5)

fig.text(0.045, 0.965, 'Phosphonate BGCs shared between genomes, over a GTDB-ordered axis',
         fontsize=13, fontweight='bold', color='#15181d')
fig.text(0.045, 0.925,
         f'Each arc joins two genomes whose BGCs fall inside the {CUTOFF} GCF cutoff. '
         f'{len(genomes):,} BGC-positive genomes ordered by GTDB lineage, so arc width is '
         f'taxonomic distance. Median BGC distance is 0.93–0.97 at every level of '
         f'separation, so the signal is the tail drawn here, not a trend.',
         fontsize=8, color='#6d7683')

out = f'{OUT}/hgt_arcs.png'
fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
print(f'{len(arcs)} arcs ({len(far)} beyond family) over {N} genomes -> {out}')
