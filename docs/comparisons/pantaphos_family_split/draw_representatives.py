# Recorded experiment, not a pipeline stage and not a reusable tool.
# Paths are written as $WORK (a scratch dir) and $REPO (the checkout) and must be
# substituted before running; these were one-off scripts and are committed for the
# record rather than for convenience. See summary.json.
"""Representative BGC from each of the three groups the ablation separated.

Genes come from BiG-SCAPE's own CDS + Pfam scan (the control database), not from
the region GenBank -- region GBKs carry no Pfam accessions, and the domains are
the whole point of the comparison.

Representative = the member with the highest mean Jaccard to the rest of its
group, i.e. the most typical architecture rather than an arbitrary pick.
"""
import sqlite3, collections, sys, statistics as st
from pathlib import Path
sys.path.insert(0, '$REPO/scripts')
sys.path.insert(0, '$REPO/scripts/figures')
from utils.domain_functions import category, name
from utils import plotting  # pins hashsalt + Agg
from utils.plotting import SVG_METADATA, canonicalise_svg
import matplotlib.pyplot as plt

AB = Path('$WORK/ablation')
COLOR = {'core':'#b3392b','tailoring':'#2f7d5d','lipid':'#7a5ea8','transport':'#1b5e7e',
         'regulation':'#c98a1b','mobile':'#8a939f','primary':'#cfd5dd','other':'#e4e8ed'}
INK, FAINT, ACCENT = '#15181d', '#6d7683', '#8d3a2c'
SHOW = {'PF13714':'pepM','PF13535':'ATP-grasp','PF00155':'Aminotran_1_2','PF00589':'integrase',
        'PF04754':'transposase','PF08775':'ParB','PF07690':'MFS','PF00682':'HMGL-like',
        'PF13673':'GNAT','PF13649':'MeT','PF00296':'Bac_luciferase','PF00005':'ABC_tran',
        'PF00330':'aconitase','PF01613':'Flavin_Reduct','PF00440':'TetR'}

def groups():
    def load(db):
        c=sqlite3.connect(db); c.row_factory=sqlite3.Row; o={}
        for r in c.execute("""SELECT g.path p, f.id fid FROM bgc_record br
            JOIN bgc_record_family bf ON bf.record_id=br.id JOIN family f ON f.id=bf.family_id
            JOIN gbk g ON g.id=br.gbk_id WHERE br.record_type='region' AND f.cutoff=0.3"""):
            o[r['p'].split('/')[-1]]=r['fid']
        return o
    ctl, abl = load(AB/'out_control/out_control.db'), load(AB/'out_ablated/out_ablated.db')
    bc=collections.Counter(ctl.values()).most_common(1)[0][0]
    ba=collections.Counter(abl.values()).most_common(1)[0][0]
    sc={k for k,v in ctl.items() if v!=bc}; sa={k for k,v in abl.items() if v!=ba}
    return {'186': set(ctl)-sc, '18': sc-sa, '11': sa}

c=sqlite3.connect(AB/'out_control/out_control.db'); c.row_factory=sqlite3.Row
meta={}
for r in c.execute("""SELECT g.id gid, g.path p, g.organism org, br.nt_start s, br.nt_stop e
    FROM gbk g JOIN bgc_record br ON br.gbk_id=g.id WHERE br.record_type='region'"""):
    meta[r['p'].split('/')[-1]]={'gid':r['gid'],'org':r['org'],'len':r['e']-r['s']}
genes=collections.defaultdict(list)
for r in c.execute("""SELECT cds.gbk_id gid, cds.id cid, cds.nt_start s, cds.nt_stop e, cds.strand st
                      FROM cds ORDER BY cds.nt_start"""):
    genes[r['gid']].append({'cid':r['cid'],'s':r['s'],'e':r['e'],'st':r['st'],'doms':[]})
best={}
for r in c.execute("SELECT cds_id, accession, bit_score FROM hsp ORDER BY bit_score DESC"):
    best.setdefault(r['cds_id'], r['accession'].split('.')[0])
dset=collections.defaultdict(set)
for gid,gs in genes.items():
    for g in gs:
        a=best.get(g['cid'])
        if a: g['doms']=[a]; dset[gid].add(a)

G=groups()
def rep(keys):
    ks=sorted(keys); ss={k:dset[meta[k]['gid']] for k in ks}
    def jac(a,b): return len(a&b)/len(a|b) if (a|b) else 0
    return max(ks, key=lambda k: st.mean([jac(ss[k],ss[o]) for o in ks if o!=k]) if len(ks)>1 else 1)

picks=[(lbl, rep(G[lbl])) for lbl in ('186','18','11')]
TITLE={'186':'186 members · mostly $\\it{P.\\ ananatis}$ · pantaphos WITH ATP-grasp',
       '18':'18 members · mostly $\\it{P.\\ allii}$ · pantaphos WITHOUT ATP-grasp',
       '11':'11 members · all $\\it{P.\\ agglomerans}$ · same pantaphos core, +27 kb of mobile cargo'}

span=max(meta[k]['len'] for _,k in picks)
fig,axes=plt.subplots(3,1,figsize=(11.5,6.2))
for ax,(lbl,k) in zip(axes,picks):
    m=meta[k]; gs=genes[m['gid']]
    ax.plot([0,m['len']],[0,0],color='#c9ced6',lw=1.0,zorder=1)
    for g in gs:
        w=g['e']-g['s']; hl=min(w*0.3, span*0.008); body=w-hl
        cat=category(g['doms'][0]) if g['doms'] else 'other'
        cat='primary' if cat=='primary metabolism' else cat
        col=COLOR.get(cat, '#e4e8ed')
        if g['st']>=0:
            v=[(g['s'],-0.3),(g['s']+body,-0.3),(g['e'],0),(g['s']+body,0.3),(g['s'],0.3)]
        else:
            v=[(g['e'],-0.3),(g['e']-body,-0.3),(g['s'],0),(g['e']-body,0.3),(g['e'],0.3)]
        ax.add_patch(plt.Polygon(v,closed=True,facecolor=col,edgecolor='#ffffff',lw=0.5,zorder=3))
    lab=[g for g in gs if g['doms'] and g['doms'][0] in SHOW]
    for i,g in enumerate(lab):
        mid=(g['s']+g['e'])/2; tier=i%3
        nm=SHOW[g['doms'][0]]
        hot = nm in ('ATP-grasp','Aminotran_1_2','pepM')
        ax.annotate(nm, xy=(mid,0.18), xytext=(mid,0.55+tier*0.42), ha='center',
                    va='bottom', fontsize=7.2, color=ACCENT if hot else FAINT,
                    fontweight='bold' if hot else 'normal',
                    arrowprops=dict(arrowstyle='-',lw=0.6,color=ACCENT if hot else '#c9ced6'))
    ax.set_xlim(-span*0.01, span*1.01); ax.set_ylim(-0.9, 2.1)
    ax.set_yticks([]); ax.spines[:].set_visible(False)
    ax.set_title(f"{TITLE[lbl]}   —   {m['org']}, {m['len']/1000:.1f} kb, {len(gs)} CDS",
                 fontsize=8.8, loc='left', color=INK, pad=4)
    ax.tick_params(labelsize=7.5, colors=FAINT)
    ax.set_xticks(range(0,int(span)+1,10000))
    ax.set_xticklabels([f'{x//1000}' for x in range(0,int(span)+1,10000)])
axes[-1].set_xlabel('kb within region', fontsize=8, color=FAINT)
used=['core','tailoring','transport','regulation','mobile','primary','other']
fig.legend(handles=[plt.Rectangle((0,0),1,1,facecolor=COLOR[c],edgecolor='#fff') for c in used],
           labels=['phosphonate core','tailoring','transport','regulation','mobile element',
                   'primary metabolism','unclassified'],
           loc='lower center', ncol=7, fontsize=7.4, frameon=False, bbox_to_anchor=(0.5,-0.02))
fig.tight_layout(rect=[0,0.045,1,1])
out=AB/'representatives.svg'
fig.savefig(out, format='svg', metadata=SVG_METADATA, bbox_inches='tight')
canonicalise_svg(out)
fig.savefig(str(out).replace('.svg','.png'), dpi=150, bbox_inches='tight')
print('wrote', out)
for lbl,k in picks:
    print(f"  {lbl:>4}: {k}  ({meta[k]['org']})")
