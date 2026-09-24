# Recorded experiment, not a pipeline stage and not a reusable tool.
# Paths are written as $WORK (a scratch dir) and $REPO (the checkout) and must be
# substituted before running; these were one-off scripts and are committed for the
# record rather than for convenience. See summary.json.
"""Core extension by CUTOFF-chaining, measured by GAP rather than by profile identity.

First attempt trimmed every secondary-profile gene from the core edges and found
a median core of 3.3 -> 0.9 kb. That over-trimmed: HMGL-like is a secondary
profile AND phosphonomethylmalate synthase, a genuine pantaphos gene, so the
measure was deleting real cluster content and calling it inflation.

The signature of CUTOFF-chaining is a LARGE GAP: an edge rule-gene separated
from the rest of the core by kilobases of non-rule sequence. That is measurable
without deciding which profiles deserve to be in the rule.
"""
import sqlite3, collections, statistics as st
DB='results_erw_verify/bigscape_results/Erwiniaceae/Erwiniaceae.db'
PRIMARY={'PF13714'}
NAMES={'PF00155':'Aminotran_1_2','PF00171':'Aldedh','PF00266':'Aminotran_5','PF00291':'PALP',
 'PF00465':'Fe-ADH','PF00682':'HMGL-like','PF01467':'CTP_transf_like','PF02775':'TPP_enzyme_C',
 'PF02776':'TPP_enzyme_N','PF02826':'2-Hacid_dh_C','PF11583':'AurF','PF12804':'NTP_transf_3',
 'PF12850':'Metallophos_2','PF14907':'NTP_transf_5','PF16383':'DUF4992','PF13714':'pepM'}
RULE=set(NAMES)
c=sqlite3.connect(DB); c.row_factory=sqlite3.Row
cds=collections.defaultdict(list)
for r in c.execute("""SELECT cds.gbk_id gid, cds.nt_start s, cds.nt_stop e, hsp.accession a
                      FROM cds JOIN hsp ON hsp.cds_id=cds.id"""):
    a=r['a'].split('.')[0]
    if a in RULE: cds[r['gid']].append((r['s'],r['e'],a))

GAP=5000
rows=[]; culprits=collections.Counter(); fam=collections.defaultdict(list)
for r in c.execute("""SELECT br.gbk_id gid, f.id fid, g.organism org FROM bgc_record br
    JOIN bgc_record_family bf ON bf.record_id=br.id JOIN family f ON f.id=bf.family_id
    JOIN gbk g ON g.id=br.gbk_id WHERE br.record_type='region' AND f.cutoff=0.3"""):
    cores=list(c.execute("SELECT nt_start s,nt_stop e FROM bgc_record WHERE gbk_id=? AND record_type='proto_core'",(r['gid'],)))
    if not cores: continue
    cs,ce=min(x['s'] for x in cores), max(x['e'] for x in cores)
    by=collections.defaultdict(set)
    for s,e,a in cds[r['gid']]:
        if s>=cs-1 and e<=ce+1: by[(s,e)].add(a)
    pos=sorted(by)
    if not pos: continue
    shed=0; culp=[]
    # peel an edge gene only when a >=GAP gap separates it from the next rule gene
    while len(pos)>1 and not (by[pos[0]] & PRIMARY) and pos[1][0]-pos[0][1] >= GAP:
        shed += pos[1][0]-pos[0][0]; culp += sorted(by[pos[0]]); pos.pop(0)
    while len(pos)>1 and not (by[pos[-1]] & PRIMARY) and pos[-1][0]-pos[-2][1] >= GAP:
        shed += pos[-1][1]-pos[-2][1]; culp += sorted(by[pos[-1]]); pos.pop()
    rows.append(shed)
    if shed>0:
        for a in set(culp): culprits[NAMES[a]]+=1
        fam[r['fid']].append((shed, r['org'] or '?'))

n=len(rows); hit=[x for x in rows if x>0]
print(f"regions examined: {n}")
print(f"core extended by CUTOFF-chaining across a >={GAP//1000} kb gap: {len(hit)}  ({len(hit)/n:.1%})")
if hit:
    print(f"  median extension {st.median(hit)/1000:.1f} kb   max {max(hit)/1000:.1f} kb")
print("\nprofile sitting on the far side of the gap:")
for nm,k in culprits.most_common():
    print(f"  {nm:16s} {k:3d} regions")
print("\nby family:")
for fid,v in sorted(fam.items(), key=lambda kv:-len(kv[1])):
    tot=sum(1 for r in c.execute("""SELECT br.id FROM bgc_record br JOIN bgc_record_family bf
        ON bf.record_id=br.id JOIN family f ON f.id=bf.family_id
        WHERE br.record_type='region' AND f.cutoff=0.3 AND f.id=?""",(fid,)))
    print(f"  GCF {fid:>2}  {len(v):3d}/{tot:<3d}  median {st.median([x[0] for x in v])/1000:5.1f} kb   {v[0][1][:36]}")
