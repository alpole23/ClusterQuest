# Recorded experiment, not a pipeline stage and not a reusable tool.
# Paths are written as $WORK (a scratch dir) and $REPO (the checkout) and must be
# substituted before running; these were one-off scripts and are committed for the
# record rather than for convenience. See summary.json.
"""Rebuild every region as antiSMASH would have drawn it without CUTOFF-chaining.

Removing just the offending Aminotran_1_2 CDS would be the wrong experiment: the
gene is not the problem, the BOUNDARY it dragged out is. The cargo that entered
through the inflated window -- integrase, transposase, ABC transporters -- would
stay, and the regions would still look different.

So: recompute the core as it would be WITHOUT the gap-crossing edge gene, put the
10 kb neighbourhood back on each side, and drop every CDS outside that window.
The GBK keeps its region feature and sequence so BiG-SCAPE can still load it;
only gene content changes, and gene content is all its distance depends on
(jaccard = shared domains, adjacency = domain order, DSS = domain similarity).

All 334 regions go through the identical code path, so an unaffected region comes
out unchanged and there is no processing asymmetry between arms.
"""
import sqlite3, collections, shutil, sys
from pathlib import Path
from Bio import SeqIO

OUT=Path('$WORK/gcf5')
ctrl, trt = OUT/'control', OUT/'trimmed'
for d in (ctrl,trt): shutil.rmtree(d,ignore_errors=True); d.mkdir(parents=True)
NB=10000; GAP=5000
PRIMARY={'PF13714'}
RULE={'PF00155','PF00171','PF00266','PF00291','PF00465','PF00682','PF01467','PF02775',
      'PF02776','PF02826','PF11583','PF12804','PF12850','PF14907','PF16383','PF13714'}
c=sqlite3.connect('results_erw_verify/bigscape_results/Erwiniaceae/Erwiniaceae.db'); c.row_factory=sqlite3.Row
rule_cds=collections.defaultdict(list)
for r in c.execute("""SELECT cds.gbk_id gid, cds.nt_start s, cds.nt_stop e, hsp.accession a
                      FROM cds JOIN hsp ON hsp.cds_id=cds.id"""):
    a=r['a'].split('.')[0]
    if a in RULE: rule_cds[r['gid']].append((r['s'],r['e'],a))

changed=0; dropped_tot=0
for r in c.execute("""SELECT br.gbk_id gid, g.path p FROM bgc_record br JOIN gbk g ON g.id=br.gbk_id
                      JOIN bgc_record_family bf ON bf.record_id=br.id JOIN family f ON f.id=bf.family_id
                      WHERE br.record_type='region' AND f.cutoff=0.3"""):
    p=Path(r['p']); stem=f'{p.parent.name}__{p.name}'
    shutil.copy(p, ctrl/stem)
    cores=list(c.execute("SELECT nt_start s,nt_stop e FROM bgc_record WHERE gbk_id=? AND record_type='proto_core'",(r['gid'],)))
    recs=list(SeqIO.parse(str(p),'genbank'))
    if cores:
        cs,ce=min(x['s'] for x in cores), max(x['e'] for x in cores)
        by=collections.defaultdict(set)
        for s,e,a in rule_cds[r['gid']]:
            if s>=cs-1 and e<=ce+1: by[(s,e)].add(a)
        pos=sorted(by)
        while len(pos)>1 and not (by[pos[0]] & PRIMARY) and pos[1][0]-pos[0][1]>=GAP: pos.pop(0)
        while len(pos)>1 and not (by[pos[-1]] & PRIMARY) and pos[-1][0]-pos[-2][1]>=GAP: pos.pop()
        if pos:
            lo,hi=max(0,pos[0][0]-NB), pos[-1][1]+NB
            if lo>cs or hi<ce:
                changed+=1
                for rec in recs:
                    keep=[]
                    for f in rec.features:
                        if f.type=='CDS' and not (int(f.location.start)>=lo and int(f.location.end)<=hi):
                            dropped_tot+=1; continue
                        keep.append(f)
                    rec.features=keep
    SeqIO.write(recs, str(trt/stem),'genbank')
print(f'regions written per arm: {len(list(ctrl.iterdir()))}')
print(f'regions whose window shrank: {changed};  CDS dropped: {dropped_tot}')
