# Recorded experiment, not a pipeline stage and not a reusable tool.
# Paths are written as $WORK (a scratch dir) and $REPO (the checkout) and must be
# substituted before running; these were one-off scripts and are committed for the
# record rather than for convenience. See summary.json.
import sqlite3, collections, itertools
def fams(db):
    c=sqlite3.connect(db); c.row_factory=sqlite3.Row
    out={}
    for r in c.execute("""SELECT g.path p, f.id fid FROM bgc_record br
        JOIN bgc_record_family bf ON bf.record_id=br.id JOIN family f ON f.id=bf.family_id
        JOIN gbk g ON g.id=br.gbk_id WHERE br.record_type='region' AND f.cutoff=0.3"""):
        out[r['p'].split('/')[-1]] = r['fid']
    return out
ctl, abl = fams('out_control/out_control.db'), fams('out_ablated/out_ablated.db')
print(f"control: {len(ctl)} regions, {len(set(ctl.values()))} families "
      f"{sorted(collections.Counter(ctl.values()).values(), reverse=True)}")
print(f"ablated: {len(abl)} regions, {len(set(abl.values()))} families "
      f"{sorted(collections.Counter(abl.values()).values(), reverse=True)}")
shared=sorted(set(ctl)&set(abl))
print(f"shared regions: {len(shared)}")
def pairs(m,k):
    return {(a,b) for a,b in itertools.combinations(k,2) if m[a]==m[b]}
pc, pa = pairs(ctl,shared), pairs(abl,shared)
print(f"\nco-membership pairs  control {len(pc)}  ablated {len(pa)}  agreed {len(pc&pa)}")
print(f"split by ablation (together in control, apart after): {len(pc-pa)}")
print(f"merged by ablation (apart in control, together after): {len(pa-pc)}")
