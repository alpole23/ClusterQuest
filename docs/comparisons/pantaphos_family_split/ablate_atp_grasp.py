# Recorded experiment, not a pipeline stage and not a reusable tool.
# Paths are written as $WORK (a scratch dir) and $REPO (the checkout) and must be
# substituted before running; these were one-off scripts and are committed for the
# record rather than for convenience. See summary.json.
"""Control (215 regions as-is) vs ablated (ATP-grasp CDS removed).

Targets are taken from BiG-SCAPE's OWN cds/hsp tables by coordinate, not from
antiSMASH's PFAM_domain features -- the two scans disagree, and the first
attempt at this stripped 40 features where BiG-SCAPE sees 206 hits. What must
be removed is what BiG-SCAPE scanned.

The control arm exists because subsetting 334 regions to 215 changes the
clustering landscape on its own; without it the ablation is confounded.
"""
import sqlite3, shutil, sys
from pathlib import Path
from Bio import SeqIO

OUT = Path('$WORK/ablation')
ctrl, abl = OUT/'control', OUT/'ablated'
for d in (ctrl, abl):
    shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)

c = sqlite3.connect('results_erw_verify/bigscape_results/Erwiniaceae/Erwiniaceae.db')
members = list(c.execute("""SELECT g.id gid, g.path FROM bgc_record br
 JOIN bgc_record_family bf ON bf.record_id=br.id JOIN family f ON f.id=bf.family_id
 JOIN gbk g ON g.id=br.gbk_id
 WHERE br.record_type='region' AND f.cutoff=0.3 AND f.id IN (1,2)"""))

targets = {}
for gid, _ in members:
    targets[gid] = [(r[0], r[1]) for r in c.execute(
        """SELECT cds.nt_start, cds.nt_stop FROM cds JOIN hsp ON hsp.cds_id=cds.id
           WHERE cds.gbk_id=? AND hsp.accession LIKE 'PF13535%'""", (gid,))]

missing = [p for _, p in members if not Path(p).exists()]
if missing:
    sys.exit(f'{len(missing)} source GBKs missing, e.g. {missing[0]}')

want = sum(len(v) for v in targets.values())
dropped_cds = dropped_any = 0
regions_hit = 0
for gid, p in members:
    p = Path(p)
    stem = f'{p.parent.name}__{p.name}'
    shutil.copy(p, ctrl/stem)
    spans = targets[gid]
    if spans:
        regions_hit += 1
    recs = list(SeqIO.parse(str(p), 'genbank'))
    for rec in recs:
        keep = []
        for f in rec.features:
            s, e = int(f.location.start), int(f.location.end)
            # BiG-SCAPE nt_start is 0-based like Biopython's location.start
            if any(s == ts and e == te for ts, te in spans):
                dropped_any += 1
                if f.type == 'CDS':
                    dropped_cds += 1
                continue
            keep.append(f)
        rec.features = keep
    SeqIO.write(recs, str(abl/stem), 'genbank')

print(f'{len(members)} regions per arm; {regions_hit} carry an ATP-grasp')
print(f'BiG-SCAPE PF13535 CDS in these regions : {want}')
print(f'CDS features dropped                   : {dropped_cds}')
print(f'all features dropped (CDS + overlapping): {dropped_any}')
if dropped_cds != want:
    print('WARNING: coordinate match incomplete -- do not trust the ablation')
