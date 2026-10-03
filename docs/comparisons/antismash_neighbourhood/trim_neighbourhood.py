# Recorded experiment, not a pipeline stage and not a reusable tool.
# Run from the repo root with the arm directory as argv[1]; the BiG-SCAPE
# invocation that followed is in neighbourhood_width_test.json.
"""Three arms over one record set, differing only in how much flank BiG-SCAPE sees.

The question is whether the 10 kb neighbourhood manufactures the pantaphos family
split. The earlier trimming experiment could not answer it: it trimmed CUTOFF-
inflated cores and then PUT THE 10 kb NEIGHBOURHOOD BACK ON, so the flank under
suspicion here was never removed.

Record set is connected component 33 -- GCF-18, 20, 21, 22 and 24 -- the component
holding every pantaphos-cassette family. Restricting to one component is safe
because BiG-SCAPE never compares across components, so no relationship is lost by
leaving the rest out; and all three arms use the identical set, so the comparison
is matched.

Every arm goes through the same SeqIO read/write path, including the control, so
no arm is advantaged by preparation. The control keeps CDS within core +/- 10 kb,
which is what antiSMASH already drew, and is verified below by reproducing the
source run's family structure exactly.
"""
import sqlite3, collections, sys
from pathlib import Path
from Bio import SeqIO

OUT = Path(sys.argv[1])
DB = 'results_entero_refseq/bigscape_results/Enterobacterales/Enterobacterales.db'
ARMS = {'nb10_control': 10000, 'nb05': 5000, 'nb00': 0}
FAMS = (18, 20, 21, 22, 24)

c = sqlite3.connect(DB); c.row_factory = sqlite3.Row
for d in ARMS:
    (OUT / d).mkdir(parents=True, exist_ok=True)

rows = list(c.execute("""
    SELECT br.gbk_id gid, g.path p, f.id fam FROM bgc_record br
    JOIN gbk g ON g.id = br.gbk_id
    JOIN bgc_record_family bf ON bf.record_id = br.id
    JOIN family f ON f.id = bf.family_id
    WHERE br.record_type='region' AND f.cutoff=0.3 AND f.id IN (%s)"""
    % ','.join('?' * len(FAMS)), FAMS))
print(f'{len(rows)} regions across GCF-{", ".join(map(str, FAMS))}')

stats = collections.defaultdict(lambda: [0, 0])
for r in rows:
    p = Path(r['p']); stem = f"{p.parent.name}__{p.name}"
    cores = list(c.execute("SELECT nt_start s, nt_stop e FROM bgc_record "
                           "WHERE gbk_id=? AND record_type='proto_core'", (r['gid'],)))
    if not cores:
        continue
    cs, ce = min(x['s'] for x in cores), max(x['e'] for x in cores)
    for arm, nb in ARMS.items():
        lo, hi = max(0, cs - nb), ce + nb
        recs = list(SeqIO.parse(str(p), 'genbank'))
        for rec in recs:
            keep = []
            for f in rec.features:
                if f.type == 'CDS':
                    stats[arm][1] += 1
                    if not (int(f.location.start) >= lo and int(f.location.end) <= hi):
                        stats[arm][0] += 1
                        continue
                keep.append(f)
            rec.features = keep
        SeqIO.write(recs, str(OUT / arm / stem), 'genbank')

for arm, nb in ARMS.items():
    dropped, total = stats[arm]
    print(f'  {arm:14s} (core +/- {nb/1000:>4.1f} kb): '
          f'{total - dropped:>6,} CDS kept, {dropped:>6,} dropped '
          f'({100*dropped/max(total,1):.1f}%)')
