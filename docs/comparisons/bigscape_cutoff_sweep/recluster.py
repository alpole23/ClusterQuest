"""Re-run BiG-SCAPE's own family assignment at several cutoffs.

Uses big_scape.network.families.generate_families unchanged, so this is the real
affinity-propagation step with the real parameters -- not a reimplementation. Only
the surrounding plumbing (connected components, persistence) is local, because the
full CLI path insists on a Pfam HMM it does not need once every CDS is scanned.
"""
import collections, sqlite3, sys, json

import big_scape.data as bs_data
from big_scape.cli.config import BigscapeConfig
from big_scape.network.families import generate_families
from pathlib import Path

DB_PATH = Path('work.db')
CUTOFFS = [float(c) for c in sys.argv[1].split(',')]
EDGE_PARAM = 1

BigscapeConfig.parse(Path(BigscapeConfig.__module__ and
    __import__('big_scape').__path__[0]) / 'config.yml', None) \
    if hasattr(BigscapeConfig, 'parse') else None
print(f'AP config: PREFERENCE={BigscapeConfig.PREFERENCE} '
      f'DENSITY={BigscapeConfig.DENSITY} DENSE_PREFERENCE={BigscapeConfig.DENSE_PREFERENCE}')

bs_data.DB.load_from_disk(DB_PATH)

raw = sqlite3.connect(f'file:{DB_PATH}?mode=ro', uri=True)
regions = {r for (r,) in raw.execute(
    "SELECT id FROM bgc_record WHERE record_type='region'")}
edges_all = raw.execute(
    "SELECT record_a_id, record_b_id, distance, jaccard, adjacency, dss, edge_param_id "
    "FROM distance WHERE edge_param_id=? ORDER BY rowid", (EDGE_PARAM,)).fetchall()
print(f'{len(regions):,} region records, {len(edges_all):,} stored comparisons')

def components(edges):
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for a, b, *_ in edges:
        union(a, b)
    out = collections.defaultdict(list)
    for n in {n for e in edges for n in e[:2]}:
        out[find(n)].append(n)
    return list(out.values())

result = {}
for cut in CUTOFFS:
    sub = [e for e in edges_all if e[2] < cut]
    comps = components(sub)
    by_node = {}
    for i, c in enumerate(comps):
        for n in c: by_node[n] = i
    buckets = collections.defaultdict(list)
    for e in sub:
        buckets[by_node[e[0]]].append(e)

    fam_of = {}
    for i, cc_edges in buckets.items():
        for rec, centre, _c, _b, _r in generate_families(cc_edges, 'mix', cut, 1):
            fam_of[rec] = centre
    for r in regions:                        # singletons: their own family
        fam_of.setdefault(r, r)
    fam_of = {r: f for r, f in fam_of.items() if r in regions}

    sizes = collections.Counter(fam_of.values())
    vals = sorted(sizes.values())
    result[cut] = fam_of
    print(f'\ncutoff {cut:g}: {len(comps)} connected components -> '
          f'{len(sizes)} families')
    print(f'   largest {vals[-1]}   singletons {sum(1 for v in vals if v == 1)}   '
          f'median {vals[len(vals)//2]}   regions in a family >1 '
          f'{sum(v for v in vals if v > 1)}')

json.dump({str(k): {str(a): b for a, b in v.items()} for k, v in result.items()},
          open('memberships.json', 'w'))
print('\nwritten memberships.json')
