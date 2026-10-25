#!/usr/bin/env python3
"""Did trimming the neighbourhood merge the genus-confined GCFs?

Three clusterings get compared:
  pipeline  what the run produced (full regions)
  control   the same full regions, re-clustered -- the NOISE FLOOR, because
            BiG-SCAPE does not reproduce itself exactly
  trim2kb   proto_core +/- 2 kb, 51% of CDS removed

The claim under test is that 69 of 85 multi-member GCFs are genus-confined
because BiG-SCAPE sees shared flanking housekeeping, not shared chemistry. If
so, trimming should merge families ACROSS genera and push the "shallowest rank
uniting all members" distribution upward. Control-vs-pipeline says how much
movement happens for no reason at all; only movement beyond that counts.
"""
import glob
import json
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from itertools import combinations

S = '/tmp/claude-1000/-home-alexp/53e768e3-330d-4a71-9869-049e46d48788/scratchpad'
D = '/mnt/bigdisk/clusterquest/work_pilot/4c/121b1cc8629a33e31769c842ab7b11'
RANKS = [('species', 's'), ('genus', 'g'), ('family', 'f'),
         ('order', 'o'), ('class', 'c'), ('phylum', 'p')]


def taxonomy():
    tx = defaultdict(dict)
    for line in open(f'{D}/gtdbtk.bac120.summary.tsv').readlines()[1:]:
        p = line.split('\t')
        g = p[0].replace('usr_', '', 1)
        for name, code in RANKS:
            m = re.search(rf'{code}__([A-Za-z0-9_.\- ]+)', p[1])
            tx[g][name] = m.group(1).strip() if m and m.group(1).strip() else None
    return tx


def membership(db):
    """region basename -> family id, at cutoff 0.30."""
    con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
    out = {}
    for path, fid in con.execute(
            "select g.path,f.family_id from bgc_record r join gbk g on g.id=r.gbk_id "
            "join bgc_record_family f on f.record_id=r.id join family y on y.id=f.family_id "
            "where r.record_type='region' and y.cutoff=0.3"):
        out[os.path.basename(path)[:-4]] = fid
    return out


def genome_of():
    g = {}
    for rec in json.load(open(f'{D}/phosphonate_metadata.json')):
        pr = os.path.normpath(rec['gbk_path']).split(os.sep)
        if 'antismash_input' in pr:
            g[rec['label']] = pr[pr.index('antismash_input') + 1]
    return g


def depth_profile(mem, gen, tx):
    """How deep you must go to unite every member of each multi-member family."""
    members = defaultdict(list)
    for lab, fid in mem.items():
        gg = gen.get(lab)
        if gg in tx:
            members[fid].append(gg)
    prof, crossers = Counter(), []
    for fid, gs in members.items():
        if len(gs) < 2:
            continue
        lvl = 'above phylum'
        for name, _ in RANKS:
            vals = {tx[g][name] for g in gs}
            if len(vals) == 1 and None not in vals:
                lvl = name
                break
        prof[lvl] += 1
        if lvl not in ('species', 'genus'):
            crossers.append((fid, len(gs), lvl, len({tx[g]['genus'] for g in gs})))
    return prof, crossers, members


def comembership(mem, shared):
    by = defaultdict(set)
    for lab in shared:
        by[mem[lab]].add(lab)
    return {frozenset(p) for s in by.values() for p in combinations(sorted(s), 2)}


tx, gen = taxonomy(), genome_of()
runs = {'pipeline': f'{D}/Bacteria.db'}
for arm in ('control', 'trim2kb'):
    hits = glob.glob(f'{S}/ablation/out_{arm}/**/*.db', recursive=True)
    if hits:
        runs[arm] = max(hits, key=os.path.getsize)

mems, profs = {}, {}
print(f'{"run":10s} {"regions":>8s} {"families":>9s} {"multi":>6s} {"singletons":>11s}')
for name, db in runs.items():
    m = membership(db)
    mems[name] = m
    sz = Counter(m.values())
    p, c, _ = depth_profile(m, gen, tx)
    profs[name] = (p, c)
    print(f'{name:10s} {len(m):8d} {len(sz):9d} '
          f'{sum(1 for v in sz.values() if v>1):6d} {sum(1 for v in sz.values() if v==1):11d}')

print('\nShallowest rank uniting all members of a multi-member family:')
order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'above phylum']
hdr = ''.join(f'{n:>14s}' for n in runs)
print(f'{"":14s}{hdr}')
for lvl in order:
    row = ''.join(f'{profs[n][0][lvl]:14d}' for n in runs)
    if any(profs[n][0][lvl] for n in runs):
        print(f'{lvl:14s}{row}')
for n in runs:
    above = sum(v for k, v in profs[n][0].items() if k not in ('species', 'genus'))
    tot = sum(profs[n][0].values())
    print(f'  {n}: {above} of {tot} families cross a genus boundary '
          f'({100*above/max(tot,1):.1f}%)')

print('\nCo-membership change (pairs of regions in the same family):')
base = set(mems['pipeline'])
for a, b in (('pipeline', 'control'), ('control', 'trim2kb')):
    if b not in mems:
        continue
    shared = sorted(base & set(mems[a]) & set(mems[b]))
    pa, pb = comembership(mems[a], shared), comembership(mems[b], shared)
    print(f'  {a:9s} -> {b:9s}  shared {len(shared)} regions  '
          f'split {len(pa-pb):6d}  merged {len(pb-pa):6d}  '
          f'jaccard {len(pa&pb)/max(len(pa|pb),1):.4f}')
