#!/usr/bin/env python3
"""Positional ablation: trim the neighbourhood, re-cluster, see if genus-confinement holds.

69 of 85 multi-member GCFs have all members in one genus. Either the chemistry is
genus-restricted, or BiG-SCAPE is clustering on the 10 kb neighbourhood, whose
housekeeping genes same-genus genomes share. The repo caught the second thing once
already: the pantaphos 186/29 split was "species-level chromosomal context, not
chemistry", proven by ablation.

A first attempt dropped CDS by domain CATEGORY (primary metabolism / mobile) and
removed 1.0% of CDS -- the curated map covers ~100 accessions, so flanking genes
mostly classify as 'other'. Useless as an ablation. Measured instead: proto_core is
11.1 kb of a 30.1 kb median region, so 63% of a region IS neighbourhood. Trim that.

Keeps proto_core plus FLANK_KB either side. The flank is not zero because
CLAUDE.md records that proto_core alone discards real tailoring chemistry on the
pantaphos cluster (LeuC/LeuD, the methyltransferase, SanS, the GNAT).

Usage: ablate_pos.py <outdir> <flank_kb>
"""
import os
import shutil
import sys
from collections import Counter

from Bio import SeqIO

SRC = '/mnt/bigdisk/clusterquest/results/antismash_results/Bacteria'
OUT, FLANK = sys.argv[1], int(float(sys.argv[2]) * 1000)
dest = f'{OUT}/trim{sys.argv[2]}kb'
os.makedirs(dest, exist_ok=True)
ctrl = f'{OUT}/control'
want_ctrl = not os.path.isdir(ctrl)
os.makedirs(ctrl, exist_ok=True)

gbks = sorted(os.path.join(r, f) for r, _, fs in os.walk(SRC)
              for f in fs if '.region' in f and f.endswith('.gbk'))
st = Counter()
lens, cdss = [], []

for i, path in enumerate(gbks, 1):
    name = os.path.basename(path)
    if want_ctrl:
        shutil.copy2(path, f'{ctrl}/{name}')
    rec = SeqIO.read(path, 'genbank')
    core = [f for f in rec.features if f.type == 'proto_core']
    if not core:
        shutil.copy2(path, f'{dest}/{name}')   # nothing to trim against
        st['no_core'] += 1
        continue
    s = max(0, min(int(f.location.start) for f in core) - FLANK)
    e = min(len(rec.seq), max(int(f.location.end) for f in core) + FLANK)

    sub = rec[s:e]                 # Biopython slices features with the sequence
    sub.id, sub.name = rec.id, rec.name
    sub.annotations = dict(rec.annotations)
    sub.annotations['molecule_type'] = 'DNA'
    # a slice drops features that straddle the cut; re-add the region/cand_cluster
    # records BiG-SCAPE's AS5 reader walks, spanning the whole trimmed record
    have = {f.type for f in sub.features}
    for ftype in ('region', 'cand_cluster', 'protocluster', 'proto_core'):
        if ftype not in have:
            src = next((f for f in rec.features if f.type == ftype), None)
            if src is not None:
                from Bio.SeqFeature import SeqFeature, SimpleLocation
                nf = SeqFeature(SimpleLocation(0, len(sub.seq), 1), type=ftype)
                nf.qualifiers = dict(src.qualifiers)
                sub.features.insert(0, nf)
    SeqIO.write(sub, f'{dest}/{name}', 'genbank')
    st['trimmed'] += 1
    lens.append((len(rec.seq), len(sub.seq)))
    cdss.append((sum(1 for f in rec.features if f.type == 'CDS'),
                 sum(1 for f in sub.features if f.type == 'CDS')))
    if i % 400 == 0:
        print(f'  {i}/{len(gbks)}')

import statistics as S
print(f"trimmed {st['trimmed']}  no proto_core {st['no_core']}  flank {FLANK/1000:g} kb")
print(f"  median bp  {S.median(a for a,_ in lens):7.0f} -> {S.median(b for _,b in lens):7.0f}"
      f"  ({100*S.median(b for _,b in lens)/S.median(a for a,_ in lens):.0f}%)")
print(f"  median CDS {S.median(a for a,_ in cdss):7.0f} -> {S.median(b for _,b in cdss):7.0f}")
print(f"  total CDS  {sum(a for a,_ in cdss):7d} -> {sum(b for _,b in cdss):7d}"
      f"  ({100*sum(b for _,b in cdss)/sum(a for a,_ in cdss):.0f}% kept)")
