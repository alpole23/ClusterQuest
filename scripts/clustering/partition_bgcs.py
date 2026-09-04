#!/usr/bin/env python3
"""Assign BGCs to BiG-SCAPE partitions by pepM sequence identity.

BiG-SCAPE compares every BGC against every other and its memory goes quadratic
above ~4,000 BGCs — `GB = 1.14 + 1.29e-7*n^2`, so ~1.9 TB at the 121,000 BGCs a
million genomes would yield. Splitting the input first by pepM identity brings
the largest job to ~84 GB while rebuilding the identical GCF network
(ARI 1.0000 on three independent sets; see CLAUDE.md).

This runs *before* BiG-SCAPE, so unlike `analysis/pepm_all_by_all.py` it cannot
read pepM sequences out of the clustering database — there is none yet. It
extracts CDS translations from the region GenBanks and finds pepM with
`hmmsearch` against PF13714, which is what antiSMASH's phosphonate rule keys on
and is therefore present in every region by construction.

**Why single linkage.** Two BGCs must land in the same partition whenever they
could possibly cluster together. Single linkage is the permissive choice: it
merges on *any* qualifying link, so it never separates a pair that a stricter
criterion would have joined. Average or complete linkage could split a family
that BiG-SCAPE would have kept.

**Choosing the threshold.** 0.60-0.80 is the verified safe window. Every
same-GCF pair measured had pepM identity >= 0.901, so 0.60 leaves wide margin;
0.90 splits real families. The cut is deliberately conservative because the
failure is silent — a split family looks like two families, not like an error.

    python scripts/clustering/partition_bgcs.py \\
        --antismash_dir antismash_input --pfam Pfam-A.hmm \\
        --threshold 0.60 --max_partition_size 20000 --out partitions.tsv
"""
import argparse
import collections
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

PEPM_ACCESSION = 'PF13714'
GAP = 0


def region_gbks(root):
    """Every antiSMASH region GenBank under root.

    Matches BiG-SCAPE's own filter, which keeps .gbk files whose name contains
    "region" or "cluster" — so the whole-genome summary GenBank is excluded here
    for the same reason it is excluded there.
    """
    return sorted(p for p in Path(root).rglob('*.gbk')
                  if 'region' in p.name.lower() or 'cluster' in p.name.lower())


def extract_proteins(gbks, dest):
    """Write every CDS translation, tagged by the region it came from."""
    from Bio import SeqIO
    n = 0
    with dest.open('w') as fh:
        for i, g in enumerate(gbks):
            try:
                for rec in SeqIO.parse(g, 'genbank'):
                    for j, feat in enumerate(rec.features):
                        if feat.type != 'CDS':
                            continue
                        aa = feat.qualifiers.get('translation', [None])[0]
                        if aa:
                            fh.write(f'>{i}|{j}\n{aa}\n')
                            n += 1
            except Exception as exc:                      # a corrupt region must
                print(f'  skipped {g.name}: {exc}', file=sys.stderr)  # not sink the run
    return n


def resolve_accession(pfam_hmm, accession):
    """Pfam stores versioned accessions; hmmfetch keys on the exact string."""
    with open(pfam_hmm) as fh:
        for line in fh:
            if line.startswith('ACC ') and line.split()[1].split('.')[0] == accession:
                return line.split()[1]
    raise SystemExit(f'{accession} not found in {pfam_hmm}')


def fetch_profile(pfam_hmm, accession, workdir, hmmfetch):
    workdir.mkdir(parents=True, exist_ok=True)
    link = workdir / 'Pfam-A.hmm'
    if not link.exists():
        os.symlink(Path(pfam_hmm).resolve(), link)
    if not (workdir / 'Pfam-A.hmm.ssi').exists():
        subprocess.run([hmmfetch, '--index', str(link)], check=True, capture_output=True)
    dest = workdir / f'{accession}.hmm'
    with dest.open('w') as fh:
        subprocess.run([hmmfetch, str(link), resolve_accession(pfam_hmm, accession)],
                       check=True, stdout=fh, stderr=subprocess.PIPE)
    return dest


def best_pepm_per_region(profile, proteins, workdir, hmmsearch, cpus):
    """Highest-scoring PF13714 hit per region, from hmmsearch's domtblout."""
    out = workdir / 'pepm.domtbl'
    subprocess.run([hmmsearch, '--cpu', str(cpus), '--domtblout', str(out),
                    '-E', '1e-5', str(profile), str(proteins)],
                   check=True, capture_output=True)
    best = {}
    for line in out.read_text().splitlines():
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 8:
            continue
        region = int(f[0].split('|')[0])
        score = float(f[7])
        if region not in best or score > best[region][1]:
            best[region] = (f[0], score)
    return {r: name for r, (name, _) in best.items()}


def align_and_identity(names, proteins, profile, workdir, hmmalign):
    """Pairwise identity with pairwise deletion, over profile match columns."""
    from Bio import SeqIO
    keep = set(names.values())
    sub = workdir / 'pepm.faa'
    with sub.open('w') as fh:
        for rec in SeqIO.parse(proteins, 'fasta'):
            if rec.id in keep:
                fh.write(f'>{rec.id}\n{rec.seq}\n')

    sto = workdir / 'pepm.sto'
    with sto.open('w') as fh:
        subprocess.run([hmmalign, '--amino', '--trim', str(profile), str(sub)],
                       check=True, stdout=fh, stderr=subprocess.PIPE)

    aln = {}
    for line in sto.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line == '//':
            continue
        name, _, chunk = line.partition(' ')
        aln[name] = aln.get(name, '') + chunk.strip()

    order = [n for n in names.values() if n in aln]
    rows = [[ord(c) if c.isupper() else (GAP if c == '-' else None) for c in aln[n]]
            for n in order]
    cols = [k for k in range(len(rows[0])) if rows[0][k] is not None]
    mat = np.zeros((len(order), len(cols)), dtype=np.uint8)
    for r, row in enumerate(rows):
        mat[r] = [row[k] if row[k] is not None else GAP for k in cols]
    return order, mat


def link_components(order, mat, threshold, region_of):
    """Single-linkage components over pepM identity."""
    parent = {i: i for i in range(len(order))}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(order) - 1):
        a, rest = mat[i], mat[i + 1:]
        both = (a != GAP) & (rest != GAP)
        compared = both.sum(axis=1)
        matches = (both & (a == rest)).sum(axis=1)
        ok = np.where((compared > 0) &
                      (matches / np.maximum(compared, 1) >= threshold))[0]
        for k in ok:
            ra, rb = find(i), find(i + 1 + int(k))
            if ra != rb:
                parent[ra] = rb

    groups = collections.defaultdict(list)
    for i, name in enumerate(order):
        groups[find(i)].append(region_of[name])
    return sorted(groups.values(), key=len, reverse=True)


def split_oversized(parts, cap):
    """Break any partition still above the memory cap into chunks.

    A component larger than the cap cannot be split without separating BGCs
    that pepM says belong together, so this *does* risk splitting a family —
    which is why it only fires above a size that would otherwise not run at all.
    Chunks are the lesser evil against an OOM kill, and the count is reported so
    the compromise is visible rather than silent.
    """
    out, split = [], 0
    for p in parts:
        if len(p) <= cap:
            out.append(p)
            continue
        split += 1
        for i in range(0, len(p), cap):
            out.append(p[i:i + cap])
    return out, split


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--antismash_dir', type=Path, required=True)
    ap.add_argument('--pfam', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--threshold', type=float, default=0.60)
    ap.add_argument('--max_partition_size', type=int, default=20000,
                    help='hard cap; ~22,000 BGCs is the 64 GB SLURM allocation')
    ap.add_argument('--min_to_partition', type=int, default=4000,
                    help='below this, one job is cheaper than the split')
    ap.add_argument('--cpus', type=int, default=4)
    ap.add_argument('--hmmfetch', default=shutil.which('hmmfetch') or 'hmmfetch')
    ap.add_argument('--hmmsearch', default=shutil.which('hmmsearch') or 'hmmsearch')
    ap.add_argument('--hmmalign', default=shutil.which('hmmalign') or 'hmmalign')
    args = ap.parse_args()

    gbks = region_gbks(args.antismash_dir)
    print(f'{len(gbks)} region GenBanks')
    if not gbks:
        raise SystemExit('no region GBKs found')

    # Partitions are staged flat into each BiG-SCAPE job, so two regions sharing
    # a filename would silently overwrite and quietly drop a BGC. Region files
    # are named by contig accession and so are unique in practice, but fail
    # loudly rather than lose data if that ever stops holding.
    dupes = collections.Counter(g.name for g in gbks)
    clashes = [n for n, c in dupes.items() if c > 1]
    if clashes:
        raise SystemExit(
            f'{len(clashes)} region filenames are not unique, e.g. {clashes[:3]}; '
            'flat staging would drop BGCs')

    # Below the quadratic regime, partitioning costs more than it saves: the
    # fixed Pfam load per partition dominates, and one job comfortably fits.
    if len(gbks) < args.min_to_partition:
        print(f'under --min_to_partition ({args.min_to_partition}); one partition')
        with args.out.open('w') as fh:
            fh.write('partition\tgbk\n')
            for g in gbks:
                fh.write(f'0\t{g}\n')
        return 0

    work = args.out.parent / '_partition'
    work.mkdir(parents=True, exist_ok=True)
    proteins = work / 'proteins.faa'
    n_prot = extract_proteins(gbks, proteins)
    print(f'{n_prot:,} CDS translations')

    profile = fetch_profile(args.pfam, PEPM_ACCESSION, work, args.hmmfetch)
    names = best_pepm_per_region(profile, proteins, work, args.hmmsearch, args.cpus)
    print(f'pepM found in {len(names)} of {len(gbks)} regions '
          f'({100 * len(names) / len(gbks):.1f}%)')

    region_of = {name: gbks[r] for r, name in names.items()}
    order, mat = align_and_identity(names, proteins, profile, work, args.hmmalign)
    print(f'alignment {len(order)} x {mat.shape[1]} match columns')

    parts = link_components(order, mat, args.threshold, region_of)

    # Any region without a pepM hit becomes its own partition rather than being
    # dropped: it cannot be placed by pepM, and silently losing a BGC from the
    # clustering would be far worse than running it alone.
    placed = {g for p in parts for g in p}
    orphans = [g for g in gbks if g not in placed]
    parts.extend([[g] for g in orphans])
    if orphans:
        print(f'{len(orphans)} regions without a pepM hit, each its own partition')

    parts, forced = split_oversized(parts, args.max_partition_size)
    if forced:
        print(f'WARNING: {forced} component(s) exceeded --max_partition_size and were '
              f'chunked; this can split a real family')

    sizes = [len(p) for p in parts]
    print(f'{len(parts)} partitions, largest {max(sizes)} '
          f'({100 * max(sizes) / len(gbks):.1f}%), '
          f'work {sum((s / len(gbks)) ** 2 for s in sizes):.0%} of one job')

    with args.out.open('w') as fh:
        fh.write('partition\tgbk\n')
        for i, p in enumerate(parts):
            for g in p:
                fh.write(f'{i}\t{g}\n')
    shutil.rmtree(work, ignore_errors=True)
    print(f'wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
