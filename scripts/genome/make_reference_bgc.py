#!/usr/bin/env python3
"""Turn a curated cluster deposit into a BiG-SCAPE reference BGC.

BiG-SCAPE's `--reference-dir` requires antiSMASH-processed GenBank: its parser reads a
`region` feature and raises InvalidGBKError without one. A raw NCBI deposit has no such
feature, and BiG-SCAPE then **drops the file silently** -- no warning, no error, it
simply does not appear among the loaded records. That is how a reference set can look
installed and contribute nothing.

Running antiSMASH would be the obvious fix and does not work here. Four characterised
phosphonate clusters -- argolaphos, bialaphos, phosphinothricin tripeptide and
phosphonothrixin -- produce ZERO regions at every strictness setting, including `loose`.
They all carry a pepM (43-70% identity to the curated references), so they are
unambiguously phosphonate clusters; antiSMASH's phosphonate rule needs pepM *plus* a
partner domain from a fixed list, and these do not supply one.

Since these deposits are curated cluster sequences from the literature -- the whole
record IS the cluster -- declaring one region spanning it is a statement of fact rather
than a guess. That is what this script writes.

Usage:
    python scripts/genome/make_reference_bgc.py --in argolaphos.gb \\
        --out assets/phosphonate_reference_bgcs/argolaphos.region001.gbk \\
        --product phosphonate
"""
import argparse
import sys
from pathlib import Path


TOOL = 'make_reference_bgc.py'

RULE = 'curated reference; no antiSMASH detection rule applies'

EDGE_NOTE = (
    'contig_edge is set True deliberately, and is NOT a claim that this record runs '
    'off a contig: it is a curated cluster, complete by construction. It is set so '
    'BiG-SCAPE\'s auto alignment mode compares this reference against antiSMASH '
    'query regions by their shared part (LCS + extension) rather than end to end. '
    'antiSMASH regions are a rule core plus a fixed neighbourhood, so they carry '
    'flanking DNA the cluster does not and can stop before the cluster ends; '
    'compared end to end, a curated reference scores far from its own genome\'s '
    'region -- measured at 0.364 for pantaphos against LMG 5342 itself, which falls '
    'to 0.000 with this set. Unrelated clusters are unaffected (0.946 either way).')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in', dest='infile', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--product', default='phosphonate',
                    help='product label BiG-SCAPE records for the region')
    ap.add_argument('--note', default='',
                    help='provenance note written into the region feature')
    ap.add_argument('--no-contig-edge', dest='contig_edge', action='store_false',
                    help='compare end to end instead of by the shared part; see the '
                         'note written into the region feature')
    ap.set_defaults(contig_edge=True)
    a = ap.parse_args()

    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    records = list(SeqIO.parse(str(a.infile), 'genbank'))
    if not records:
        sys.exit(f'{a.infile}: no GenBank records')
    if len(records) > 1:
        sys.exit(f'{a.infile}: {len(records)} records; a reference BGC must be one '
                 f'contiguous cluster, so split it first')
    rec = records[0]

    cds = [f for f in rec.features if f.type == 'CDS']
    if not cds:
        sys.exit(f'{a.infile}: no CDS features. BiG-SCAPE needs genes, and a deposit '
                 f'without them cannot act as a reference.')
    if any(f.type == 'region' for f in rec.features):
        sys.exit(f'{a.infile}: already has a region feature — use it directly.')

    # Span the coding extent rather than the whole record: trailing sequence with no
    # genes contributes nothing and would misreport the cluster's length.
    lo = min(int(f.location.start) for f in cds)
    hi = max(int(f.location.end) for f in cds)

    quals = {
        'region_number': ['1'],
        # BiG-SCAPE's AS5 reader requires this even though it carries no information
        # for a single-region record; without it the file is rejected outright.
        'candidate_cluster_numbers': ['1'],
        'contig_edge': [str(a.contig_edge)],
        'product': [a.product],
        'tool': [TOOL],
        'rules': [RULE],
        'note': [EDGE_NOTE if a.contig_edge else
                 'contig_edge False: compared end to end against query regions.'],
    }
    if a.note:
        quals['note'].append(a.note)
    # antiSMASH writes a four-level hierarchy and BiG-SCAPE's AS5 reader walks all of
    # it: region -> cand_cluster -> protocluster -> proto_core. Supplying only the
    # region fails, then only the cand_cluster fails, each with its own error — so the
    # whole chain is written here rather than discovered one rejection at a time.
    loc = FeatureLocation(lo, hi)
    feats = [
        SeqFeature(loc, type='region', qualifiers=quals),
        SeqFeature(loc, type='cand_cluster', qualifiers={
            'candidate_cluster_number': ['1'],
            'contig_edge': [str(a.contig_edge)],
            'product': [a.product],
            'kind': ['single'],
            'protoclusters': ['1'],
            'detection_rules': [RULE],
            'tool': [TOOL],
        }),
        SeqFeature(loc, type='protocluster', qualifiers={
            'protocluster_number': ['1'],
            'contig_edge': [str(a.contig_edge)],
            'product': [a.product],
            'aStool': ['rule-based-clusters'],
            'category': ['other'],
            'core_location': [f'[{lo}:{hi}](+)'],
            'cutoff': ['0'],
            'neighbourhood': ['0'],
            'detection_rule': [RULE],
            'tool': [TOOL],
        }),
        SeqFeature(loc, type='proto_core', qualifiers={
            'protocluster_number': ['1'],
            'product': [a.product],
            'aStool': ['rule-based-clusters'],
            'cutoff': ['0'],
            'neighbourhood': ['0'],
            'detection_rule': [RULE],
            'tool': [TOOL],
        }),
    ]
    for i, f in enumerate(feats):
        rec.features.insert(i, f)

    # BiG-SCAPE picks its parser from structured_comment['antiSMASH-Data']['Version']
    # and falls back to the antiSMASH-4 reader, which looks for `cluster` features
    # rather than `region`, when that key is missing. Without this the file is rejected
    # as "does not contain an antiSMASH cluster or region feature" despite having a
    # perfectly good region.
    #
    # The version declares the FILE FORMAT, not that antiSMASH ran — it did not, and
    # the NOTE beside it says so, so a human reading the deposit is not misled.
    sc = rec.annotations.setdefault('structured_comment', {})
    sc['antiSMASH-Data'] = {
        'Version': '5.0',
        'NOTE': ('Region feature declared by scripts/genome/make_reference_bgc.py. '
                 'antiSMASH did NOT run on this record and detects no region in it; '
                 'the deposit is a curated cluster sequence, so the whole record is '
                 'the region. Version above denotes the AS5 file format BiG-SCAPE '
                 'parses, not a tool that produced this file.'),
    }

    a.out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([rec], str(a.out), 'genbank')
    print(f'{a.infile.name}: {len(cds)} CDS, region {lo:,}-{hi:,} '
          f'({(hi-lo)/1000:.1f} kb) -> {a.out}')


if __name__ == '__main__':
    main()
