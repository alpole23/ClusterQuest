#!/usr/bin/env python3
"""
Generate an iTOL coupling-enzyme colorstrip for phosphonate BGC trees.

Reads antiSMASH JSON files and BiG-SCAPE DB to classify each phosphonate BGC
by the coupling enzyme acting on phosphonopyruvate (the step immediately
downstream of PEP mutase). Outputs an iTOL DATASET_COLORSTRIP file compatible
with any tree whose leaf labels are {contig}.region{NNN} (as produced by
bgc_pfam_tree.py or bgc_synteny_tree.py).

Coupling enzyme classes detected (checked in priority order):
  Synthase      SMCOG1271          Phosphonomethylmalate synthase (HMGL superfamily)
                                   phosphonopyruvate + acetyl-CoA → phosphonomethylmalate
                                   → phosphinothricin-type products (refs: FrbC, HvrC)
  Ppd-CDP       SMCOG1055          Same ThDP decarboxylation as Ppd, but BGC additionally
                 + NTP_transf_3    encodes cytidylyltransferase(s) (NTP_transf_3) for
                                   CDP-activation → phosphonolipid pathway.
                                   Checked before plain Ppd because both share SMCOG1055.
  Ppd           SMCOG1055          Phosphonopyruvate decarboxylase (ThDP-dependent)
                                   → 2-phosphonoacetaldehyde; BGC lacks cytidylyltransferase.
                                   Checked before Reductase: some BGCs contain an unrelated
                                   Fe-ADH gene elsewhere in the region that would otherwise
                                   mask the SMCOG1055-annotated coupling enzyme (e.g. GCF11).
  Reductase     Fe-ADH rule        Phosphonopyruvate reductase (iron-containing ADH)
                                   → phosphonolactate (ref: VlpB)
  Transaminase  SMCOG1013          Phosphonopyruvate transaminase (Aminotran_3, fold type IV PLP)
                                   phosphonopyruvate → L-phosphonoalanine (ref: PnaA)
                                   Co-occurs with sulfhydrylase (SMCOG1168) in all GCF-7 BGCs.
                                   Note: SMCOG1013 also appears downstream in Reductase clusters
                                   (GCF-4); Reductase is checked first to avoid false positives.
  Unknown       —                  No coupling enzyme identified

Usage:
    python scripts/bgc_coupling_annotation.py \\
        --antismash_dir results/antismash_results/Pantoea \\
        --metadata results/bgc_trees/Pantoea/phosphonate_metadata.json \\
        --outfile results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \\
        [--bgc_type phosphonate]

The metadata JSON must be the output of bgc_pfam_tree.py or bgc_synteny_tree.py
(contains 'label' and 'gbk_path' for each BGC).
"""

import argparse
import json
import os
import re
from collections import defaultdict


# ─── Coupling enzyme class definitions ──────────────────────────────────────

CLASSES = [
    # (class_id, display_label, hex_color)
    ('Synthase',                        'Synthase — phosphonomethylmalate synthase (PnPyr + AcCoA)',            '#e41a1c'),
    ('Decarboxylase',                   'Decarboxylase — phosphonopyruvate decarboxylase (ThDP-dependent)',     '#377eb8'),
    ('Decarboxylase-Nucleotidyltransferase', 'Decarboxylase-Nucleotidyltransferase — phosphonopyruvate decarboxylase + CDP-activation', '#984ea3'),
    ('Reductase',                       'Reductase — phosphonopyruvate reductase (Fe-ADH)',                    '#4daf4a'),
    ('Transaminase',                    'Transaminase — phosphonopyruvate transaminase, Aminotran_3 (→ PnAla)','#ff7f00'),
    ('Unknown',                         'Unknown / not detected',                                               '#aaaaaa'),
]

CLASS_COLORS = {cid: color for cid, _, color in CLASSES}


def classify_bgc(json_path, contig_id, region_num):
    """
    Open an antiSMASH JSON, find the record matching contig_id and region_num,
    and classify the coupling enzyme based on gene_functions and sec_met_domain
    annotations of the CDSes within the region.

    Returns a class_id string.
    """
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception:
        return 'Unknown'

    for rec in data['records']:
        if contig_id not in rec.get('id', ''):
            continue

        # Find the matching phosphonate region
        region_match = None
        for feat in rec.get('features', []):
            if feat.get('type') == 'region':
                rnum = str(feat.get('qualifiers', {}).get('region_number', ['?'])[0]).zfill(3)
                products = feat.get('qualifiers', {}).get('product', [])
                if rnum == str(region_num).zfill(3) and any('phosphonate' in p for p in products):
                    region_match = feat
                    break
        if region_match is None:
            continue

        # Parse region boundaries so we only inspect CDSes within this region
        region_loc = region_match.get('location', '')
        region_m = re.search(r'\[(\d+):(\d+)\]', region_loc)
        region_start = int(region_m.group(1)) if region_m else 0
        region_end   = int(region_m.group(2)) if region_m else float('inf')

        # Collect biosynthetic rule hits and SMCOG annotations from CDSes in region
        rule_hits = set()
        smcog_hits = set()

        for feat in rec.get('features', []):
            if feat.get('type') != 'CDS':
                continue
            # Filter to CDSes within the region boundaries
            cds_loc = feat.get('location', '')
            cds_m = re.search(r'\[(\d+):(\d+)\]', cds_loc)
            if cds_m:
                cds_start = int(cds_m.group(1))
                cds_end   = int(cds_m.group(2))
                if cds_end <= region_start or cds_start >= region_end:
                    continue
            quals = feat.get('qualifiers', {})
            for gf in quals.get('gene_functions', []):
                m = re.search(r'(SMCOG\d+)', gf)
                if m:
                    smcog_hits.add(m.group(1))
            for sd in quals.get('sec_met_domain', []):
                # e.g. "Fe-ADH (E-value: ...)"
                domain_name = sd.split('(')[0].strip()
                rule_hits.add(domain_name)
            # Also parse rule-based-clusters from gene_functions
            for gf in quals.get('gene_functions', []):
                if 'rule-based-clusters' in gf:
                    # extract domain name after the last ':'
                    parts = gf.split(':')
                    if len(parts) >= 3:
                        rule_hits.add(parts[-1].strip())

        # Classification (checked in priority order)
        if 'SMCOG1271' in smcog_hits:
            return 'Synthase'
        # Ppd-CDP before plain Ppd: both share SMCOG1055, but Ppd-CDP additionally
        # encodes cytidylyltransferase(s) (NTP_transf_3) for CDP-activation.
        has_tpp = 'TPP_enzyme_C' in rule_hits or 'TPP_enzyme_M' in rule_hits
        has_ntp = 'NTP_transf_3' in rule_hits or 'NTP_transf_2' in rule_hits
        if has_tpp and has_ntp:
            return 'Decarboxylase-Nucleotidyltransferase'
        if 'SMCOG1055' in smcog_hits:
            return 'Decarboxylase'
        # Fallback: TPP_enzyme_C alone is sufficient evidence for a decarboxylase
        # coupling enzyme. Some BGCs have ThDP enzymes too divergent to score against
        # the SMCOG1055 HMM but still carry the TPP_enzyme_C domain in antiSMASH's
        # rule-based scan (e.g. GCF-1, GCF-12 singletons in Pantoea).
        if has_tpp:
            return 'Decarboxylase'
        # Reductase after Ppd: some BGCs contain an unrelated Fe-ADH gene elsewhere
        # in the antiSMASH region that would mask an SMCOG1055-annotated coupling enzyme.
        # True Reductase BGCs (GCF4/9) carry Fe-ADH but no SMCOG1055.
        if 'Fe-ADH' in rule_hits:
            return 'Reductase'
        # Transaminase: SMCOG1013 also appears downstream in Reductase clusters (GCF-4),
        # so Reductase is checked first.
        if 'SMCOG1013' in smcog_hits:
            return 'Transaminase'

        return 'Unknown'

    return 'Unknown'


# ─── Build genome → antiSMASH JSON index ─────────────────────────────────────

def build_json_index(antismash_dir, bgc_type='phosphonate'):
    """
    Walk antismash_dir and return a dict:
        genome_name → json_path
    where genome_name is the subdirectory name (e.g. 'Pantoea_ananatis_LMG2665').
    """
    index = {}
    for genome in os.listdir(antismash_dir):
        json_path = os.path.join(antismash_dir, genome, f'{genome}.json')
        if os.path.exists(json_path):
            index[genome] = json_path
    return index


def genome_from_gbk_path(gbk_path):
    """
    Extract the genome folder name from a gbk_path stored in metadata.
    Paths look like:
      /…/antismash_input/Pantoea_ananatis_LMG2665/CONTIG.region001.gbk
    or results/antismash_results/Pantoea/Pantoea_ananatis_LMG2665/…
    """
    # Walk up from the .gbk file — genome is the directory one level up
    return os.path.basename(os.path.dirname(gbk_path))


def parse_bgc_label(label):
    """
    Extract contig_id and region_number from a BGC label like:
      JBBJSA010000012.1.region001  → ('JBBJSA010000012.1', '001')
      JBBJSA010000012.1.region001_1 (duplicate suffix) → strip suffix first
    """
    # Strip numeric duplicate suffix added by make_labels_unique
    label = re.sub(r'_\d+$', '', label)
    m = re.search(r'^(.+?)\.region(\d+)$', label)
    if m:
        return m.group(1), m.group(2)
    return label, '001'


# ─── iTOL output ─────────────────────────────────────────────────────────────

def write_colorstrip(metadata, classifications, outpath, bgc_type):
    counts = defaultdict(int)
    for cls in classifications.values():
        counts[cls] += 1

    with open(outpath, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tCoupling enzyme ({bgc_type})\n')
        f.write('COLOR\t#333333\n')
        f.write('STRIP_WIDTH\t40\n')
        f.write('SHOW_BORDER\t1\n')
        f.write('BORDER_WIDTH\t0.5\n')

        # Legend — only include classes that appear in the data
        present = [c for c in CLASSES if counts.get(c[0], 0) > 0]
        legend_shapes  = '\t'.join('1' for _ in present)
        legend_colors  = '\t'.join(c[2] for c in present)
        legend_labels  = '\t'.join(
            f'{c[1]} (n={counts[c[0]]})' for c in present
        )
        f.write(f'LEGEND_TITLE\tCoupling enzyme\n')
        f.write(f'LEGEND_SHAPES\t{legend_shapes}\n')
        f.write(f'LEGEND_COLORS\t{legend_colors}\n')
        f.write(f'LEGEND_LABELS\t{legend_labels}\n')
        f.write('DATA\n')

        for bgc in metadata:
            lbl = bgc['label']
            cls = classifications.get(lbl, 'Unknown')
            color = CLASS_COLORS[cls]
            f.write(f'{lbl}\t{color}\t{cls}\n')

    n_classified = sum(1 for c in classifications.values() if c != 'Unknown')
    print(f'  Coupling enzyme strip: {outpath}')
    print(f'  Classified: {n_classified}/{len(metadata)} BGCs')
    for cid, lbl, _ in CLASSES:
        if counts[cid]:
            print(f'    {cid:<20} {counts[cid]:>4}')


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Generate iTOL coupling-enzyme colorstrip for phosphonate BGC trees',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--antismash_dir', required=True,
                        help='antiSMASH results directory (contains one subdir per genome)')
    parser.add_argument('--metadata',      required=True,
                        help='BGC metadata JSON from bgc_pfam_tree.py or bgc_synteny_tree.py')
    parser.add_argument('--outfile',       required=True,
                        help='Output iTOL colorstrip file path')
    parser.add_argument('--bgc_type',      default='phosphonate',
                        help='BGC product type label for display (default: phosphonate)')
    args = parser.parse_args()

    with open(args.metadata) as f:
        metadata = json.load(f)

    print(f'Building antiSMASH JSON index from: {args.antismash_dir}')
    json_index = build_json_index(args.antismash_dir, args.bgc_type)
    print(f'  Found {len(json_index)} genome JSON files')

    print(f'Classifying coupling enzymes for {len(metadata)} BGCs...')
    classifications = {}
    missing_json = 0

    for bgc in metadata:
        label = bgc['label']
        gbk_path = bgc.get('gbk_path', '')

        # Derive genome folder name from gbk_path
        genome = genome_from_gbk_path(gbk_path) if gbk_path else None

        # Parse contig and region from label
        contig_id, region_num = parse_bgc_label(label)

        if genome and genome in json_index:
            cls = classify_bgc(json_index[genome], contig_id, region_num)
        else:
            # Fall back: search all JSONs for a record matching contig_id
            cls = 'Unknown'
            for gen, jpath in json_index.items():
                c = classify_bgc(jpath, contig_id, region_num)
                if c != 'Unknown':
                    cls = c
                    break
            if cls == 'Unknown':
                missing_json += 1

        classifications[label] = cls

    if missing_json:
        print(f'  Warning: {missing_json} BGCs could not be matched to an antiSMASH JSON')

    print(f'Writing iTOL annotation to: {args.outfile}')
    os.makedirs(os.path.dirname(args.outfile) or '.', exist_ok=True)
    write_colorstrip(metadata, classifications, args.outfile, args.bgc_type)
    print('\nDone. Upload this file to iTOL alongside the .nwk tree.')


if __name__ == '__main__':
    main()
