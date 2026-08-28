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
  Transaminase  SMCOG1019          Phosphonopyruvate transaminase, PalB-like (Aminotran_1_2 /
                                   PF00155, AAT superfamily, fold type I PLP)
                                   phosphonopyruvate → L-phosphonoalanine (ref: PnaA)
                                   Co-occurs with sulfhydrylase (SMCOG1168) in all GCF-7 BGCs.
                                   Note: SMCOG1013 (Aminotran_3, fold type IV) was used here until
                                   2026-08-25. It is a different aminotransferase class from PalB
                                   and produced 5 false Transaminase calls in 6 on Pantoea.
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
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import itol
from utils.coupling_confidence import load_references, support, percent_identity
from utils.antismash_parser import (build_json_index, cds_in_segments, find_region_feature,
                                    genome_from_gbk_path, parse_bgc_label,
                                    parse_location_segments)


# ─── Coupling enzyme class definitions ──────────────────────────────────────

CLASSES = [
    # (class_id, display_label, hex_color)
    ('Synthase',                        'Synthase — phosphonomethylmalate synthase (PnPyr + AcCoA)',            '#e41a1c'),
    ('Decarboxylase',                   'Decarboxylase — phosphonopyruvate decarboxylase (ThDP-dependent)',     '#377eb8'),
    ('Decarboxylase-Nucleotidyltransferase', 'Decarboxylase-Nucleotidyltransferase — phosphonopyruvate decarboxylase + CDP-activation', '#984ea3'),
    ('Reductase',                       'Reductase — phosphonopyruvate reductase (Fe-ADH)',                    '#4daf4a'),
    ('Transaminase',                    'Transaminase — phosphonopyruvate transaminase, PalB-like Aminotran_1_2 (→ PnAla)','#ff7f00'),
    ('Unknown',                         'Unknown / not detected',                                               '#aaaaaa'),
]

CLASS_COLORS = {cid: color for cid, _, color in CLASSES}


def classify_bgc(json_path, contig_id, region_num):
    """Classify a BGC's coupling enzyme, and return the evidence behind the call.

    Returns (class_id, deciding_marker, marker_seqs) where marker_seqs maps every
    marker seen in the region to the protein sequences carrying it — including
    `PEP_mutase`, so callers can score pepM divergence without a second pass.
    """
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
        return 'Unknown', None, {}

    for rec in data['records']:
        if contig_id not in rec.get('id', ''):
            continue

        # Find the matching phosphonate region
        region_match = find_region_feature(rec, region_num, product_filter='phosphonate')
        if region_match is None:
            continue

        # Segments, not a single (start, end): an origin-spanning region is a genuine
        # join of two disjoint intervals. The first pair alone drops the second segment
        # — which is where the coupling enzyme sat in all 5 such Pantoea BGCs, leaving
        # them misclassified `Unknown` — while (min_start, max_end) invents a span
        # covering most of the replicon.
        region_segs = parse_location_segments(region_match.get('location', ''))

        # Collect biosynthetic rule hits and SMCOG annotations from CDSes in region
        rule_hits = set()
        smcog_hits = set()
        # marker -> protein sequences carrying it, so the CDS that drives the call can
        # be scored against the characterised references afterwards
        marker_seqs = {}

        for feat in rec.get('features', []):
            if feat.get('type') != 'CDS':
                continue
            # Filter to CDSes overlapping any segment of the region
            if not cds_in_segments(feat, region_segs):
                continue
            quals = feat.get('qualifiers', {})
            translation = (quals.get('translation') or [''])[0]
            for gf in quals.get('gene_functions', []):
                m = re.search(r'(SMCOG\d+)', gf)
                if m:
                    smcog_hits.add(m.group(1))
                    if translation:
                        marker_seqs.setdefault(m.group(1), []).append(translation)
            for sd in quals.get('sec_met_domain', []):
                # e.g. "Fe-ADH (E-value: ...)"
                domain_name = sd.split('(')[0].strip()
                rule_hits.add(domain_name)
                if translation:
                    marker_seqs.setdefault(domain_name, []).append(translation)
            # Also parse rule-based-clusters from gene_functions
            for gf in quals.get('gene_functions', []):
                if 'rule-based-clusters' in gf:
                    # extract domain name after the last ':'
                    parts = gf.split(':')
                    if len(parts) >= 3:
                        rule_hits.add(parts[-1].strip())

        # Classification (checked in priority order)
        if 'SMCOG1271' in smcog_hits:
            return 'Synthase', 'SMCOG1271', marker_seqs
        # Ppd-CDP before plain Ppd: both share SMCOG1055, but Ppd-CDP additionally
        # encodes cytidylyltransferase(s) (NTP_transf_3) for CDP-activation.
        has_tpp = 'TPP_enzyme_C' in rule_hits or 'TPP_enzyme_M' in rule_hits
        has_ntp = 'NTP_transf_3' in rule_hits or 'NTP_transf_2' in rule_hits
        if has_tpp and has_ntp:
            return 'Decarboxylase-Nucleotidyltransferase', 'TPP_enzyme_C', marker_seqs
        if 'SMCOG1055' in smcog_hits:
            return 'Decarboxylase', 'SMCOG1055', marker_seqs
        # Fallback: TPP_enzyme_C alone is sufficient evidence for a decarboxylase
        # coupling enzyme. Some BGCs have ThDP enzymes too divergent to score against
        # the SMCOG1055 HMM but still carry the TPP_enzyme_C domain in antiSMASH's
        # rule-based scan (e.g. GCF-1, GCF-12 singletons in Pantoea).
        if has_tpp:
            return 'Decarboxylase', 'TPP_enzyme_C', marker_seqs
        # Reductase after Ppd: some BGCs contain an unrelated Fe-ADH gene elsewhere
        # in the antiSMASH region that would mask an SMCOG1055-annotated coupling enzyme.
        # True Reductase BGCs (GCF4/9) carry Fe-ADH but no SMCOG1055.
        if 'Fe-ADH' in rule_hits:
            return 'Reductase', 'Fe-ADH', marker_seqs
        # Transaminase (PalB-like). SMCOG1019 = Aminotran_1_2 / PF00155, the AAT
        # superfamily (fold type I PLP) that PalB belongs to. This previously tested
        # SMCOG1013 (Aminotran_3, fold type IV) — a different enzyme class entirely,
        # which on Pantoea called 6 BGCs Transaminase where only 1 carries SMCOG1019.
        # Reductase is still checked first: an unrelated Fe-ADH elsewhere in the region
        # should not be overridden by an aminotransferase hit.
        if 'SMCOG1019' in smcog_hits:
            return 'Transaminase', 'SMCOG1019', marker_seqs

        return 'Unknown', None, marker_seqs

    return 'Unknown', None, {}


# ─── Build genome → antiSMASH JSON index ─────────────────────────────────────

# ─── iTOL output ─────────────────────────────────────────────────────────────

def write_support_tsv(support_rows, outpath):
    """Write per-BGC reference support.

    Advisory metadata: it records how similar the enzyme that drove each call is to the
    nearest characterised reference of every class. It does not gate or reorder
    anything. A low value is genuinely ambiguous — the protein may not belong to the
    assigned class, or it may be a novel variant unlike the one characterised example.
    Both call for manual inspection, which is the point of publishing the number.

    `n_refs` is included per class because a score against a single reference (as for
    Transaminase and Reductase) says "similar to PnaA/VlpB specifically", not "similar
    to that enzyme class".
    """
    classes = sorted({c for row in support_rows.values() for c in row[3]})
    with open(outpath, 'w') as f:
        f.write('# Reference support for coupling enzyme assignments — ADVISORY ONLY.\n')
        f.write('# The class comes from antiSMASH SMCOG/domain markers, which are broad by\n')
        f.write('# design: characterised phosphonate coupling enzymes are scarce, and a\n')
        f.write('# narrow reference-driven classifier would only recover known chemistry.\n')
        f.write('# pct_id_<class> = percent identity to the closest reference of that class.\n')
        f.write('# Empirical poles on Pantoea: 93.8-100%% orthologue, 21.8-30.8%% superfamily\n')
        f.write('# background. Low support warrants manual review, NOT automatic rejection.\n')
        header = ['bgc', 'assigned_class', 'deciding_marker', 'protein_len',
                  'assigned_pct_id', 'assigned_ref', 'assigned_n_refs', 'runner_up', 'margin',
                  'pepm_pct_id', 'pepm_ref', 'pepm_len']
        header += [f'pct_id_{c}' for c in classes]
        f.write('\t'.join(header) + '\n')
        for label in sorted(support_rows):
            cls, marker, plen, sup, pepm_pid, pepm_ref, pepm_len = support_rows[label]
            own = sup.get(cls, {})
            others = sorted(((c, v['pct_id']) for c, v in sup.items() if c != cls),
                            key=lambda x: -x[1])
            runner, runner_pid = (others[0] if others else ('-', 0.0))
            row = [label, cls, marker or '-', str(plen),
                   f"{own.get('pct_id', 0.0):.1f}", str(own.get('best_ref', '-')),
                   str(own.get('n_refs', 0)), runner,
                   f"{own.get('pct_id', 0.0) - runner_pid:.1f}",
                   f"{pepm_pid:.1f}", pepm_ref, str(pepm_len)]
            row += [f"{sup[c]['pct_id']:.1f}" for c in classes]
            f.write('\t'.join(row) + '\n')


def write_colorstrip(metadata, classifications, outpath, bgc_type):
    counts = defaultdict(int)
    for cls in classifications.values():
        counts[cls] += 1

    # Legend — only include classes that appear in the data
    present = [c for c in CLASSES if counts.get(c[0], 0) > 0]
    legend_items = [(f'{lbl} (n={counts[cid]})', color) for cid, lbl, color in present]

    itol.write_colorstrip(
        outpath, f'Coupling enzyme ({bgc_type})',
        entries=[(bgc['label'],
                  CLASS_COLORS[classifications.get(bgc['label'], 'Unknown')],
                  classifications.get(bgc['label'], 'Unknown'))
                 for bgc in metadata],
        legend=('Coupling enzyme', itol.simple_legend(legend_items)),
        color='#333333',
        options=[('STRIP_WIDTH', 40), ('SHOW_BORDER', 1), ('BORDER_WIDTH', 0.5)],
    )

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
    parser.add_argument('--reference_faa', default=None,
                        help='Characterised coupling enzyme FASTA; enables the reference-support TSV')
    parser.add_argument('--reference_pepm', default=None,
                        help='Characterised pepM FASTA; adds the pepM divergence axis')
    parser.add_argument('--bgc_type',      default='phosphonate',
                        help='BGC product type label for display (default: phosphonate)')
    args = parser.parse_args()

    with open(args.metadata) as f:
        metadata = json.load(f)

    print(f'Building antiSMASH JSON index from: {args.antismash_dir}')
    json_index = build_json_index(args.antismash_dir)
    print(f'  Found {len(json_index)} genome JSON files')

    print(f'Classifying coupling enzymes for {len(metadata)} BGCs...')
    classifications = {}
    support_rows = {}
    # Reference support is advisory metadata for manual review — see
    # utils/coupling_confidence. Absent references simply skip it.
    ref_path = args.reference_faa
    references = load_references(ref_path) if ref_path and os.path.exists(ref_path) else {}
    if not references:
        print('  (no coupling enzyme references given — skipping reference support)')
    pepm_references = []
    if args.reference_pepm and os.path.exists(args.reference_pepm):
        from Bio import SeqIO as _SeqIO
        pepm_references = [(r.id.split('|')[2] if r.id.count('|') >= 2 else r.id, str(r.seq))
                           for r in _SeqIO.parse(args.reference_pepm, 'fasta')]
    missing_json = 0

    for bgc in metadata:
        label = bgc['label']
        gbk_path = bgc.get('gbk_path', '')

        # Derive genome folder name from gbk_path
        genome = genome_from_gbk_path(gbk_path) if gbk_path else None

        # Parse contig and region from label
        contig_id, region_num = parse_bgc_label(label)

        if genome and genome in json_index:
            cls, marker, marker_seqs = classify_bgc(json_index[genome], contig_id, region_num)
        else:
            # Fall back: search all JSONs for a record matching contig_id
            cls, marker, marker_seqs = 'Unknown', None, {}
            for gen, jpath in json_index.items():
                c, mk, ms = classify_bgc(jpath, contig_id, region_num)
                if c != 'Unknown':
                    cls, marker, marker_seqs = c, mk, ms
                    break
            if cls == 'Unknown':
                missing_json += 1

        classifications[label] = cls
        if references and marker and marker_seqs.get(marker):
            # Score the CDS that drove the call against every characterised class.
            # Advisory only — nothing here changes `cls`.
            seq = max(marker_seqs[marker], key=len)
            # pepM is the hallmark gene, present in every phosphonate BGC, so its
            # divergence is the natural second axis: how unusual is the scaffold gene,
            # independently of how unusual the coupling chemistry is.
            pepm_seqs = marker_seqs.get('PEP_mutase') or marker_seqs.get('PEP_mutase_1') or []
            pepm_pid, pepm_ref, pepm_len = 0.0, '-', 0
            if pepm_seqs and pepm_references:
                pseq = max(pepm_seqs, key=len)
                pepm_len = len(pseq)
                best = max(((percent_identity(pseq, rs), rn) for rn, rs in pepm_references),
                           default=(0.0, '-'))
                pepm_pid, pepm_ref = round(best[0], 1), best[1]
            support_rows[label] = (cls, marker, len(seq), support(seq, references),
                                   pepm_pid, pepm_ref, pepm_len)

    if missing_json:
        print(f'  Warning: {missing_json} BGCs could not be matched to an antiSMASH JSON')

    if support_rows:
        conf_path = os.path.splitext(args.outfile)[0].replace('_itol_coupling', '') + '_coupling_support.tsv'
        print(f'Writing reference support to: {conf_path}')
        write_support_tsv(support_rows, conf_path)

    print(f'Writing iTOL annotation to: {args.outfile}')
    os.makedirs(os.path.dirname(args.outfile) or '.', exist_ok=True)
    write_colorstrip(metadata, classifications, args.outfile, args.bgc_type)
    print('\nDone. Upload this file to iTOL alongside the .nwk tree.')


if __name__ == '__main__':
    main()
