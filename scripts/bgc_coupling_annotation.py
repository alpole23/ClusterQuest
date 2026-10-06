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
  Decarboxylase SMCOG1055          Phosphonopyruvate decarboxylase (ThDP-dependent)
                or TPP_enzyme_C    → 2-phosphonoacetaldehyde.
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
import collections
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import itol
from utils.coupling_confidence import (load_references, support, percent_identity,
                                        BACKGROUND_CEILING_PCT)
from utils.antismash_parser import (build_json_index, cds_in_segments, find_region_feature,
                                    genome_from_gbk_path, parse_bgc_label,
                                    parse_location_segments)


# ─── Coupling enzyme class definitions ──────────────────────────────────────

CLASSES = [
    # (class_id, display_label, hex_color)
    ('Synthase',                        'Synthase — phosphonomethylmalate synthase (PnPyr + AcCoA)',            '#e41a1c'),
    ('Decarboxylase',                   'Decarboxylase — phosphonopyruvate decarboxylase (ThDP-dependent)',     '#377eb8'),
    ('Reductase',                       'Reductase — phosphonopyruvate reductase (Fe-ADH)',                    '#4daf4a'),
    ('Transaminase',                    'Transaminase — phosphonopyruvate transaminase, PalB-like Aminotran_1_2 (→ PnAla)','#ff7f00'),
    ('Unknown',                         'Unknown / not detected',                                               '#aaaaaa'),
]

CLASS_COLORS = {cid: color for cid, _, color in CLASSES}


# Vocabulary per class, for reading the product name. NOT a keyword classifier:
# an earlier attempt to classify on product text inverted on both lab-confirmed
# clusters, because "serine hydroxymethyltransferase" matches `methyltransferase`.
# This is used only two ways, both narrow.
#
# POSITIVE, and only for Decarboxylase: 95% of the 917 TPP-carrying CDS in the
# Enterobacterales run are annotated "phosphonopyruvate decarboxylase" outright.
# PGAP commits to the substrate for that enzyme, so the name is real evidence.
#
# NEGATIVE, for every class: a name that specifies a DIFFERENT reaction is a veto.
# GCF-37's "decarboxylase" is annotated 3D-(3,5/4)-trihydroxycyclohexane-1,2-dione
# acylhydrolase -- an inositol-catabolism enzyme that has no business being called
# the coupling step. 7 of 917 Decarboxylase candidates are of that kind.
#
# A vague name is NOT evidence either way, and must not be penalised: phosphonomethyl-
# malate synthase has no PGAP name at all, so its real instances read "homocitrate
# synthase/isopropylmalate synthase family protein" or "beta/alpha barrel domain-
# containing protein". Penalising vagueness would reject the true synthases.
CLASS_VOCAB = {
    'Synthase':      ('homocitrate', 'isopropylmalate', 'malate synthase', 'hmgl',
                      'citramalate'),
    'Decarboxylase': ('decarboxylase', 'thiamine pyrophosphate', 'thdp', 'pyruvate'),
    'Reductase':     ('alcohol dehydrogenase', 'reductase', 'fe-adh', 'dehydrogenase'),
    'Transaminase':  ('transaminase', 'aminotransferase'),
}
VAGUE_NAME = ('domain-containing protein', 'hypothetical', 'uncharacteri', 'duf',
              'family protein', 'putative', 'unknown', 'barrel')


def name_verdict(product, cls):
    """'names-substrate' | 'consistent' | 'vague' | 'contradicts' for one product."""
    p = (product or '').lower()
    if not p:
        return 'vague'
    if 'phosphono' in p:
        return 'names-substrate'
    if any(w in p for w in CLASS_VOCAB.get(cls, ())):
        return 'consistent'
    if any(w in p for w in VAGUE_NAME):
        return 'vague'
    return 'contradicts'


# Below this many points the top two coupling candidates are indistinguishable
# and the legacy priority order decides instead. Measured on Enterobacterales:
# of 253 regions carrying two scoreable candidates the margin is a median 13.7
# points and only 2 fall under 5, so this is a guard against a coin-flip rather
# than a knob that moves calls. Flagged in the support file either way, because a
# close call is exactly what someone should look at by hand.
AMBIGUOUS_MARGIN = 5.0


def classify_bgc(json_path, contig_id, region_num):
    """Classify a BGC's coupling enzyme, and return the evidence behind the call.

    Returns (class_id, deciding_marker, marker_seqs, candidates).

    `marker_seqs` maps every marker seen in the region to the protein sequences
    carrying it, including `PEP_mutase`, so callers can score pepM divergence
    without a second pass. `candidates` is every (class, marker) the region has
    evidence for, in the legacy priority order; the class returned here is the
    first of them, and the caller re-decides on reference identity where it can.
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
        # Per CDS, so a candidate is one GENE rather than one marker, and so the
        # cascade can read its product name and its distance to pepM.
        cds_markers = collections.defaultdict(set)
        cds_product, cds_pos, cds_seq = {}, {}, {}
        pepm_pos = None

        for feat in rec.get('features', []):
            if feat.get('type') != 'CDS':
                continue
            # Filter to CDSes overlapping any segment of the region
            if not cds_in_segments(feat, region_segs):
                continue
            quals = feat.get('qualifiers', {})
            translation = (quals.get('translation') or [''])[0]
            tag = (quals.get('locus_tag') or ['?'])[0]
            cds_product[tag] = (quals.get('product') or [''])[0]
            if translation:
                cds_seq[tag] = translation
            _segs = parse_location_segments(feat.get('location', ''))
            if _segs:
                cds_pos[tag] = (min(a for a, _ in _segs) + max(b for _, b in _segs)) // 2
            for gf in quals.get('gene_functions', []):
                m = re.search(r'(SMCOG\d+)', gf)
                if m:
                    smcog_hits.add(m.group(1))
                    cds_markers[tag].add(m.group(1))
                    if translation:
                        marker_seqs.setdefault(m.group(1), []).append(translation)
            for sd in quals.get('sec_met_domain', []):
                # e.g. "Fe-ADH (E-value: ...)"
                domain_name = sd.split('(')[0].strip()
                rule_hits.add(domain_name)
                cds_markers[tag].add(domain_name)
                if domain_name in ('PEP_mutase', 'PEP_mutase_1'):
                    pepm_pos = cds_pos.get(tag)
                if translation:
                    marker_seqs.setdefault(domain_name, []).append(translation)
            # Also parse rule-based-clusters from gene_functions
            for gf in quals.get('gene_functions', []):
                if 'rule-based-clusters' in gf:
                    # extract domain name after the last ':'
                    parts = gf.split(':')
                    if len(parts) >= 3:
                        rule_hits.add(parts[-1].strip())
                        cds_markers[tag].add(parts[-1].strip())
                        if parts[-1].strip() in ('PEP_mutase', 'PEP_mutase_1'):
                            pepm_pos = cds_pos.get(tag)

        # Every class the region has evidence for, with the marker that found it.
        # Returning only the first used to discard the rest: 268 of 1,303 regions
        # here carry two or more, so "the coupling enzyme" was being decided by a
        # hand-ordered list with nothing recorded about what it ruled out.
        # One candidate per CDS, not per marker. A protein carrying both
        # TPP_enzyme_C and Aminotran_1_2 is ONE enzyme and one candidate: counting
        # it twice made 37 regions look like a Decarboxylase/Transaminase conflict
        # when the two "candidates" were the same 2.4 kb gene 47 bp downstream of
        # pepM, annotated phosphonopyruvate decarboxylase.
        MARKER_CLASS = [('SMCOG1271', 'Synthase'), ('SMCOG1055', 'Decarboxylase'),
                        ('TPP_enzyme_C', 'Decarboxylase'), ('TPP_enzyme_M', 'Decarboxylase'),
                        ('Fe-ADH', 'Reductase'), ('SMCOG1019', 'Transaminase')]
        cands = []
        for tag, marks in cds_markers.items():
            classes = [(mk, cl) for mk, cl in MARKER_CLASS if mk in marks]
            if not classes:
                continue
            # A CDS with two markers takes the class its product name supports;
            # failing that, marker precedence within the CDS.
            prod = cds_product.get(tag, '')
            named = [(mk, cl) for mk, cl in classes
                     if name_verdict(prod, cl) in ('names-substrate', 'consistent')]
            marker, cls = (named or classes)[0]
            cands.append({'cls': cls, 'marker': marker, 'tag': tag, 'product': prod,
                          'pos': cds_pos.get(tag), 'seq': cds_seq.get(tag, '')})
        for c in cands:
            c['kb'] = (abs(c['pos'] - pepm_pos) / 1000.0
                       if c['pos'] is not None and pepm_pos is not None else None)
        order = [c for _, c in MARKER_CLASS]
        cands.sort(key=lambda c: order.index(c['cls']))
        if not cands:
            return 'Unknown', None, marker_seqs, []
        return cands[0]['cls'], cands[0]['marker'], marker_seqs, cands

    return 'Unknown', None, {}, []


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
        f.write('# decided_by names the signal that chose the class; evidence is a separate\n')
        f.write('# check on the chosen call, and "weak" means its identity is at background.\n')
        header = ['bgc', 'assigned_class', 'deciding_marker', 'protein_len',
                  'assigned_pct_id', 'assigned_ref', 'assigned_ref_organism',
                  'assigned_n_refs', 'runner_up', 'margin',
                  'candidate_classes', 'kb_from_pepM', 'conservation', 'decided_by',
                  'evidence',
                  'pepm_pct_id', 'pepm_ref', 'pepm_len']
        header += [f'pct_id_{c}' for c in classes]
        f.write('\t'.join(header) + '\n')
        for label in sorted(support_rows):
            (cls, marker, plen, sup, pepm_pid, pepm_ref, pepm_len,
             amb) = support_rows[label]
            decided_by, amb_classes, kb, cons = amb
            # Two independent columns, deliberately. `decided_by` names the signal
            # that chose the class; `evidence` says whether the chosen call has any
            # reference support. A single-candidate call is decided trivially and
            # can still have no evidence behind it.
            verdict = decided_by
            margin_txt = f'{kb:.1f}' if kb is not None else '-'
            cons_txt = f'{cons:.2f}' if cons is not None else '-'
            # Reference identity as a confidence check on EVERY call, including the
            # single-candidate ones. Those never face a tie-break, so without this
            # the one signal that can say "this is not a coupling enzyme" is never
            # consulted for them -- and 172 of them score at or below the
            # superfamily background, 168 being Reductase at 18.5% to VlpB.
            #
            # It flags, it does not overturn. With 1-3 references per class a low
            # score cannot separate "wrong class" from "novel variant unlike the one
            # characterised example", and Reductase has exactly one reference, so a
            # genuine Enterobacterales enzyme has nothing close to score against.
            own = sup.get(cls, {})
            own_pid = float(own.get('pct_id', 0.0))
            evidence = ('weak - at superfamily background'
                        if own_pid <= BACKGROUND_CEILING_PCT else 'ok')
            others = sorted(((c, v['pct_id']) for c, v in sup.items() if c != cls),
                            key=lambda x: -x[1])
            runner, runner_pid = (others[0] if others else ('-', 0.0))
            row = [label, cls, marker or '-', str(plen),
                   f"{own.get('pct_id', 0.0):.1f}", str(own.get('best_ref', '-')),
                   str(own.get('best_ref_organism', '-') or '-'),
                   str(own.get('n_refs', 0)), runner,
                   f"{own.get('pct_id', 0.0) - runner_pid:.1f}",
                   '+'.join(amb_classes) or '-', margin_txt, cons_txt, verdict,
                   evidence,
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
    pending = {}
    by_family = collections.defaultdict(list)
    fid_of = {}
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
            cls, marker, marker_seqs, cands = classify_bgc(
                json_index[genome], contig_id, region_num)
        else:
            # Fall back: search all JSONs for a record matching contig_id
            cls, marker, marker_seqs, cands = 'Unknown', None, {}, []
            for gen, jpath in json_index.items():
                c, mk, ms, cd = classify_bgc(jpath, contig_id, region_num)
                if c != 'Unknown':
                    cls, marker, marker_seqs, cands = c, mk, ms, cd
                    break
            if cls == 'Unknown':
                missing_json += 1

        fams = bgc.get('families') or []
        fid = fams[0].get('family_id') if fams else None
        fid_of[label] = fid
        by_family[fid].append(label)
        pending[label] = dict(cls=cls, marker=marker, seqs=marker_seqs, cands=cands,
                              fid=fid, ncds=bgc.get('n_domains', 0))
        classifications[label] = cls
    # ── Pass 2: rank each BGC's candidates ──────────────────────────────────
    #
    # The order is conservation, then distance, then identity, then marker
    # precedence, with a name veto applied first. Each step is there because the
    # one before it cannot decide every case, and each was calibrated on this run:
    #
    #   name veto     a product naming a DIFFERENT specific reaction is not the
    #                 coupling enzyme. GCF-37's "decarboxylase" is annotated
    #                 3D-(3,5/4)-trihydroxycyclohexane-1,2-dione acylhydrolase.
    #                 A vague name is never penalised -- phosphonomethylmalate
    #                 synthase has no PGAP name, so its real instances read
    #                 "beta/alpha barrel domain-containing protein".
    #   conservation  fraction of COMPLETE family members carrying the class.
    #                 Fragments are excluded: a truncated region that stops after
    #                 pepM has no information about the coupling enzyme, and
    #                 counting it as a dissenting vote counts ignorance as
    #                 disagreement. Ties on about half of multi-candidate regions.
    #   distance      to pepM. In 1,016 regions carrying exactly one candidate the
    #                 coupling enzyme is a median 1.2 kb away and 100% are within
    #                 5 kb, so this is tightly calibrated. Breaks most ties.
    #   identity      to the class's own characterised references. Last of the
    #                 measured signals because the reference set is 1-3 proteins
    #                 per class, so a class can lose for reasons about the
    #                 reference set rather than the protein.
    #   precedence    the legacy marker order, and the call is flagged.
    med_cds = {}
    for fid, labels in by_family.items():
        n = sorted(len(pending[l]['cands']) and pending[l].get('ncds', 0) for l in labels)
        med_cds[fid] = n[len(n) // 2] if n else 0
    conservation = collections.defaultdict(lambda: collections.defaultdict(float))
    for fid, labels in by_family.items():
        full = [l for l in labels if pending[l].get('ncds', 0) >= 0.6 * med_cds[fid]] or labels
        for c in ('Synthase', 'Decarboxylase', 'Reductase', 'Transaminase'):
            hit = sum(1 for l in full if any(k['cls'] == c for k in pending[l]['cands']))
            if hit:
                conservation[fid][c] = hit / len(full)

    for label, st in pending.items():
        cands, fid = st['cands'], st['fid']
        decided_by = 'single candidate' if len(cands) < 2 else None
        if len(cands) > 1:
            live = [c for c in cands
                    if name_verdict(c['product'], c['cls']) != 'contradicts'] or cands
            if len(live) < len(cands):
                decided_by = 'product name (veto)'
            # A product that NAMES THE PHOSPHONATE SUBSTRATE outranks one that does
            # not. This is the step that was missing, and a lab-characterised
            # cluster caught it: Winslowiella iniecta B149 has its explicitly
            # annotated "phosphonopyruvate decarboxylase" 9.8 kb from pepM and a
            # generic "pyridoxal phosphate-dependent aminotransferase" at 1.1 kb.
            # Distance picks the aminotransferase and is wrong. 95% of the 917
            # TPP-carrying CDS in this run are named for the substrate, so when one
            # candidate has that name and another does not, the name decides.
            named = [c for c in live
                     if name_verdict(c['product'], c['cls']) == 'names-substrate']
            if named and len(named) < len(live):
                live, decided_by = named, 'product names the substrate'
            if len(live) > 1:
                cons = conservation.get(fid, {})
                top = max(live, key=lambda c: cons.get(c['cls'], 0.0))
                tied = [c for c in live
                        if abs(cons.get(c['cls'], 0.0) - cons.get(top['cls'], 0.0)) < 1e-9]
                if len(tied) == 1:
                    live, decided_by = tied, 'conservation'
                else:
                    near = [c for c in tied if c['kb'] is not None]
                    if near:
                        best = min(near, key=lambda c: c['kb'])
                        close = [c for c in near if c['kb'] - best['kb'] < 0.5]
                        if len(close) == 1:
                            live, decided_by = close, 'distance to pepM'
                        else:
                            live = close
                    if decided_by is None and references and len(live) > 1:
                        scored = [(support(c['seq'], references).get(c['cls'], {})
                                   .get('pct_id', 0.0), c) for c in live if c['seq']]
                        if len(scored) > 1:
                            scored.sort(key=lambda x: -x[0])
                            m = scored[0][0] - scored[1][0]
                            if m >= AMBIGUOUS_MARGIN and scored[0][0] > BACKGROUND_CEILING_PCT:
                                live, decided_by = [scored[0][1]], 'reference identity'
            if decided_by is None:
                decided_by = 'AMBIGUOUS (marker precedence)'
            st['cls'], st['marker'] = live[0]['cls'], live[0]['marker']
            st['chosen'] = live[0]
        else:
            st['chosen'] = cands[0] if cands else None
        st['decided_by'] = decided_by or 'single candidate'
        classifications[label] = st['cls']

        chosen = st['chosen']
        if references and chosen and chosen.get('seq'):
            seqs = st['seqs']
            pepm_seqs = seqs.get('PEP_mutase') or seqs.get('PEP_mutase_1') or []
            pepm_pid, pepm_ref, pepm_len = 0.0, '-', 0
            if pepm_seqs and pepm_references:
                pseq = max(pepm_seqs, key=len)
                pepm_len = len(pseq)
                best = max(((percent_identity(pseq, rs), rn) for rn, rs in pepm_references),
                           default=(0.0, '-'))
                pepm_pid, pepm_ref = round(best[0], 1), best[1]
            support_rows[label] = (
                st['cls'], st['marker'], len(chosen['seq']),
                support(chosen['seq'], references), pepm_pid, pepm_ref, pepm_len,
                (st['decided_by'], sorted({c['cls'] for c in cands}),
                 chosen.get('kb'), conservation.get(fid, {}).get(st['cls'])))

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
