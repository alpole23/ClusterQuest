#!/usr/bin/env python3
"""
Build phylogenetic trees for phosphonate BGC coupling enzyme analysis.

Two tree types:
  A: Combined pepM tree — universal marker (SMCOG1231 / PEP_mutase), all BGCs + reference
     anchors. Annotated with coupling enzyme class, GCF family, and source (ref vs query).
     Unknown BGCs fall naturally into the closest clade.

  B: Per-class coupling enzyme trees — one tree per class using the class-defining marker
     gene. Ppd and Ppd-CDP share one tree (same enzyme; class distinction is an annotation
     layer). References anchor each class tree.

HMM strategy (no external MSA tool required):
  1. hmmbuild from single seed reference    → initial HMM
  2. hmmalign all references to initial HMM → aligned references
  3. hmmbuild from aligned references       → refined HMM
  4. hmmalign all (refs + queries)          → final alignment
  5. FastTree (LG model)                   → .nwk tree

Usage:
    python scripts/bgc_coupling_tree.py \\
        --antismash_dir  results/antismash_results/Pantoea \\
        --metadata       results/bgc_trees/Pantoea/phosphonate_metadata.json \\
        --coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \\
        --ref_pepm_faa   assets/reference_sequences/reference_pepM.faa \\
        --ref_coupling_faa assets/reference_sequences/reference_coupling_enzymes.faa \\
        --outdir         results/bgc_trees/Pantoea/coupling_enzyme_trees \\
        --hmmbuild       /path/to/hmmbuild \\
        --hmmalign       /path/to/hmmalign \\
        --hmmsearch      /path/to/hmmsearch \\
        --fasttree       /path/to/FastTree \\
        --tree           both
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict

from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import (COUPLING_COLORS as _BASE_COUPLING_COLORS, LEGACY_CLASS_NAMES,
                              GCF_PALETTE, load_coupling_classes)

# Add the Reference class used only in coupling trees
COUPLING_COLORS = {**_BASE_COUPLING_COLORS, 'Reference': '#333333'}

# Coupling enzyme class for each reference pepM BGC
REF_PEPM_CLASS = {
    'BGC0000904':          'Synthase',      # FR-900098, FrbD pepM
    'BGC0000897':          'Decarboxylase', # Dehydrophos, DhpE pepM
    'BGC0000938':          'Decarboxylase', # Fosfomycin, Fom1 pepM
    'BGC0000806':          'Decarboxylase', # 2-AEP, Glycomyces pepM
    'Phosphonoalamide_BGC':'Transaminase',  # PnaD pepM
    'Valinophos_BGC':      'Reductase',     # VlpA pepM
    'Pantaphos_BGC':       'Synthase',      # Pantaphos, HvrA pepM
}

# Coupling enzyme class for each reference coupling enzyme entry (by protein name)
REF_COUPLING_CLASS = {
    'FrbC': 'Synthase',
    'HvrC': 'Synthase',
    'DhpF': 'Decarboxylase',
    'Fom2': 'Decarboxylase',
    'Ppd':  'Decarboxylase',
    'VlpB': 'Reductase',
    'PnaA': 'Transaminase',
}

# Markers for coupling enzyme CDS extraction.
# TPP_enzyme_C (Pfam clusterhmmer) is used for Decarboxylase and
# Decarboxylase-Nucleotidyltransferase rather than SMCOG1055: all Decarboxylase BGCs
# carry both annotations, but Decarboxylase-Nucleotidyltransferase BGCs have
# TPP_enzyme_C only — their ThDP decarboxylases are too divergent from the SMCOG1055
# seed to get a hit. Both classes share one combined tree keyed as 'Decarboxylase'.
CLASS_MARKERS = {
    'Synthase':                             ('smcog',  'SMCOG1271'),
    'Decarboxylase':                        ('domain', 'TPP_enzyme_C'),
    'Decarboxylase-Nucleotidyltransferase': ('domain', 'TPP_enzyme_C'),
    'Reductase':                            ('domain', 'Fe-ADH'),
    'Transaminase':                         ('smcog',  'SMCOG1013'),
}

# ─── Data loading ─────────────────────────────────────────────────────────────




def load_metadata(path):
    """Read BGC metadata JSON → {label: {organism, gcf, gbk_path}}."""
    with open(path) as f:
        records = json.load(f)
    meta = {}
    for r in records:
        lbl = r['label']
        gcf = next(
            (f['family_id'] for f in r.get('families', []) if f['cutoff'] == 0.3),
            None,
        )
        meta[lbl] = {
            'organism': r.get('organism', ''),
            'gcf':      gcf,
            'gbk_path': r.get('gbk_path', ''),
        }
    return meta


def build_json_index(antismash_dir):
    """Walk antismash_dir → {genome_name: json_path}."""
    index = {}
    for genome in os.listdir(antismash_dir):
        jp = os.path.join(antismash_dir, genome, f'{genome}.json')
        if os.path.exists(jp):
            index[genome] = jp
    return index


def genome_from_gbk_path(gbk_path):
    return os.path.basename(os.path.dirname(gbk_path))


def parse_label(label):
    """'CONTIG.regionNNN' → (contig_id, zero-padded region str)."""
    label = re.sub(r'_\d+$', '', label)
    m = re.search(r'^(.+?)\.region(\d+)$', label)
    if m:
        return m.group(1), m.group(2).zfill(3)
    return label, '001'


def parse_location_bounds(loc_str):
    """Extract (min_start, max_end) from an antiSMASH location string."""
    coords = re.findall(r'\[(\d+):(\d+)\]', str(loc_str))
    if not coords:
        return None, None
    starts = [int(s) for s, _ in coords]
    ends   = [int(e) for _, e in coords]
    return min(starts), max(ends)


# ─── Sequence extraction ──────────────────────────────────────────────────────

def _cds_in_region(feat, region_start, region_end):
    """Return True if a CDS feature overlaps the region."""
    if region_start is None:
        return True
    start, end = parse_location_bounds(feat.get('location', ''))
    if start is None:
        return True
    return not (end < region_start or start > region_end)


def _get_region_bounds(rec, contig_id, region_num):
    """Find the phosphonate region and return its (start, end)."""
    for feat in rec.get('features', []):
        if feat.get('type') != 'region':
            continue
        rnum = str(feat.get('qualifiers', {}).get('region_number', ['?'])[0]).zfill(3)
        products = feat.get('qualifiers', {}).get('product', [])
        if rnum == region_num and any('phosphonate' in p for p in products):
            return parse_location_bounds(feat.get('location', ''))
    return None, None


def extract_cds_from_json(json_path, contig_id, region_num,
                           is_pepm=False, smcog=None, domain=None):
    """
    Find and return (translation, gene_name) for a CDS within the specified
    phosphonate region matching the given marker.

    Specify one of:
      is_pepm=True       — match SMCOG1231 or PEP_mutase sec_met_domain
      smcog='SMCOG1271'  — match a specific SMCOG annotation
      domain='Fe-ADH'    — match a specific sec_met_domain or rule-based-cluster
    """
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception:
        return None, None

    for rec in data['records']:
        if contig_id not in rec.get('id', ''):
            continue

        r_start, r_end = _get_region_bounds(rec, contig_id, region_num)

        for feat in rec.get('features', []):
            if feat.get('type') != 'CDS':
                continue
            if not _cds_in_region(feat, r_start, r_end):
                continue

            quals = feat.get('qualifiers', {})
            translation = quals.get('translation', [''])[0]
            if not translation:
                continue

            gene_functions = quals.get('gene_functions', [])
            sec_met        = quals.get('sec_met_domain', [])

            matched = False
            if is_pepm:
                matched = (
                    any('SMCOG1231' in gf for gf in gene_functions) or
                    any('PEP_mutase' in sd for sd in sec_met)
                )
            elif smcog:
                matched = any(smcog in gf for gf in gene_functions)
            elif domain:
                matched = (
                    any(domain in sd for sd in sec_met) or
                    any('rule-based-clusters' in gf and domain in gf
                        for gf in gene_functions)
                )

            if matched:
                gene = quals.get('gene', [''])[0] or quals.get('locus_tag', [''])[0]
                return translation, gene

    return None, None


def extract_all_cds_from_region(json_path, contig_id, region_num):
    """Return all (translation, locus_tag) pairs for every CDS in the phosphonate region."""
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception:
        return []
    results = []
    for rec in data['records']:
        if contig_id not in rec.get('id', ''):
            continue
        r_start, r_end = _get_region_bounds(rec, contig_id, region_num)
        for feat in rec.get('features', []):
            if feat.get('type') != 'CDS':
                continue
            if not _cds_in_region(feat, r_start, r_end):
                continue
            quals = feat.get('qualifiers', {})
            translation = quals.get('translation', [''])[0]
            if not translation:
                continue
            gene = quals.get('gene', [''])[0] or quals.get('locus_tag', [''])[0]
            results.append((translation, gene))
    return results


# ─── HMM build and alignment ──────────────────────────────────────────────────

def write_fasta(records, path):
    with open(path, 'w') as f:
        SeqIO.write(records, f, 'fasta')


def run_hmmbuild(input_faa, out_hmm, hmmbuild_bin, name='profile'):
    cmd = [hmmbuild_bin, '--amino', '-n', name, out_hmm, input_faa]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  hmmbuild failed:\n{r.stderr[:600]}', file=sys.stderr)
        return False
    return True


def run_hmmalign(hmm, seqs_faa, out_afa, hmmalign_bin):
    cmd = [hmmalign_bin, '--amino', '--trim', '--outformat', 'afa', hmm, seqs_faa]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  hmmalign failed:\n{r.stderr[:600]}', file=sys.stderr)
        return False
    with open(out_afa, 'w') as f:
        f.write(r.stdout)
    return True


def run_hmmsearch(hmm_path, seqs_faa, out_tbl, hmmsearch_bin, evalue=1e-3):
    cmd = [hmmsearch_bin, '--tblout', out_tbl, '--noali',
           '-E', str(evalue), hmm_path, seqs_faa]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  hmmsearch failed:\n{r.stderr[:600]}', file=sys.stderr)
        return False
    return True


def parse_hmmsearch_tbl(tbl_path):
    """Parse hmmsearch --tblout → {seq_id: (bitscore, evalue)}."""
    hits = {}
    with open(tbl_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            seq_id = parts[0]
            try:
                evalue   = float(parts[4])
                bitscore = float(parts[5])
            except ValueError:
                continue
            if seq_id not in hits or evalue < hits[seq_id][1]:
                hits[seq_id] = (bitscore, evalue)
    return hits


def hmm_build_refs(ref_records, workdir, hmmbuild_bin, hmmalign_bin, label):
    """
    Build a refined HMM from reference sequences only.
    Returns path to final.hmm, or None on failure.
    """
    os.makedirs(workdir, exist_ok=True)
    seed_faa  = os.path.join(workdir, 'seed.faa')
    refs_faa  = os.path.join(workdir, 'refs.faa')
    hmm_init  = os.path.join(workdir, 'init.hmm')
    refs_afa  = os.path.join(workdir, 'refs_aligned.faa')
    hmm_final = os.path.join(workdir, 'final.hmm')

    write_fasta([ref_records[0]], seed_faa)
    write_fasta(ref_records, refs_faa)

    print(f'    [1/4] hmmbuild from seed: {ref_records[0].id}')
    if not run_hmmbuild(seed_faa, hmm_init, hmmbuild_bin, label + '_init'):
        return None

    print(f'    [2/4] hmmalign {len(ref_records)} references...')
    if not run_hmmalign(hmm_init, refs_faa, refs_afa, hmmalign_bin):
        return None

    print(f'    [3/4] hmmbuild refined HMM from aligned references...')
    if not run_hmmbuild(refs_afa, hmm_final, hmmbuild_bin, label):
        return None

    return hmm_final


def hmm_align_seqs(hmm_path, records, out_afa, hmmalign_bin):
    """Align records to an existing HMM. Returns MultipleSeqAlignment or None."""
    tmp_faa = out_afa + '.input.faa'
    write_fasta(records, tmp_faa)
    if not run_hmmalign(hmm_path, tmp_faa, out_afa, hmmalign_bin):
        return None
    alignment = AlignIO.read(out_afa, 'fasta')
    print(f'    Alignment: {len(alignment)} seqs × {alignment.get_alignment_length()} cols')
    return alignment


def hmm_align_all(ref_records, query_records, workdir, hmmbuild_bin, hmmalign_bin, label):
    """Build refined HMM from refs then align all sequences. Returns (alignment, out_afa) or (None, None)."""
    hmm_final = hmm_build_refs(ref_records, workdir, hmmbuild_bin, hmmalign_bin, label)
    if hmm_final is None:
        return None, None
    all_records = ref_records + query_records
    out_afa = os.path.join(workdir, f'{label}_all.afa')
    print(f'    [4/4] hmmalign {len(all_records)} sequences...')
    alignment = hmm_align_seqs(hmm_final, all_records, out_afa, hmmalign_bin)
    if alignment is None:
        return None, None
    return alignment, out_afa


# ─── Extraction methods ───────────────────────────────────────────────────────

def extract_via_annotation(bgc_labels, metadata, json_index, marker_type, marker_value):
    """Extract coupling enzyme CDS using antiSMASH annotation markers.
    Returns (found_records, missing_labels).
    """
    found, missing = [], []
    for lbl in bgc_labels:
        meta = metadata.get(lbl, {})
        contig_id, region_num = parse_label(lbl)
        genome    = genome_from_gbk_path(meta.get('gbk_path', '')) if meta.get('gbk_path') else None
        json_path = json_index.get(genome)
        seq, _ = extract_cds_from_json(
            json_path, contig_id, region_num,
            smcog=marker_value  if marker_type == 'smcog'  else None,
            domain=marker_value if marker_type == 'domain' else None,
        ) if json_path else (None, None)
        if seq:
            found.append(SeqRecord(Seq(seq), id=lbl, description=''))
        else:
            missing.append(lbl)
    return found, missing


def extract_via_hmmsearch(bgc_labels, metadata, json_index, hmm_path, workdir, hmmsearch_bin,
                           evalue=1e-3):
    """Extract coupling enzyme CDS by hmmsearch against the class HMM.
    All CDS in each region are searched; the highest-scoring hit per BGC is selected.
    Returns (found_records, missing_labels).
    """
    os.makedirs(workdir, exist_ok=True)
    all_cds_faa = os.path.join(workdir, 'all_cds.faa')
    tbl_path    = os.path.join(workdir, 'hmmsearch.tbl')

    cds_to_bgc  = {}   # cds_id → bgc_label
    cds_seq_map = {}   # cds_id → sequence string
    for lbl in bgc_labels:
        meta = metadata.get(lbl, {})
        contig_id, region_num = parse_label(lbl)
        genome    = genome_from_gbk_path(meta.get('gbk_path', '')) if meta.get('gbk_path') else None
        json_path = json_index.get(genome)
        if not json_path:
            continue
        for idx, (seq, _) in enumerate(extract_all_cds_from_region(json_path, contig_id, region_num)):
            cds_id = f'{lbl}__cds__{idx}'
            cds_to_bgc[cds_id]  = lbl
            cds_seq_map[cds_id] = seq

    if not cds_seq_map:
        return [], list(bgc_labels)

    write_fasta(
        [SeqRecord(Seq(s), id=i, description='') for i, s in cds_seq_map.items()],
        all_cds_faa,
    )
    if not run_hmmsearch(hmm_path, all_cds_faa, tbl_path, hmmsearch_bin, evalue):
        return [], list(bgc_labels)

    hits = parse_hmmsearch_tbl(tbl_path)

    best_hit = {}  # bgc_label → (bitscore, cds_id)
    for cds_id, (bitscore, _) in hits.items():
        lbl = cds_to_bgc.get(cds_id)
        if lbl is None:
            continue
        if lbl not in best_hit or bitscore > best_hit[lbl][0]:
            best_hit[lbl] = (bitscore, cds_id)

    found, missing = [], []
    for lbl in bgc_labels:
        if lbl in best_hit:
            cds_id = best_hit[lbl][1]
            found.append(SeqRecord(Seq(cds_seq_map[cds_id]), id=lbl, description=''))
        else:
            missing.append(lbl)
    return found, missing


def extract_with_hmm_fallback(bgc_labels, metadata, json_index, marker_type, marker_value,
                               hmm_path, workdir, hmmsearch_bin):
    """Annotation-first extraction with HMM rescue for annotation-missing BGCs.

    Uses antiSMASH annotation markers (SMCOG/domain) as the primary source — zero extra
    compute since antiSMASH already ran those HMMs.  For any BGCs that annotation missed
    (e.g. divergent sequences that fall below SMCOG thresholds), runs hmmsearch against
    the class reference HMM and selects the highest-scoring CDS per BGC.
    """
    ann_found, ann_missing = extract_via_annotation(
        bgc_labels, metadata, json_index, marker_type, marker_value)

    if not ann_missing or not hmmsearch_bin:
        return ann_found, ann_missing

    print(f'      HMM fallback for {len(ann_missing)} annotation-missing BGCs...')
    hmm_rescued, still_missing = extract_via_hmmsearch(
        ann_missing, metadata, json_index, hmm_path,
        os.path.join(workdir, 'fallback_search'), hmmsearch_bin)

    if hmm_rescued:
        print(f'      Rescued: {len(hmm_rescued)}  |  Still missing: {len(still_missing)}')
    return ann_found + hmm_rescued, still_missing


# ─── Tree building ─────────────────────────────────────────────────────────────

def build_fasttree(alignment, out_nwk, fasttree_bin):
    # hmmalign uses '.' for insert-state gaps and lowercase for insert-state residues.
    # FastTree requires uppercase residues and '-' only.
    cleaned = MultipleSeqAlignment([
        SeqRecord(Seq(str(r.seq).upper().replace('.', '-')), id=r.id, description='')
        for r in alignment
    ])
    tmp_faa = out_nwk + '.input.faa'
    write_fasta(list(cleaned), tmp_faa)
    print(f'    Running FastTree ({len(cleaned)} seqs × {cleaned.get_alignment_length()} cols)...')
    cmd = [fasttree_bin, '-quiet', '-lg', tmp_faa]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  FastTree failed:\n{r.stderr[:600]}', file=sys.stderr)
        return None
    with open(out_nwk, 'w') as f:
        f.write(r.stdout)
    print(f'    Written: {out_nwk}')
    return out_nwk


# ─── iTOL annotation writers ──────────────────────────────────────────────────

def write_itol_colorstrip(labels, color_fn, display_fn, dataset_label, legend_items, out_path):
    """
    Generic DATASET_COLORSTRIP writer.
    color_fn(label)   → hex color string
    display_fn(label) → display label string (shown in strip tooltip)
    legend_items      → [(label_str, color_str), ...]
    """
    with open(out_path, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{dataset_label}\n')
        f.write('COLOR\t#333333\n')
        f.write('STRIP_WIDTH\t40\n')
        f.write('SHOW_BORDER\t1\n')
        f.write('BORDER_WIDTH\t0.5\n')
        if legend_items:
            f.write(f'LEGEND_TITLE\t{dataset_label}\n')
            f.write('LEGEND_SHAPES\t' + '\t'.join('1' for _ in legend_items) + '\n')
            f.write('LEGEND_COLORS\t' + '\t'.join(c for _, c in legend_items) + '\n')
            f.write('LEGEND_LABELS\t' + '\t'.join(l for l, _ in legend_items) + '\n')
        f.write('DATA\n')
        for lbl in labels:
            f.write(f'{lbl}\t{color_fn(lbl)}\t{display_fn(lbl)}\n')


def write_itol_text(labels, text_fn, dataset_label, out_path):
    """DATASET_TEXT writer for leaf label annotations (e.g. organism names)."""
    with open(out_path, 'w') as f:
        f.write('DATASET_TEXT\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{dataset_label}\n')
        f.write('COLOR\t#333333\n')
        f.write('DATA\n')
        for lbl in labels:
            text = text_fn(lbl)
            if text:
                # node_id, text, position (1=after), color, style, size_factor
                f.write(f'{lbl}\t{text}\t1\t#333333\tnormal\t1\n')


def write_tree_itol(seq_labels, get_class_fn, metadata, ref_records, outdir):
    """Write iTOL annotation files for a coupling enzyme tree.

    get_class_fn(label) → class string.  Caller supplies the reference-labeling
    logic: Tree A maps REF| labels to their known coupling class via REF_PEPM_CLASS;
    Tree B marks all REF| labels as 'Reference'.
    ref_records: original reference SeqRecords without the REF| prefix.
    """
    present = sorted(set(get_class_fn(l) for l in seq_labels))
    coupling_legend = [(cls, COUPLING_COLORS.get(cls, '#aaaaaa')) for cls in present]

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: COUPLING_COLORS.get(get_class_fn(l), '#aaaaaa'),
        display_fn = get_class_fn,
        dataset_label = 'Coupling enzyme class',
        legend_items  = coupling_legend,
        out_path = os.path.join(outdir, 'itol_coupling_class.txt'),
    )

    gcf_set = sorted(set(
        metadata[l]['gcf'] for l in seq_labels
        if l in metadata and metadata[l]['gcf'] is not None
    ))
    gcf_color = {g: GCF_PALETTE[i % len(GCF_PALETTE)] for i, g in enumerate(gcf_set)}

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: gcf_color.get(metadata[l]['gcf'], '#dddddd')
                               if l in metadata and metadata[l]['gcf'] is not None
                               else '#333333' if l.startswith('REF|') else '#dddddd',
        display_fn = lambda l: f'GCF-{metadata[l]["gcf"]}'
                               if l in metadata and metadata[l]['gcf'] is not None
                               else 'Reference' if l.startswith('REF|') else 'No GCF',
        dataset_label = 'GCF family',
        legend_items  = [(f'GCF-{g}', gcf_color[g]) for g in gcf_set],
        out_path = os.path.join(outdir, 'itol_gcf.txt'),
    )

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: '#333333' if l.startswith('REF|') else '#cccccc',
        display_fn = lambda l: 'Reference' if l.startswith('REF|') else 'Query',
        dataset_label = 'Source',
        legend_items  = [('Reference', '#333333'), ('Query (Pantoea)', '#cccccc')],
        out_path = os.path.join(outdir, 'itol_source.txt'),
    )

    ref_map = {}
    for r in ref_records:
        parts = r.id.split('|')  # BGC|acc|name|function|organism
        org = parts[4].replace('_', ' ') if len(parts) > 4 else r.id
        ref_map[f'REF|{r.id}'] = org

    write_itol_text(
        seq_labels,
        text_fn = lambda l: ref_map.get(l) or metadata.get(l, {}).get('organism', ''),
        dataset_label = 'Organism',
        out_path = os.path.join(outdir, 'itol_organism.txt'),
    )


# ─── Tree A ────────────────────────────────────────────────────────────────────

def build_tree_a(args, metadata, coupling_classes, json_index, ref_pepm_records):
    print('\n=== Tree A: Combined pepM tree ===')
    outdir = os.path.join(args.outdir, 'tree_A')
    os.makedirs(outdir, exist_ok=True)

    print(f'Extracting pepM sequences from {len(metadata)} BGCs...')
    query_records = []
    missing = []

    for lbl, meta in metadata.items():
        contig_id, region_num = parse_label(lbl)
        genome    = genome_from_gbk_path(meta['gbk_path']) if meta['gbk_path'] else None
        json_path = json_index.get(genome)

        seq, _ = extract_cds_from_json(json_path, contig_id, region_num, is_pepm=True) \
                 if json_path else (None, None)

        if seq:
            query_records.append(SeqRecord(Seq(seq), id=lbl, description=''))
        else:
            missing.append(lbl)

    print(f'  Extracted: {len(query_records)}  |  Missing pepM: {len(missing)}')
    if missing:
        missing_log = os.path.join(outdir, 'missing_pepm.txt')
        with open(missing_log, 'w') as f:
            f.write('\n'.join(missing))
        print(f'  Missing labels written to: {missing_log}')

    ref_records = [
        SeqRecord(Seq(str(r.seq)), id=f'REF|{r.id}', description=r.description)
        for r in ref_pepm_records
    ]

    alignment, _ = hmm_align_all(
        ref_records, query_records,
        os.path.join(outdir, 'hmm_work'),
        args.hmmbuild, args.hmmalign, 'pepm',
    )
    if alignment is None:
        print('Tree A: alignment failed.', file=sys.stderr)
        return None

    seq_labels = [r.id for r in alignment]
    nwk_path = os.path.join(outdir, 'pepm_tree.nwk')
    build_fasttree(alignment, nwk_path, args.fasttree)

    print('  Writing iTOL annotations...')
    get_class = lambda l: REF_PEPM_CLASS.get(l.split('|')[1], 'Reference') \
                          if l.startswith('REF|') else coupling_classes.get(l, 'Unknown')
    write_tree_itol(
        seq_labels,
        get_class_fn = get_class,
        metadata     = metadata,
        ref_records  = ref_pepm_records,
        outdir       = outdir,
    )
    print(f'Tree A complete → {outdir}')

    newick = open(nwk_path).read() if os.path.exists(nwk_path) else None
    leaf_meta = {
        lbl: {
            'class': get_class(lbl),
            'gcf':   metadata.get(lbl, {}).get('gcf'),
            'is_ref': lbl.startswith('REF|'),
        }
        for lbl in seq_labels
    }
    return {'n_bgcs': len(query_records), 'newick': newick, 'leaf_metadata': leaf_meta}


# ─── Tree B ────────────────────────────────────────────────────────────────────

def build_tree_b(args, metadata, coupling_classes, json_index, ref_coupling_records):
    print('\n=== Tree B: Per-class coupling enzyme trees ===')
    outdir_b = os.path.join(args.outdir, 'tree_B')
    os.makedirs(outdir_b, exist_ok=True)

    # Group BGCs by class; merge Ppd + Ppd-CDP into one tree.
    # Skip BiG-SCAPE sub-record duplicates (_1, _2, _3 suffixes) — they share the
    # same CDS with the main region label and would produce duplicate sequences.
    class_bgcs = defaultdict(list)
    for lbl, cls in coupling_classes.items():
        if re.search(r'_\d+$', lbl):
            continue
        key = 'Decarboxylase' if cls in ('Decarboxylase', 'Decarboxylase-Nucleotidyltransferase') else cls
        class_bgcs[key].append(lbl)

    # Group reference sequences by class; merge both Decarboxylase classes into one tree
    ref_by_class = defaultdict(list)
    for r in ref_coupling_records:
        parts = r.id.split('|')
        name = parts[2] if len(parts) > 2 else ''
        cls  = REF_COUPLING_CLASS.get(name)
        if cls:
            key = 'Decarboxylase' if cls in ('Decarboxylase', 'Decarboxylase-Nucleotidyltransferase') else cls
            ref_by_class[key].append(r)

    tree_b_results = []

    for class_key, bgc_labels in sorted(class_bgcs.items()):
        if class_key == 'Unknown':
            print(f'\n  Skipping Unknown class (no coupling enzyme to extract)')
            continue

        print(f'\n  Class: {class_key} ({len(bgc_labels)} BGCs)')
        class_outdir = os.path.join(outdir_b, class_key.replace('+', '_'))
        os.makedirs(class_outdir, exist_ok=True)

        refs = ref_by_class.get(class_key, [])
        if not refs:
            print(f'    No references for {class_key} — skipping.')
            continue

        ref_records = [
            SeqRecord(Seq(str(r.seq)), id=f'REF|{r.id}', description=r.description)
            for r in refs
        ]

        # Build HMM from references (steps 1–3; used by all extraction methods)
        hmm_workdir = os.path.join(class_outdir, 'hmm_work')
        label       = class_key.lower().replace('+', '_')
        hmm_final   = hmm_build_refs(ref_records, hmm_workdir, args.hmmbuild, args.hmmalign, label)
        if hmm_final is None:
            print(f'    Tree B/{class_key}: HMM build failed.', file=sys.stderr)
            continue

        # Extract coupling enzyme CDS: annotation-first with HMM fallback for misses
        marker_type, marker_value = CLASS_MARKERS[class_key]
        query_records, missing = extract_with_hmm_fallback(
            bgc_labels, metadata, json_index, marker_type, marker_value,
            hmm_final, hmm_workdir, getattr(args, 'hmmsearch', None))
        print(f'    Extracted: {len(query_records)}  |  Missing: {len(missing)}')

        if not query_records:
            print(f'    No sequences extracted for {class_key} — skipping tree.')
            continue

        # Align refs + query sequences against the class HMM (step 4)
        all_records = ref_records + query_records
        out_afa     = os.path.join(hmm_workdir, f'{label}_all.afa')
        print(f'    [4/4] hmmalign {len(all_records)} sequences...')
        alignment = hmm_align_seqs(hmm_final, all_records, out_afa, args.hmmalign)
        if alignment is None:
            print(f'    Tree B/{class_key}: alignment failed.', file=sys.stderr)
            continue

        seq_labels = [r.id for r in alignment]
        out_nwk = os.path.join(class_outdir, f'{class_key.replace("+", "_")}_tree.nwk')
        build_fasttree(alignment, out_nwk, args.fasttree)

        print(f'    Writing iTOL annotations...')
        get_class = lambda l: 'Reference' if l.startswith('REF|') \
                              else coupling_classes.get(l, 'Unknown')
        write_tree_itol(
            seq_labels,
            get_class_fn = get_class,
            metadata     = metadata,
            ref_records  = refs,
            outdir       = class_outdir,
        )
        print(f'    Tree B/{class_key} complete → {class_outdir}')

        newick = open(out_nwk).read() if os.path.exists(out_nwk) else None
        leaf_meta = {
            lbl: {
                'class': get_class(lbl),
                'gcf':   metadata.get(lbl, {}).get('gcf'),
                'is_ref': lbl.startswith('REF|'),
            }
            for lbl in seq_labels
        }
        tree_b_results.append({
            'class_key':     class_key,
            'n_bgcs':        len(query_records),
            'newick':        newick,
            'leaf_metadata': leaf_meta,
        })

    return tree_b_results


# ─── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Build pepM and coupling enzyme phylogenetic trees for phosphonate BGCs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--antismash_dir',       required=True,
                        help='antiSMASH results directory')
    parser.add_argument('--metadata',            required=True,
                        help='BGC metadata JSON from bgc_pfam_tree.py')
    parser.add_argument('--coupling_annotation', required=True,
                        help='iTOL colorstrip from bgc_coupling_annotation.py')
    parser.add_argument('--ref_pepm_faa',        required=True,
                        help='Reference pepM sequences (FASTA)')
    parser.add_argument('--ref_coupling_faa',    required=True,
                        help='Reference coupling enzyme sequences (FASTA)')
    parser.add_argument('--outdir',              required=True,
                        help='Output directory')
    parser.add_argument('--hmmbuild',            required=True,
                        help='Path to hmmbuild binary')
    parser.add_argument('--hmmalign',            required=True,
                        help='Path to hmmalign binary')
    parser.add_argument('--hmmsearch',           default=None,
                        help='Path to hmmsearch binary (enables HMM fallback for annotation-missing BGCs)')
    parser.add_argument('--fasttree',            required=True,
                        help='Path to FastTree binary')
    parser.add_argument('--tree', default='both', choices=['A', 'B', 'both'],
                        help='Which tree(s) to build (default: both)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print('Loading metadata and annotations...')
    metadata        = load_metadata(args.metadata)
    coupling_classes = load_coupling_classes(args.coupling_annotation)
    print(f'  {len(metadata)} BGCs in metadata')
    print(f'  {len(coupling_classes)} coupling class assignments')

    print('Building antiSMASH JSON index...')
    json_index = build_json_index(args.antismash_dir)
    print(f'  {len(json_index)} genome JSON files')

    ref_pepm_records     = list(SeqIO.parse(args.ref_pepm_faa,     'fasta'))
    ref_coupling_records = list(SeqIO.parse(args.ref_coupling_faa, 'fasta'))
    print(f'  {len(ref_pepm_records)} pepM reference sequences')
    print(f'  {len(ref_coupling_records)} coupling enzyme reference sequences')

    manifest = {}

    if args.tree in ('A', 'both'):
        tree_a = build_tree_a(args, metadata, coupling_classes, json_index, ref_pepm_records)
        if tree_a:
            manifest['tree_a'] = tree_a

    if args.tree in ('B', 'both'):
        tree_b = build_tree_b(args, metadata, coupling_classes, json_index, ref_coupling_records)
        if tree_b:
            manifest['tree_b'] = tree_b

    if manifest:
        manifest_path = os.path.join(args.outdir, 'coupling_trees_manifest.json')
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f)
        print(f'\nManifest written: {manifest_path}')

    print('\nDone.')


if __name__ == '__main__':
    main()
