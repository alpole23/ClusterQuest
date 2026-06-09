#!/usr/bin/env python3
"""Shared constants for BGC analysis pipeline."""

# BGC type colors — comprehensive palette covering all antiSMASH product types
BGC_COLORS = {
    # PKS types - Blue family
    'T1PKS': '#1f77b4',
    'T2PKS': '#4a90d9',
    'T3PKS': '#7eb3ed',
    'transAT-PKS': '#2c5aa0',
    'PKS-like': '#5dade2',
    'hglE-KS': '#85c1e9',
    # NRPS types - Red/Orange family
    'NRPS': '#d62728',
    'NRPS-like': '#ff6b6b',
    'thioamide-NRP': '#e74c3c',
    'NAPAA': '#c0392b',
    # Hybrid types - Purple family
    'PKS-NRPS_Hybrids': '#9b59b6',
    'NRPS-PKS_Hybrids': '#8e44ad',
    # RiPPs - Green family
    'RiPP': '#27ae60',
    'lanthipeptide': '#2ecc71',
    'lanthipeptide-class-i': '#27ae60',
    'lanthipeptide-class-ii': '#229954',
    'lanthipeptide-class-iii': '#1e8449',
    'lanthipeptide-class-iv': '#196f3d',
    'lanthipeptide-class-v': '#145a32',
    'thiopeptide': '#58d68d',
    'LAP': '#82e0aa',
    'lassopeptide': '#abebc6',
    'sactipeptide': '#d5f5e3',
    'bottromycin': '#a9dfbf',
    'cyanobactin': '#73c6b6',
    'microviridin': '#45b39d',
    'proteusin': '#16a085',
    'RRE-containing': '#138d75',
    'fungal-RiPP': '#117a65',
    'ranthipeptide': '#0e6655',
    'redox-cofactor': '#0b5345',
    'RiPP-like': '#7dcea0',
    'thioamitides': '#52be80',
    'epipeptide': '#48c9b0',
    'guanidinotides': '#1abc9c',
    'glycocin': '#17a589',
    'triceptide': '#148f77',
    'spliceotide': '#117864',
    'methanobactin': '#0e6251',
    'cyclic-lactone-autoinducer': '#85929e',
    'darobactin': '#76d7c4',
    'rcdpeptide': '#45b39d',
    # Terpenes - Yellow/Gold family
    'terpene': '#f39c12',
    # Saccharides - Brown family
    'saccharide': '#a0522d',
    'oligosaccharide': '#cd853f',
    'polysaccharide': '#8b4513',
    'amglyccycl': '#d2691e',
    # Alkaloids/Other nitrogen - Teal family
    'alkaloid': '#008080',
    'indole': '#20b2aa',
    'NI-siderophore': '#40e0d0',
    # Siderophores - Cyan family
    'siderophore': '#00ced1',
    'NAGGN': '#00bfff',
    # Fatty acids / Lipids - Olive family
    'fatty_acid': '#808000',
    'PUFA': '#9acd32',
    'ladderane': '#6b8e23',
    'hserlactone': '#556b2f',
    'acyl_amino_acids': '#8fbc8f',
    'N-acyl amino acid': '#8fbc8f',
    # Phosphonates - Pink family
    'phosphonate': '#ff69b4',
    'phosphonate-like': '#ffb6c1',
    # Aromatic compounds - Magenta family
    'arylpolyene': '#ff00ff',
    'resorcinol': '#da70d6',
    'stilbene': '#ee82ee',
    'phenazine': '#dda0dd',
    'aminocoumarin': '#ba55d3',
    # Nucleosides - Gray family
    'nucleoside': '#708090',
    # Bacteriocins - Dark cyan
    'bacteriocin': '#008b8b',
    'RaS-RiPP': '#5f9ea0',
    # Beta-lactams - Coral
    'betalactam': '#ff7f50',
    'beta-lactam': '#ff7f50',
    # Ectoine - Light purple
    'ectoine': '#dda0dd',
    'ectoine-like': '#e6e6fa',
    # Melanin - Dark gray
    'melanin': '#2f4f4f',
    # Butyrolactone - Peach
    'butyrolactone': '#ffdab9',
    # Blactam - Salmon
    'blactam': '#fa8072',
    # CDPS - Lavender
    'CDPS': '#e6e6fa',
    # Furan - Wheat
    'furan': '#f5deb3',
    # Prodigiosin - Crimson
    'prodigiosin': '#dc143c',
    # Cyanide - Light steel blue
    'cyanide': '#b0c4de',
    'hydrogen-cyanide': '#b0c4de',
    # Linaridin - Medium purple
    'linaridin': '#9370db',
    'linear-azol(in)e-containing-peptide': '#9370db',
    # Opine - Thistle
    'opine-like-metallophore': '#d8bfd8',
    # Other/Unknown - Gray
    'other': '#95a5a6',
    'Other': '#95a5a6',
    'unknown': '#bdc3c7',
    'Unknown': '#bdc3c7',
    'NA': '#ecf0f1',
    # Legacy names kept for backwards compatibility
    'hybrid': '#9b59b6',
    'NRP-metallophore': '#c5b0d5',
    'betalactone': '#7b4173',
}

# Gene function colors (matches antiSMASH conventions)
GENE_COLORS = {
    'biosynthetic': '#e74c3c',           # Red - core biosynthetic
    'biosynthetic-additional': '#e67e22', # Orange - additional biosynthetic
    'regulatory': '#27ae60',              # Green - regulatory
    'transport': '#3498db',               # Blue - transport
    'resistance': '#9b59b6',              # Purple - resistance
    'other': '#95a5a6',                   # Gray - other/unknown
}

# GCF color palette for coupling enzyme tree iTOL annotations (10-color cycling palette)
GCF_PALETTE = [
    '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
    '#e6ab02', '#a6761d', '#666666', '#a6cee3', '#b2df8a',
]

# GCF color palette for tree visualization
GCF_COLORS = [
    '#e41a1c',  # Red
    '#377eb8',  # Blue
    '#4daf4a',  # Green
    '#984ea3',  # Purple
    '#ff7f00',  # Orange
    '#ffff33',  # Yellow
    '#a65628',  # Brown
    '#f781bf',  # Pink
    '#999999',  # Gray
    '#66c2a5',  # Teal
    '#fc8d62',  # Coral
    '#8da0cb',  # Light blue
    '#e78ac3',  # Magenta
    '#a6d854',  # Lime
    '#ffd92f',  # Gold
    '#e5c494',  # Tan
]

# KCB similarity thresholds
KCB_THRESHOLDS = {
    'high': 75,
    'medium': 50,
    'low': 15,
}

# KCB similarity colors
KCB_COLORS = {
    'high': '#28a745',    # Green
    'medium': '#fd7e14',  # Orange
    'low': '#6c757d',     # Gray
}

# Coupling enzyme class colors (used across bgc_gcf_tree, bgc_gcf_heatmap,
# bgc_all_bgcs_tree, bgc_coupling_tree, bgc_coupling_annotation)
COUPLING_COLORS = {
    'Synthase':                        '#e41a1c',
    'Decarboxylase':                   '#377eb8',
    'Decarboxylase-Nucleotidyltransferase': '#984ea3',
    'Reductase':                       '#4daf4a',
    'Transaminase':                    '#ff7f00',
    'Unknown':                         '#aaaaaa',
}

# Display order for coupling enzyme classes in legends
COUPLING_ORDER = ['Synthase', 'Reductase', 'Decarboxylase', 'Decarboxylase-Nucleotidyltransferase', 'Transaminase', 'Unknown']

# Normalize legacy coupling class names from older annotation files
LEGACY_CLASS_NAMES = {
    'Fe-ADH':    'Reductase',
    'TPP+NTP':   'Decarboxylase-Nucleotidyltransferase',
    'PalB':      'Transaminase',
    'FrbC':      'Synthase',
    # pre-rename class IDs
    'FrbC-like': 'Synthase',
    'VlpB-like': 'Reductase',
    'PalB-like': 'Transaminase',
    'Ppd':       'Decarboxylase',
    'Ppd-CDP':   'Decarboxylase-Nucleotidyltransferase',
}

# Pfam accession → short human-readable name.
# Verified against antiSMASH clusterhmmer output on phosphonate BGCs.
# Add entries here as needed; unknown domains fall back to their accession.
DOMAIN_NAMES = {
    'PF13714': 'PEP_mutase',
    'PF00296': 'HEPD',           # 2-hydroxyethylphosphonate dioxygenase
    'PF02775': 'ThDP_C',         # phosphonopyruvate decarboxylase (C-term)
    'PF02776': 'ThDP_N',         # phosphonopyruvate decarboxylase (N-term)
    'PF00266': 'Aminotrans_V',   # 2-AEP transaminase (class V)
    'PF00155': 'Aminotrans_I',   # aminotransferase class I/II
    'PF00682': 'FrbC_HMGL',      # FrbC-like / phosphonomethylmalate synthase (HMGL superfamily)
    'PF13649': 'Radical_SAM',
    'PF08241': 'Methyltransf',
    'PF00330': 'Aconitase',
    'PF00694': 'Aconitase_C',
    'PF00149': 'Metallophos',    # calcineurin-like phosphoesterase
    'PF00571': 'CBS',
    'PF07690': 'MFS_transporter',
    'PF00005': 'ABC_ATPase',
    'PF01118': 'SemiAldDH_NAD',  # semialdehyde dehydrogenase
    'PF02774': 'SemiAldDH_C',
    'PF00464': 'SHMT',           # serine hydroxymethyltransferase
    'PF12804': 'NTP_transf',
    'PF13673': 'Acetyltransf',
    'PF07228': 'DUF1453',
    'PF00733': 'Asn_synthase',
    'PF13537': 'GATase',
    'PF00581': 'Rhodanese',
    'PF01613': 'Trp_syntA',
    'PF00202': 'Aminotrans_III',
    'PF00892': 'EamA_transporter',
    'PF22617': 'PF22617',
    'PF05321': 'PF05321',
}


def load_coupling_classes(path, region_only=False):
    """Parse iTOL coupling colorstrip → {gbk_basename: class_id}.

    Args:
        path:        Path to phosphonate_itol_coupling.txt
        region_only: Skip sub-record labels (those ending in _0.._9).
                     Use True when building trees where only region-level
                     BGC records are leaves.
    """
    classes = {}
    in_data = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == 'DATA':
                in_data = True
                continue
            if in_data and line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    label = parts[0]
                    cls = LEGACY_CLASS_NAMES.get(parts[2], parts[2])
                    if region_only and any(label.endswith(f'_{i}') for i in range(10)):
                        continue
                    classes[label] = cls
    return classes
