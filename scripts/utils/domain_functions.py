#!/usr/bin/env python3
"""Pfam accession -> biosynthetic role, for describing what a BGC actually does.

Why accessions and not product text. A keyword metric over NCBI product names was
tried first and gave an answer that inverted on the two lab-confirmed clusters. The
cause was not subtle: "serine hydroxymethyltransferase" matched `methyltransferase`
and "aspartate-semialdehyde dehydrogenase" matched `dehydrogenase`, so two core
amino-acid metabolism genes were counted as tailoring chemistry. Pfam separates them
by construction -- SHMT is PF00464, a methyltransferase is PF13649 -- and antiSMASH
scans every region with clusterhmmer regardless of whether NCBI annotated the
assembly, so coverage is uniform (77.6% of CDS on Erwiniaceae) rather than tracking
annotation quality (40.2%).

The categories answer "what role does this play in the cluster", not "what fold is
this". PRIMARY_METABOLISM exists because antiSMASH region boundaries sweep in
chromosomal neighbours: peptidyl-tRNA hydrolase and DnaB are not tailoring enzymes no matter how
many of them sit next to a BGC, and counting them is what made the first attempt
fail. MOBILE is separate because a transposase says something about how the cluster
got there, not about its chemistry.

Note what CORE means here: core to *phosphonate chemistry in general* -- the pepM
hallmark and the coupling enzymes that every phosphonate pathway branches through.
It does not mean "core to this BGC". A cluster's own defining biosynthetic genes
mostly live in TAILORING, which is where the chemistry that distinguishes one
phosphonate product from another is counted.

This map is CURATED AND INCOMPLETE. It covers the ~100 domains that account for 93%
of observed hits on Erwiniaceae; everything else returns 'other'. 'other' means "not
classified here", never "not a biosynthetic gene" -- do not read an absence as
evidence. Each entry carries its Pfam name so the assignment can be argued with.
"""

# ─── Core phosphonate pathway ────────────────────────────────────────────────
# The hallmark and the coupling enzymes that set which pathway runs downstream.
CORE = {
    'PF13714': 'PEP_mutase',          # phosphoenolpyruvate phosphomutase (pepM)
    'PF02775': 'TPP_enzyme_C',        # Ppd, phosphonopyruvate decarboxylase
    'PF02776': 'TPP_enzyme_N',
    'PF00682': 'HMGL-like',           # synthase route (FrbC/HvrC)
    # Always a second domain on the SAME protein as PF00682, never alone: 236 of 236
    # HCS_D2 hits co-occur with HMGL-like. It adds no independent signal and is here
    # only so a divergent synthase that loses the HMGL hit is still called core.
    'PF22617': 'HCS_D2',
    'PF00465': 'Fe-ADH',              # reductase route (VlpB)
    'PF25137': 'ADH_Fe_C',
    'PF00155': 'Aminotran_1_2',       # transaminase route (PalB/PnaA)
}

# ─── Tailoring: chemistry past the coupling step ─────────────────────────────
# What turns a simple phosphonate into an elaborated molecule.
TAILORING = {
    # methylation
    'PF13649': 'Methyltransf_25', 'PF08241': 'Methyltransf_11',
    'PF08242': 'Methyltransf_12', 'PF13847': 'Methyltransf_31',
    'PF00891': 'Methyltransf_2',  'PF01209': 'Ubie_methyltran',
    'PF13489': 'Methyltransf_23',
    # oxidation
    'PF00296': 'Bac_luciferase',      # LLM-class flavin-dependent monooxygenase
    'PF00067': 'p450', 'PF01494': 'FAD_binding_3', 'PF00743': 'FMO-like',
    'PF01613': 'Flavin_Reduct', 'PF00724': 'Oxidored_FMN',
    # redox
    'PF00107': 'ADH_zinc_N', 'PF08240': 'ADH_N', 'PF13561': 'adh_short_C2',
    'PF02826': '2-Hacid_dh_C', 'PF13454': 'NAD_binding_9',
    'PF01408': 'GFO_IDH_MocA', 'PF22725': 'GFO_IDH_MocA_C3',
    # isomerisation / dehydration — the pantaphos route runs through these
    'PF00330': 'Aconitase', 'PF00694': 'Aconitase_C',
    'PF01370': 'Epimerase', 'PF01177': 'Asp_Glu_race',
    # acyl / amide / ligation
    'PF13673': 'Acetyltransf_10', 'PF13508': 'Acetyltransf_7',
    'PF00583': 'Acetyltransf_1',
    'PF00733': 'Asn_synthase', 'PF13522': 'GATase_6',
    'PF13535': 'ATP-grasp_4', 'PF02655': 'ATP-grasp_3',
    'PF00551': 'Formyl_trans_N', 'PF02911': 'Formyl_trans_C',
    # nucleotide activation
    'PF00483': 'NTP_transferase', 'PF12804': 'NTP_transf_3',
    'PF01909': 'NTP_transf_2',
    # glycosylation
    'PF00534': 'Glycos_transf_1', 'PF00535': 'Glycos_transf_2',
    'PF13439': 'Glyco_transf_4',
    # halogenation
    'PF04820': 'Trp_halogenase',
    # carrier protein — marks assembly-line chemistry
    'PF00550': 'PP-binding',
    # other PLP transaminases acting downstream rather than as the coupling enzyme
    'PF00202': 'Aminotran_3', 'PF00266': 'Aminotran_5',
    'PF01041': 'DegT_DnrJ_EryC1', 'PF01053': 'Cys_Met_Meta_PP',
    'PF03756': 'AfsA',                # hotdog-fold, gamma-butyrolactone-type
    'PF00293': 'NUDIX',
}

# ─── Lipid handling ──────────────────────────────────────────────────────────
# Kept as its own category rather than folded into tailoring, because whether a
# phosphonate BGC makes a lipid or a small molecule is an open question here and
# these are the domains that would bear on it. NOTE the confirmed phosphonolipid
# cluster (P. ananatis LMG 5342 region 2) carries NONE of these: the headgroup is
# made by the BGC and conjugated by the cell's general lipid machinery. Presence is
# informative; absence is not evidence against a lipid product.
# NOTE only PF01066 is exercised by the Erwiniaceae data (42 hits). The other three
# have ZERO hits there — they are forward-looking entries for datasets that carry
# them, and nothing about this category has been tested beyond the one domain.
LIPID = {
    'PF01066': 'CDP-OH_P_transf',     # CDP-alcohol phosphatidyltransferase; 42 hits
    'PF00975': 'Thioesterase',        # unexercised on Erwiniaceae
    'PF01553': 'Acyltransferase',     # unexercised on Erwiniaceae
    'PF03279': 'Lipid_A_acyltrans',   # unexercised on Erwiniaceae
}

# ─── Transport ───────────────────────────────────────────────────────────────
TRANSPORT = {
    'PF07690': 'MFS_1', 'PF05232': 'BTP', 'PF00892': 'EamA', 'PF01810': 'LysE',
    'PF00005': 'ABC_tran', 'PF00664': 'ABC_membrane', 'PF01032': 'FecCD',
    'PF00593': 'TonB_dep_Rec', 'PF07715': 'Plug', 'PF00497': 'SBP_bac_3',
    'PF13407': 'Peripla_BP_4', 'PF01769': 'MgtE',
}

# ─── Regulation ──────────────────────────────────────────────────────────────
REGULATION = {
    'PF00356': 'LacI', 'PF00440': 'TetR_N', 'PF00126': 'HTH_1',
    'PF12833': 'HTH_18', 'PF13744': 'HTH_37', 'PF13276': 'HTH_21',
    'PF18607': 'HTH_54', 'PF01258': 'zf-dskA_traR',
}

# ─── Mobility ────────────────────────────────────────────────────────────────
# How the cluster arrived, not what it makes. Counted, never as chemistry.
MOBILE = {
    'PF13683': 'rve_3', 'PF00665': 'rve', 'PF01527': 'HTH_Tnp_1',
    'PF00589': 'Phage_integrase', 'PF22022': 'Phage_int_M',
    'PF13356': 'Arm-DNA-bind_3', 'PF04754': 'Transposase_31',
    'PF00717': 'Peptidase_S24', 'PF08775': 'ParB',
}

# ─── Primary metabolism ──────────────────────────────────────────────────────
# Named for what the ENZYME does (primary/central metabolism), not for its role in
# the cluster — "primary" alone reads as "primary to this BGC", which is the exact
# opposite of the meaning. Chromosomal neighbours antiSMASH's region boundary sweeps
# in. Listed explicitly so they are EXCLUDED from tailoring rather than silently
# inflating it — the exact failure that broke the product-keyword version.
#
# Corroborated by GENE-ORDER PROXIMITY to pepM — deliberately not "operon position".
# BGCs are not reliably single operons: measured here, the Ppd coupling enzyme
# (TPP_enzyme_C) is on the OPPOSITE strand from pepM in 48% of the clusters that carry
# it, so anything requiring co-orientation would discard half of a known pathway gene's
# calls. The measure is strand-agnostic for that reason.
#
# Rhodanese sits a median 14 genes / 12.3 kb from pepM and HAD-like hydrolase 12 genes
# / 10.5 kb, neither ever within 3 genes — what a flanking chromosomal gene looks like.
#
# THE TEST IS ASYMMETRIC AND WEAK IN ONE DIRECTION. Two of the four known coupling
# enzymes fail it: Fe-ADH is never within 3 genes of pepM (0%) and Aminotran_1_2 sits a
# median 12 genes away — the phosphonoalamide cluster famously places PalB far from
# pepM. Distance is therefore not evidence AGAINST pathway membership; it only fails to
# supply evidence for it. These assignments rest on the enzymes' known primary-metabolic
# function, with position as corroboration where it happens to agree.
PRIMARY_METABOLISM = {
    'PF00464': 'SHMT',                # serine hydroxymethyltransferase (one-carbon)
    # OPEN QUESTION. Semialdhyde_dh is the one entry here whose position argues against
    # the assignment: within 3 genes of pepM in 55% of cases, SAME strand 100% of the
    # time, median 989 bp intergenic — indistinguishable from the HMGL-like coupling
    # enzyme (same strand, 1,361 bp). It is a classic Lys/Thr biosynthesis enzyme, so
    # the adjacency may be genome organisation rather than pathway membership, but it
    # sits beside the pathway in GCF-18, the confirmed phosphonolipid. Moving it to
    # TAILORING would raise that cluster's elaboration from 2.00.
    'PF01118': 'Semialdhyde_dh',      # aspartate-semialdehyde dehydrogenase (Lys/Thr)
    'PF22698': 'Semialdhyde_dhC_1', 'PF02774': 'Semialdhyde_dhC',
    'PF00162': 'PGK',                 # phosphoglycerate kinase (glycolysis)
    'PF00180': 'Iso_dh',              # isocitrate dehydrogenase (TCA)
    'PF01195': 'Pept_tRNA_hydro',     # essential translation factor
    'PF01926': 'MMR_HSR1', 'PF06071': 'YchF-GTPase_C',
    'PF00772': 'DnaB', 'PF03796': 'DnaB_C', 'PF00817': 'IMS',
    'PF00849': 'PseudoU_synth_2', 'PF00085': 'Thioredoxin', 'PF01323': 'DSBA',
        # SpoIIE and CBS (one protein — a CBS-domain ser/thr phosphatase) sit within 3
    # genes of pepM in 86% and 75% of cases, which looked like pathway membership until
    # strand and spacing were checked: they are on the OPPOSITE strand 100% of the time
    # with a median 3,722 bp intergenic gap. A bidirectional promoter shares a short
    # intergenic region, typically under ~500 bp; 3.7 kb across a strand switch reads as
    # a separate transcriptional unit. They stay here.
    'PF00571': 'CBS', 'PF07228': 'SpoIIE', 'PF00149': 'Metallophos',
    'PF00702': 'Hydrolase', 'PF00370': 'FGGY_N', 'PF00581': 'Rhodanese',
    'PF04264': 'YceI', 'PF05899': 'Cupin_3',
}

CATEGORIES = {
    'core': CORE, 'tailoring': TAILORING, 'lipid': LIPID,
    'transport': TRANSPORT, 'regulation': REGULATION, 'mobile': MOBILE,
    'primary metabolism': PRIMARY_METABOLISM,
}

# Counted as "elaboration" — chemistry the cluster performs on its own product.
# Transport and regulation are deliberately excluded: they say a molecule is handled,
# not that it is elaborated, and the one confirmed bioactive small molecule here
# (pantaphos) encodes no transporter at all.
ELABORATION = ('tailoring', 'lipid')

_LOOKUP = {acc: cat for cat, d in CATEGORIES.items() for acc in d}
_NAMES = {acc: nm for d in CATEGORIES.values() for acc, nm in d.items()}


def bare(accession):
    """`PF07228.15` -> `PF07228`. antiSMASH writes versioned accessions, BiG-SCAPE
    writes bare ones, and the two have to join."""
    return (accession or '').split('.')[0].strip()


def category(accession):
    """Biosynthetic role for one Pfam accession; 'other' when not curated.

    'other' is not 'not biosynthetic' — the map covers ~93% of observed hits, and
    the tail is genuinely unclassified rather than known-irrelevant.
    """
    return _LOOKUP.get(bare(accession), 'other')


def name(accession):
    return _NAMES.get(bare(accession), bare(accession))


def profile(accessions):
    """{category: count} over an iterable of Pfam accessions."""
    out = {c: 0 for c in CATEGORIES}
    out['other'] = 0
    for a in accessions:
        out[category(a)] += 1
    return out


def elaboration(accessions):
    """How much chemistry the cluster does to its own product.

    Counts distinct tailoring/lipid domains, not occurrences: three copies of the
    same methyltransferase domain in one cluster is one kind of chemistry, and
    counting occurrences would let a tandem duplication look like elaboration.
    """
    return len({bare(a) for a in accessions if category(a) in ELABORATION})
