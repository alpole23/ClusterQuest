#!/usr/bin/env python3
"""Rank gene cluster families by how much they warrant laboratory follow-up.

The report currently marks a BGC "Potentially Novel" when KnownClusterBlast finds
no match. On Erwiniaceae that is **0 of 333 regions**, so every family carries the
flag and it orders nothing. MIBiG holds only a handful of phosphonate clusters, all
from actinomycetes, so a Pantoea BGC cannot match one — the axis is structurally
empty for this chemistry, not accidentally so.

Distance to the *characterised enzymes* the pipeline aligns against was the first
attempt, and it does not work either. Six of the seven references are Streptomyces;
the one Enterobacterial reference (HvrC) is the pantaphos synthase. The resulting
identities are bimodal with an empty gap:

    22-45%  every Decarboxylase, Reductase and Transaminase family
    (nothing between 45.4% and 93.8%)
    94-100% every Synthase family

That is a binary readout of "does a same-taxon reference exist", not a novelty
gradient. It put five Reductase families at the top of the ranking purely because
VlpB is the most distant reference in the set, and within a class it is constant:
every Reductase member scored 0.75-0.78 regardless of its own sequence.

It was replaced by ISOLATION -- the median distance from each member to the nearest
region outside its family -- and that has now been withdrawn as well. The composite
was ISOLATION x EVIDENCE, with isolation zeroed whenever a family's coupling enzyme
matched a reference at >=60% identity, and both halves turned out to be unsound on
this data:

  the gate      CHARACTERISED_PCT = 60 was chosen because measured identities were
                "bimodal with nothing between 45.4% and 93.8%" -- true of Erwiniaceae,
                not of Enterobacterales, where they run 18.5-100% with the widest gap
                at 75-94%. Five families sit at 64-75%, inside the supposedly empty
                zone, and were zeroed out of follow-up while no characterised CLUSTER
                sits within the family cutoff of any of them (0.38-0.63 away).
  isolation     a within-run measurement, so it says whether anything else in THIS
                dataset resembles the family -- not whether the family is novel. It is
                not comparable between runs of different scope, and adding genomes can
                only lower it.

**Isolation was then dropped from the output entirely (2026-10-29).** Keeping it as a
published column assumed a reader would carry the caveat above alongside the number,
and in practice a column headed "isolation" sitting beside the reference identity
reads as a novelty score whatever the tooltip says. It measured something real but
not the thing anyone wanted it for, and it cost the process an input -- the whole
BiG-SCAPE database -- that nothing else here needs. The distances it was computed
from are not lost: they are in the database, and the GCF network in the report draws
the same relationships in a form that cannot be mistaken for a ranking.

What remains is a table of per-family facts, emitted one row per family and sorted by
family id. Reference identity is published as context, labelled with the organism it
came from; nothing multiplies anything into a rank. Revisit a composite only with
runs broad enough to calibrate one against.

  EVIDENCE is still reported -- how confident we are the family is real rather than an
  assembly artefact -- because it qualifies every other column on the row.
"""
# GTDB ranks, shallowest first. The order IS the semantics: span_rank reports the
# first rank at which every member of a family agrees, so "genus" means the family
# does not leave one genus and "phylum" means it crosses orders and classes.
RANK_CODES = [('species', 's'), ('genus', 'g'), ('family', 'f'),
              ('order', 'o'), ('class', 'c'), ('phylum', 'p')]

import argparse, collections, csv, json, math, os, re, sys
from pathlib import Path

# record_id is spelled inconsistently across the pipeline's own outputs — some rows
# carry the accession version (JARNMU010000002.1), some do not (JABDZE010000001).
# Normalise both sides of every join rather than trusting either spelling.
_VERSION = re.compile(r'\.\d+$')


def bare_accession(acc):
    return _VERSION.sub('', acc)

# Identity above which a reference is LABELLED as characterising the family's
# chemistry. It gates nothing -- see the module docstring. The value was chosen when
# identities looked bimodal with nothing between 45.4% and 93.8%, which held on
# Erwiniaceae and does not hold here: Enterobacterales runs 18.5-100% with its widest
# gap at 75-94%, and five families sit at 64-75%. So treat the label as advisory and
# read `reference_pct_id` itself.
CHARACTERISED_PCT = 60.0
# Evidence weights: independent observation dominates, structural integrity next.
W_GENOMES, W_GENERA, W_INTACT, W_CLASS = 0.35, 0.25, 0.30, 0.10


def load(reps_path, tab_path, sup_path):
    reps = json.loads(Path(reps_path).read_text())
    tab = list(csv.DictReader(Path(tab_path).open(), delimiter='\t'))
    # The coupling support is optional upstream: GCF_BIOSYNTHETIC_TREE emits it
    # `optional: true`, so this may be a placeholder. Degrade to no coupling
    # evidence rather than crashing -- every other column still has its data,
    # and this table is what replaced the GCF x species heatmap.
    sp = Path(sup_path)
    if not sp.exists() or sp.name.startswith('NO_'):
        sys.stderr.write(f'note: no coupling support ({sp.name}); '
                         'coupling_class and the divergence columns will be blank\n')
        return reps, tab, {}
    lines = [l for l in sp.open() if not l.startswith('#')]
    # Support rows are keyed "ACCESSION.VERSION.regionNNN", with extra "_N" rows for the
    # sub-records (protocluster, cand_cluster, proto_core). Keep region level only, and
    # index on the *unversioned* accession: region_tabulation.tsv has already lost the
    # version, and it is not recoverable from there.
    sup = {}
    for r in csv.DictReader(lines, delimiter='\t'):
        bgc = r['bgc']
        if '.region' not in bgc:
            continue
        head, region = bgc.rsplit('.region', 1)
        if '_' in region:                       # a sub-record, not the region itself
            continue
        acc = bare_accession(head)              # drop the version if present
        sup[(acc, int(region))] = r
    return reps, tab, sup


def num(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def build(reps, tab, sup, gtdbtk_summary=None):
    """One record per family, joining the three sources on region identity."""
    # tabulation is the only place carrying both spellings of a region
    by_region = {}
    for r in tab:
        # Do NOT strip an extension here: `file` carries none, and 25 of 333 genome
        # names legitimately contain a dot (accession-style ..._GCA_963520565.1).
        genome = r['file']
        key = f"{genome}|{r['region_name']}"
        by_region[key] = {'sup': sup.get((bare_accession(r['record_id']), int(r['region']))),
                          'edge': r['contig_edge'].strip().lower() in ('true', '1', 'yes'),
                          'genus': genome.split('_')[0]}

    # GTDB taxonomy per genome, when GTDB-Tk ran. This REPLACES nothing: the
    # existing `genera` count stays as it is, derived from the genome name's
    # first token, because it feeds the evidence score and changing its source
    # would silently move a published number.
    #
    # It is kept alongside rather than trusted, because that shortcut is wrong
    # for any name whose first token is not a genus -- every "Candidatus ..."
    # collapses to the genus "Candidatus". Measured: 125 of 1,530 Enterobacterales
    # genomes (8.2%), and the shortcut finds 417 genera where GTDB finds 476.
    gtdb = {}
    if gtdbtk_summary and os.path.exists(gtdbtk_summary):
        with open(gtdbtk_summary) as fh:
            next(fh, None)
            for line in fh:
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                name = parts[0][4:] if parts[0].startswith('usr_') else parts[0]
                lin = {}
                for rank, code in RANK_CODES:
                    m = re.search(rf'{code}__([A-Za-z0-9_.\- ]+)', parts[1])
                    lin[rank] = m.group(1).strip() if m and m.group(1).strip() else None
                gtdb[name] = lin

    fam = collections.defaultdict(lambda: {
        'n': 0, 'genomes': set(), 'genera': set(), 'intact': 0, 'seen': 0,
        'd_coupling': [], 'd_pepm': [], 'classes': collections.Counter(),
        'ref_pct': [], 'ref_org': collections.Counter(),
        # Which reference produced the identity. The report shows "75.3% (HvrC)";
        # without the name a reader cannot tell a Pantoea match from a Streptomyces
        # one, and four of the five families at 64-75% match HvrC while the fifth
        # matches FrbC.
        'ref_name': collections.Counter(),
        'lineages': []})

    for key, v in reps['bgc_to_gcf'].items():
        f = fam[v['family_id']]
        f['n'] += 1
        f['genomes'].add(key.rsplit('|', 1)[0])
        meta = by_region.get(key)
        if not meta:
            continue                      # region absent from the tabulation join
        f['seen'] += 1
        f['genera'].add(meta['genus'])
        lin = gtdb.get(key.rsplit('|', 1)[0])
        if lin:
            f['lineages'].append(lin)
        f['intact'] += (not meta['edge'])
        s = meta['sup']
        if s:
            a, p = num(s.get('assigned_pct_id')), num(s.get('pepm_pct_id'))
            if a is not None:
                f['d_coupling'].append(1 - a / 100)
            if p is not None:
                f['d_pepm'].append(1 - p / 100)
            if a is not None:
                f['ref_pct'].append(a)
                if s.get('assigned_ref_organism'):
                    f['ref_org'][s['assigned_ref_organism']] += 1
                if s.get('assigned_ref'):
                    f['ref_name'][s['assigned_ref']] += 1
            f['classes'][s.get('assigned_class', 'Unknown')] += 1
    return fam


def median(xs, default=0.0):
    return sorted(xs)[len(xs) // 2] if xs else default


def score(f):
    """(evidence, divergences, reference facts) for one family.

    No composite: see the module docstring. `characterised` is still returned, as a
    LABEL on the reference-identity column rather than a gate on anything -- a reader
    seeing 100% to HvrC should be told that means pantaphos, but it must not gate the
    family out of anything.
    """
    best_ref = max(f['ref_pct']) if f['ref_pct'] else None
    characterised = best_ref is not None and best_ref >= CHARACTERISED_PCT

    d_c = median(f['d_coupling']) if f['d_coupling'] else None
    d_p = median(f['d_pepm']) if f['d_pepm'] else None

    # EVIDENCE — each term on 0..1
    n_gen = len(f['genomes'])
    t_genomes = min(1.0, math.log10(1 + n_gen) / math.log10(11))   # 1 -> .29, 10+ -> 1
    t_genera = min(1.0, max(0.0, (len(f['genera']) - 1) / 2))       # 3+ genera -> 1
    t_intact = f['intact'] / f['seen'] if f['seen'] else 0.0
    top = f['classes'].most_common(1)
    t_class = 1.0 if top and top[0][0] not in ('', 'Unknown') else 0.0
    evidence = (W_GENOMES * t_genomes + W_GENERA * t_genera +
                W_INTACT * t_intact + W_CLASS * t_class)

    return (evidence, d_c, d_p, t_intact, top[0][0] if top else '-', best_ref,
            f['ref_org'].most_common(1)[0][0] if f['ref_org'] else '',
            f['ref_name'].most_common(1)[0][0] if f['ref_name'] else '',
            characterised)


def _span(lineages):
    """Shallowest GTDB rank at which every member of a family agrees.

    Blank for a single-member family, where the question is meaningless. Reported
    rather than a boolean because "crosses a genus boundary" and "crosses a phylum
    boundary" are very different claims, and the ablation in
    docs/comparisons/ shows the genus-level answer is sensitive to how much
    flanking sequence BiG-SCAPE is shown -- so the rank carries the nuance a flag
    would throw away.
    """
    if len(lineages) < 2:
        return ''
    for rank, _ in RANK_CODES:
        vals = {l[rank] for l in lineages}
        if len(vals) == 1 and None not in vals:
            return rank
    return 'above phylum'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gcf_representatives', required=True)
    ap.add_argument('--tabulation', required=True)
    ap.add_argument('--coupling_support', required=True)
    ap.add_argument('--gtdbtk_summary', default=None,
                    help='GTDB-Tk bac120 summary; adds the gtdb_* taxonomy columns. '
                         'Optional -- the table is written without them if absent.')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    reps, tab, sup = load(a.gcf_representatives, a.tabulation, a.coupling_support)
    fam = build(reps, tab, sup, a.gtdbtk_summary)

    rows = []
    for fid, f in fam.items():
        (ev, d_c, d_p, intact, cls, best_ref,
         ref_org, ref_name, characterised) = score(f)
        # Reference columns are CONTEXT. The organism is there so a reader can see at
        # a glance that a 23% identity means "the only reference is a Streptomyces",
        # not "novel chemistry".
        if best_ref is None:
            ref_status = 'none for this class'
        elif characterised:
            ref_status = 'characterised'
        else:
            ref_status = 'distant only'
        rows.append(dict(
            gcf=fid, members=f['n'], genomes=len(f['genomes']),
            genera=len(f['genera']), intact=round(intact, 3), coupling_class=cls,
            evidence=round(ev, 4),
            reference_pct_id=round(best_ref, 1) if best_ref is not None else '',
            reference_name=ref_name, reference_organism=ref_org,
            reference_status=ref_status,
            gtdb_genera=len({l['genus'] for l in f['lineages'] if l['genus']}) or '',
            gtdb_phyla=len({l['phylum'] for l in f['lineages'] if l['phylum']}) or '',
            span_rank=_span(f['lineages']),
            gtdb_genus_list='; '.join(sorted({l['genus'] for l in f['lineages']
                                              if l['genus']})) or '',
            coupling_divergence=round(d_c, 4) if d_c is not None else '',
            pepm_divergence=round(d_p, 4) if d_p is not None else ''))

    # By family id. There is no score to sort on, and the report sorts client-side.
    rows.sort(key=lambda r: int(r['gcf']) if str(r['gcf']).isdigit() else 0)

    # `genera` is the first-token count that feeds the evidence score; the gtdb_*
    # columns are the corrected taxonomy and are deliberately separate, so a reader
    # can see both rather than wonder which one a single column meant.
    cols = ['gcf', 'evidence', 'members', 'genomes', 'genera',
            'gtdb_genera', 'gtdb_phyla', 'span_rank', 'intact', 'coupling_class',
            'reference_status', 'reference_pct_id', 'reference_name',
            'reference_organism', 'coupling_divergence', 'pepm_divergence',
            'gtdb_genus_list']
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

    print(f'{len(rows)} families')
    big = sorted(rows, key=lambda r: -r['members'])[:3]
    for r in big:
        print(f'  largest: GCF-{r["gcf"]} {r["members"]} members, '
              f'{r["genomes"]} genomes, {r["coupling_class"] or "unclassified"}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
