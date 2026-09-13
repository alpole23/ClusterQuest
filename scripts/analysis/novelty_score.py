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

What replaces it is ISOLATION. BiG-SCAPE already computes the full all-pairs
distance matrix over the run -- 55,278 pairs for 333 regions, no reference set
involved. Isolation asks the question the reference set cannot bias: is there
anything else in this run like this family? Measured on Erwiniaceae the two axes are
essentially uncorrelated (Spearman rho = -0.17), and isolation is not a family-size
artefact (singletons mean 0.490, multi-member families 0.505, singletons spanning
the full range).

Three axes, deliberately kept apart:

  DISTANCE   isolation in BiG-SCAPE space, zeroed when the family's chemistry is
             already characterised
  EVIDENCE   how confident we are the family is real, not an assembly artefact
  REACH      how easy it would be to obtain a strain (reported, not ranked on)

Rank on DISTANCE x EVIDENCE. They multiply rather than add because both are
necessary: a maximally divergent single truncated region is not a lead, and adding
would let novelty compensate for having no evidence behind it.
"""
import argparse, collections, csv, json, math, re, sys
from pathlib import Path

# record_id is spelled inconsistently across the pipeline's own outputs — some rows
# carry the accession version (JARNMU010000002.1), some do not (JABDZE010000001).
# Normalise both sides of every join rather than trusting either spelling.
_VERSION = re.compile(r'\.\d+$')


def bare_accession(acc):
    return _VERSION.sub('', acc)

# Identity above which a reference counts as characterising the family's chemistry.
# The exact value is unimportant and that is the point: measured identities are
# bimodal with NOTHING between 45.4% and 93.8%, so any threshold in that gap gives
# identical results. This is a robust separator, not a tuned one.
CHARACTERISED_PCT = 60.0
# Evidence weights: independent observation dominates, structural integrity next.
W_GENOMES, W_GENERA, W_INTACT, W_CLASS = 0.35, 0.25, 0.30, 0.10


def load(reps_path, tab_path, sup_path):
    reps = json.loads(Path(reps_path).read_text())
    tab = list(csv.DictReader(Path(tab_path).open(), delimiter='\t'))
    lines = [l for l in Path(sup_path).open() if not l.startswith('#')]
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


def isolation_by_family(db_path, cutoff):
    """{family_id: median nearest-cross-family BiG-SCAPE distance}.

    BiG-SCAPE writes the complete all-pairs matrix, so this needs no extra compute:
    55,278 rows for 333 regions. For each member we take the smallest distance to any
    region *outside* its family, then the MEDIAN of those over the family -- median
    rather than min so one atypical member cannot make a family look connected, and
    not max so one cannot make it look isolated.
    """
    import sqlite3
    db = sqlite3.connect(db_path)
    fam_of = dict(db.execute(
        """SELECT bf.record_id, bf.family_id FROM bgc_record_family bf
           JOIN family f ON f.id = bf.family_id WHERE f.cutoff = ?""", (cutoff,)))
    if not fam_of:
        db.close()
        return {}
    nearest = {r: 1.0 for r in fam_of}
    for a, b, d in db.execute("SELECT record_a_id, record_b_id, distance FROM distance"):
        fa, fb = fam_of.get(a), fam_of.get(b)
        if fa is None or fb is None or fa == fb:
            continue
        if d < nearest[a]:
            nearest[a] = d
        if d < nearest[b]:
            nearest[b] = d
    db.close()
    per = collections.defaultdict(list)
    for rec, fid in fam_of.items():
        per[fid].append(nearest[rec])
    return {fid: median(v) for fid, v in per.items()}


def num(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def build(reps, tab, sup):
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

    fam = collections.defaultdict(lambda: {
        'n': 0, 'genomes': set(), 'genera': set(), 'intact': 0, 'seen': 0,
        'd_coupling': [], 'd_pepm': [], 'classes': collections.Counter(),
        'ref_pct': [], 'ref_org': collections.Counter()})

    for key, v in reps['bgc_to_gcf'].items():
        f = fam[v['family_id']]
        f['n'] += 1
        f['genomes'].add(key.rsplit('|', 1)[0])
        meta = by_region.get(key)
        if not meta:
            continue                      # region absent from the tabulation join
        f['seen'] += 1
        f['genera'].add(meta['genus'])
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
            f['classes'][s.get('assigned_class', 'Unknown')] += 1
    return fam


def median(xs, default=0.0):
    return sorted(xs)[len(xs) // 2] if xs else default


def score(f, isolation):
    """(distance, evidence, priority, ...) for one family.

    DISTANCE is isolation, gated by whether the chemistry is already characterised.
    A family sitting 100% identical to HvrC is pantaphos: isolated or not, it is a
    solved cluster and scores zero. Below the gap it is not characterised at any
    useful resolution and isolation carries the signal alone.

    Note what this fixes beyond the bias. The previous version could not rank a
    family whose coupling enzyme matched no known class at all -- it had no distance
    to compute, so those families went to an unranked bucket. Isolation needs no
    reference, so they rank normally now; on Erwiniaceae the single MOST isolated
    family in the run (GCF-16, isolation 0.861) was one of them.
    """
    intact = f['intact'] / f['seen'] if f['seen'] else 0.0
    if isolation is None:
        # No clustering distances for this family: nothing honest to rank on.
        return None, None, None, None, None, intact, '(no distances)'

    best_ref = max(f['ref_pct']) if f['ref_pct'] else None
    characterised = best_ref is not None and best_ref >= CHARACTERISED_PCT
    distance = 0.0 if characterised else isolation

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

    return (distance, evidence, distance * evidence, d_c, d_p, t_intact,
            top[0][0] if top else '-', best_ref,
            f['ref_org'].most_common(1)[0][0] if f['ref_org'] else '',
            characterised)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gcf_representatives', required=True)
    ap.add_argument('--tabulation', required=True)
    ap.add_argument('--coupling_support', required=True)
    ap.add_argument('--bigscape_db', required=True,
                    help='BiG-SCAPE sqlite; supplies the all-pairs distances '
                         'that isolation is computed from')
    ap.add_argument('--cutoff', type=float, default=0.3)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    reps, tab, sup = load(a.gcf_representatives, a.tabulation, a.coupling_support)
    fam = build(reps, tab, sup)
    iso = isolation_by_family(a.bigscape_db, a.cutoff)
    if not iso:
        print(f'warning: no families at cutoff {a.cutoff} in {a.bigscape_db}; '
              f'nothing can be ranked', file=sys.stderr)

    rows, unclassified = [], []
    for fid, f in fam.items():
        (dist, ev, pri, d_c, d_p, intact, cls,
         best_ref, ref_org, characterised) = score(f, iso.get(fid))
        rec = dict(gcf=fid, members=f['n'], genomes=len(f['genomes']),
                   genera=len(f['genera']), intact=round(intact, 3),
                   coupling_class=cls)
        if pri is None:
            rec.update(priority='', distance='', isolation='', evidence='',
                       status='no distances', reference_pct_id='',
                       reference_organism='', reference_status='',
                       coupling_divergence='', pepm_divergence='')
            unclassified.append(rec)
        else:
            # Reference columns are published as CONTEXT, not as ranking terms. The
            # organism is there so a reader can see at a glance that a 23% identity
            # means "the only reference is a Streptomyces", not "novel chemistry".
            if best_ref is None:
                ref_status = 'none for this class'
            elif characterised:
                ref_status = 'characterised'
            else:
                ref_status = 'distant only'
            rec.update(priority=round(pri, 4), distance=round(dist, 4),
                       isolation=round(iso.get(fid, 0.0), 4),
                       evidence=round(ev, 4), status='ranked',
                       reference_pct_id=round(best_ref, 1) if best_ref is not None else '',
                       reference_organism=ref_org, reference_status=ref_status,
                       coupling_divergence=round(d_c, 4) if d_c is not None else '',
                       pepm_divergence=round(d_p, 4) if d_p is not None else '')
            rows.append(rec)

    rows.sort(key=lambda r: -r['priority'])
    unclassified.sort(key=lambda r: -r['members'])

    cols = ['rank', 'gcf', 'status', 'priority', 'distance', 'isolation', 'evidence',
            'members', 'genomes', 'genera', 'intact', 'coupling_class',
            'reference_status', 'reference_pct_id', 'reference_organism',
            'coupling_divergence', 'pepm_divergence']
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        # Unclassifiable families lead the file: they are not ranked, and burying them
        # under a score they do not have is how they would get overlooked.
        for r in unclassified:
            w.writerow({**r, 'rank': ''})
        for i, r in enumerate(rows, 1):
            w.writerow({**r, 'rank': i})

    print(f'{len(rows)} families ranked, {len(unclassified)} unclassifiable')
    if rows:
        t = rows[0]
        print(f'  top: GCF-{t["gcf"]} priority {t["priority"]:.3f} '
              f'(distance {t["distance"]:.2f} x evidence {t["evidence"]:.2f}), '
              f'{t["members"]} members, {t["coupling_class"]}')
    for r in unclassified:
        print(f'  unclassifiable: GCF-{r["gcf"]}, {r["members"]} members '
              f'across {r["genomes"]} genomes')
    return 0


if __name__ == '__main__':
    sys.exit(main())
