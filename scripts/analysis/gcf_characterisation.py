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

**The composite priority score has been withdrawn.** It was ISOLATION x EVIDENCE,
with isolation zeroed whenever a family's coupling enzyme matched a reference at
>=60% identity, and both halves turned out to be unsound on this data:

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

What remains is a table of per-family facts, emitted one row per family and sorted by
family id. Isolation and reference identity are still computed and published, each on
its own terms; nothing multiplies them into a rank. Revisit a composite only with
runs broad enough to calibrate one against.

  EVIDENCE is still reported -- how confident we are the family is real rather than an
  assembly artefact -- because it qualifies every other column on the row.
"""
import argparse, collections, csv, json, math, re, sys
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
    """{family_id: median nearest-cross-family BiG-SCAPE distance}, or None.

    BiG-SCAPE writes the complete all-pairs matrix, so this needs no extra compute:
    55,278 rows for 333 regions. For each member we take the smallest distance to any
    region *outside* its family, then the MEDIAN of those over the family -- median
    rather than min so one atypical member cannot make a family look connected, and
    not max so one cannot make it look isolated.

    **A member with no cross-family distance recorded is dropped, not scored 1.0.**
    It used to be initialised to 1.0 and left there, which is a sentinel wearing the
    costume of a measurement: 1.0 is the most isolated a family can look. That is
    harmless unpartitioned, where the matrix is complete -- measured on 1,302
    Enterobacterales regions, 0 of 1,302 members lack a cross-family pair. On a
    PARTITIONED run the merged table holds only within-partition distances, so a
    member whose nearest other family sits in another partition has none: 25 of 1,302
    members, inflating 61 members by a median 0.356, and taking four families to a
    flat 1.000 when the true maximum across all 72 is 0.668. Dropping them reports
    the median over what was actually measured, and a family with nothing measured
    returns None, which `score()` already handles as "no distances".
    """
    import sqlite3
    db = sqlite3.connect(db_path)
    fam_of = dict(db.execute(
        """SELECT bf.record_id, bf.family_id FROM bgc_record_family bf
           JOIN family f ON f.id = bf.family_id WHERE f.cutoff = ?""", (cutoff,)))
    if not fam_of:
        db.close()
        return {}
    nearest = {}
    for a, b, d in db.execute("SELECT record_a_id, record_b_id, distance FROM distance"):
        fa, fb = fam_of.get(a), fam_of.get(b)
        if fa is None or fb is None or fa == fb:
            continue
        for rec in (a, b):
            if d < nearest.get(rec, 2.0):
                nearest[rec] = d
    db.close()
    per = collections.defaultdict(list)
    missing = collections.Counter()
    for rec, fid in fam_of.items():
        if rec in nearest:
            per[fid].append(nearest[rec])
        else:
            missing[fid] += 1
    total_missing = sum(missing.values())
    if total_missing:
        # Loud, because the usual cause is a partitioned run whose merged distance
        # table is incomplete by design, and the symptom is families that look novel.
        print(f'warning: {total_missing} of {len(fam_of)} members have no '
              f'cross-family distance and are excluded from isolation; '
              f'{sum(1 for f in per if not per[f]) + len(set(missing) - set(per))} '
              f'famil(ies) have none at all. On a partitioned run this is expected '
              f'-- the merged table holds only within-partition pairs.',
              file=sys.stderr)
    out = {fid: median(v) for fid, v in per.items() if v}
    for fid in missing:
        out.setdefault(fid, None)
    return out


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
        'ref_pct': [], 'ref_org': collections.Counter(),
        # Which reference produced the identity. The report shows "75.3% (HvrC)";
        # without the name a reader cannot tell a Pantoea match from a Streptomyces
        # one, and four of the five families at 64-75% match HvrC while the fifth
        # matches FrbC.
        'ref_name': collections.Counter()})

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
                if s.get('assigned_ref'):
                    f['ref_name'][s['assigned_ref']] += 1
            f['classes'][s.get('assigned_class', 'Unknown')] += 1
    return fam


def median(xs, default=0.0):
    return sorted(xs)[len(xs) // 2] if xs else default


def score(f, isolation):
    """(evidence, divergences, reference facts) for one family.

    No composite: see the module docstring. `characterised` is still returned, as a
    LABEL on the reference-identity column rather than a gate on anything -- a reader
    seeing 100% to HvrC should be told that means pantaphos, but it must not silently
    erase the family's isolation.
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

    rows = []
    for fid, f in fam.items():
        (ev, d_c, d_p, intact, cls, best_ref,
         ref_org, ref_name, characterised) = score(f, iso.get(fid))
        # Reference columns are CONTEXT. The organism is there so a reader can see at
        # a glance that a 23% identity means "the only reference is a Streptomyces",
        # not "novel chemistry".
        if best_ref is None:
            ref_status = 'none for this class'
        elif characterised:
            ref_status = 'characterised'
        else:
            ref_status = 'distant only'
        # None, not 0.0: a family with no cross-family distance measured has no
        # isolation, and 0 would read as "sits on top of another family".
        isolation = iso.get(fid)
        rows.append(dict(
            gcf=fid, members=f['n'], genomes=len(f['genomes']),
            genera=len(f['genera']), intact=round(intact, 3), coupling_class=cls,
            isolation='' if isolation is None else round(isolation, 4),
            evidence=round(ev, 4),
            status='measured' if isolation is not None else 'no distances',
            reference_pct_id=round(best_ref, 1) if best_ref is not None else '',
            reference_name=ref_name, reference_organism=ref_org,
            reference_status=ref_status,
            coupling_divergence=round(d_c, 4) if d_c is not None else '',
            pepm_divergence=round(d_p, 4) if d_p is not None else ''))

    # By family id. There is no score to sort on, and the report sorts client-side.
    rows.sort(key=lambda r: int(r['gcf']) if str(r['gcf']).isdigit() else 0)

    cols = ['gcf', 'status', 'isolation', 'evidence', 'members', 'genomes', 'genera',
            'intact', 'coupling_class', 'reference_status', 'reference_pct_id',
            'reference_name', 'reference_organism', 'coupling_divergence',
            'pepm_divergence']
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

    no_iso = [r for r in rows if r['status'] != 'measured']
    print(f'{len(rows)} families, {len(rows) - len(no_iso)} with isolation measured')
    if no_iso:
        print(f'  no cross-family distances: '
              f'{", ".join(f"GCF-{r['gcf']}" for r in no_iso)}')
    far = sorted((r for r in rows if r['isolation'] != ''),
                 key=lambda r: -r['isolation'])[:3]
    for r in far:
        print(f'  most isolated: GCF-{r["gcf"]} {r["isolation"]:.3f}, '
              f'{r["members"]} members, {r["coupling_class"]}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
