#!/usr/bin/env python3
"""Rank gene cluster families by how much they warrant laboratory follow-up.

The report currently marks a BGC "Potentially Novel" when KnownClusterBlast finds
no match. On Erwiniaceae that is **0 of 333 regions**, so every family carries the
flag and it orders nothing. MIBiG holds only a handful of phosphonate clusters, all
from actinomycetes, so a Pantoea BGC cannot match one — the axis is structurally
empty for this chemistry, not accidentally so.

What does vary is distance to the *characterised enzymes* the pipeline already
aligns against: assigned_pct_id spans 22.1-100.0 and pepm_pct_id 39.7-100.0.

Three axes, deliberately kept apart:

  DISTANCE   how far the chemistry sits from anything characterised
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

# The coupling enzyme is the branch point that sets the downstream pathway, so its
# divergence carries more weight than the pepM, which is the shared hallmark.
W_COUPLING, W_PEPM = 0.70, 0.30
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
        'd_coupling': [], 'd_pepm': [], 'classes': collections.Counter()})

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
            f['classes'][s.get('assigned_class', 'Unknown')] += 1
    return fam


def median(xs, default=0.0):
    return sorted(xs)[len(xs) // 2] if xs else default


def score(f):
    """(distance, evidence, priority, ...) — or None for distance when unclassifiable.

    A family whose coupling enzyme matches no known class has **no distance**, not a
    distance of zero. The first draft defaulted the missing median to 0.0 and ranked
    those families dead last, which is precisely backwards: an enzyme that resembles
    nothing characterised is the strongest novelty signal the data can carry. They are
    routed to a separate bucket instead of being scored.
    """
    if not f['d_coupling']:
        return None, None, None, None, None, (f['intact'] / f['seen'] if f['seen'] else 0.0), '(unclassified)'
    # DISTANCE — median over members, so one odd region cannot carry a family
    d_c, d_p = median(f['d_coupling']), median(f['d_pepm'])
    distance = W_COUPLING * d_c + W_PEPM * d_p

    # EVIDENCE — each term on 0..1
    n_gen = len(f['genomes'])
    t_genomes = min(1.0, math.log10(1 + n_gen) / math.log10(11))   # 1 -> .29, 10+ -> 1
    t_genera = min(1.0, max(0.0, (len(f['genera']) - 1) / 2))       # 3+ genera -> 1
    t_intact = f['intact'] / f['seen'] if f['seen'] else 0.0
    top = f['classes'].most_common(1)
    t_class = 1.0 if top and top[0][0] not in ('', 'Unknown') else 0.0
    evidence = (W_GENOMES * t_genomes + W_GENERA * t_genera +
                W_INTACT * t_intact + W_CLASS * t_class)

    return distance, evidence, distance * evidence, d_c, d_p, t_intact, top[0][0] if top else '-'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gcf_representatives', required=True)
    ap.add_argument('--tabulation', required=True)
    ap.add_argument('--coupling_support', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    reps, tab, sup = load(a.gcf_representatives, a.tabulation, a.coupling_support)
    fam = build(reps, tab, sup)

    rows, unclassified = [], []
    for fid, f in fam.items():
        dist, ev, pri, d_c, d_p, intact, cls = score(f)
        rec = dict(gcf=fid, members=f['n'], genomes=len(f['genomes']),
                   genera=len(f['genera']), intact=round(intact, 3),
                   coupling_class=cls)
        if pri is None:
            rec.update(distance='', evidence='', priority='', status='unclassified',
                       coupling_divergence='', pepm_divergence='')
            unclassified.append(rec)
        else:
            rec.update(distance=round(dist, 4), evidence=round(ev, 4),
                       priority=round(pri, 4), status='ranked',
                       coupling_divergence=round(d_c, 4), pepm_divergence=round(d_p, 4))
            rows.append(rec)

    rows.sort(key=lambda r: -r['priority'])
    unclassified.sort(key=lambda r: -r['members'])

    cols = ['rank', 'gcf', 'status', 'priority', 'distance', 'evidence',
            'coupling_divergence', 'pepm_divergence', 'members', 'genomes',
            'genera', 'intact', 'coupling_class']
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
