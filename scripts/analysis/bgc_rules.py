#!/usr/bin/env python3
"""Apply the phosphonate BGC rule set, and check the rules against what is known.

Rules live in `assets/phosphonate_rules.json` and are edited there, not here. This
script reads them, applies them to each region, and -- the part that matters --
verifies that no rule contradicts a cluster whose chemistry is known from lab work.

Why that check exists. Three separate rules for predicting phosphonolipids have been
proposed for this pipeline and all three were INVERTED on the two characterised
clusters:

  1. TPP + NTP_transf_3 => lipid.  The confirmed lipid carries no NTP_transf at all.
  2. Low elaboration    => lipid.  Reversed once serine hydroxymethyltransferase
                                   stopped counting as a tailoring methyltransferase.
  3. PF01066 present    => lipid.  It sits in the confirmed NON-lipid's core operon.

Each looked right until tested. `--validate` runs every rule against `known_clusters`
and fails loudly on a contradiction, so the fourth one cannot ship quietly.

Usage:
    python scripts/analysis/bgc_rules.py --validate
    python scripts/analysis/bgc_rules.py --antismash results/antismash_results/Taxon \\
        --db results/bigscape_results/Taxon/Taxon.db --out rules_report.tsv
"""
import argparse
import collections
import csv
import json
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.domain_functions import category  # noqa: E402

DEFAULT_RULES = Path(__file__).resolve().parents[2] / 'assets' / 'phosphonate_rules.json'


# ─── Rule file ───────────────────────────────────────────────────────────────

def load_rules(path):
    spec = json.loads(Path(path).read_text())
    for key in ('domain_sets', 'rules', 'known_clusters'):
        if key not in spec:
            sys.exit(f'{path}: missing required key "{key}"')
    for r in spec['rules']:
        for field in ('id', 'tier', 'verdict', 'summary', 'rationale', 'source', 'added'):
            if not r.get(field):
                sys.exit(f'{path}: rule "{r.get("id", "?")}" is missing "{field}". '
                         f'Every rule must record why it is believed and where that '
                         f'belief comes from.')
    ids = [r['id'] for r in spec['rules']]
    dupes = [i for i, c in collections.Counter(ids).items() if c > 1]
    if dupes:
        sys.exit(f'{path}: duplicate rule ids {dupes}')
    return spec


def compile_sets(spec):
    """{name: (set-of-accessions, compiled product regex or None)}"""
    out = {}
    for name, d in spec['domain_sets'].items():
        pats = d.get('product_patterns') or []
        rx = re.compile('|'.join(f'(?:{p})' for p in pats), re.I) if pats else None
        out[name] = (set(d.get('accessions') or []), set(d.get('smcog') or []), rx)
    return out


# ─── Matching ────────────────────────────────────────────────────────────────

def gene_matches(gene, sets, names):
    """Does this gene match any of the named domain sets?

    A set matches on a Pfam accession, an SMCOG id, OR a product-name pattern. The
    product path is not redundant: the IS630 transposase bounding the LMG 5342
    phosphonolipid cluster carries no Pfam domain at all, and on accessions alone the
    boundary is missed and three housekeeping genes are pulled inside the cluster.
    """
    if isinstance(names, str):
        names = [names]
    for n in names:
        accs, smcogs, rx = sets.get(n, (set(), set(), None))
        if accs & set(gene['accessions']):
            return True
        if smcogs & set(gene.get('smcog') or []):
            return True
        if rx and gene.get('product') and rx.search(gene['product']):
            return True
    return False


def evaluate(region, spec, sets):
    """Apply every rule to one region. Returns {rule_id: True/False/None}."""
    genes = region['genes']
    res = {}
    for r in spec['rules']:
        t = r.get('test') or {}
        verdict = None
        if 'any_domain' in t:
            verdict = any(gene_matches(g, sets, t['any_domain']) for g in genes)
        elif 'no_domain' in t:
            verdict = not any(gene_matches(g, sets, t['no_domain']) for g in genes)
        elif 'all_domain_sets' in t:
            verdict = all(any(gene_matches(g, sets, n) for g in genes)
                          for n in t['all_domain_sets'])
        elif 'min_role_count' in t:
            counts = collections.Counter(g['role'] for g in region.get('inside', genes))
            verdict = all(counts.get(k, 0) >= v for k, v in t['min_role_count'].items())
        elif 'not_dominated_by' in t:
            n = sum(1 for g in genes if gene_matches(g, sets, t['not_dominated_by']))
            verdict = n < max(1, len(genes) // 2)
        elif 'contig_edge' in t:
            ce = region.get('contig_edge')
            verdict = None if ce is None else (ce == t['contig_edge'])
        res[r['id']] = verdict
    return res


def apply_boundary(region, sets):
    """Narrow a region to its biosynthetic span.

    Mobile elements first: phosphonate BGCs here are horizontally acquired and the
    elements that moved them sit at the edges, which is a tighter and more meaningful
    boundary than antiSMASH's generous region. Falls back to excluding beyond the
    first essential housekeeping gene, then to the whole region.
    """
    genes = region['genes']
    core = [i for i, g in enumerate(genes) if g['role'] == 'core']
    if not core:
        region['inside'] = genes
        region['boundary'] = 'none (no core gene)'
        return region
    lo_i, hi_i = min(core), max(core)

    mob = [i for i, g in enumerate(genes) if gene_matches(g, sets, 'mobile_element')]
    left = max((i for i in mob if i < lo_i), default=None)
    right = min((i for i in mob if i > hi_i), default=None)
    how = []
    if left is not None or right is not None:
        how.append('mobile element')
    else:
        ess = [i for i, g in enumerate(genes) if gene_matches(g, sets, 'essential_housekeeping')]
        left = max((i for i in ess if i < lo_i), default=None)
        right = min((i for i in ess if i > hi_i), default=None)
        if left is not None or right is not None:
            how.append('essential housekeeping')

    a = (left + 1) if left is not None else 0
    b = (right - 1) if right is not None else len(genes) - 1
    region['inside'] = genes[a:b + 1]
    region['boundary'] = (' / '.join(how) if how else 'antiSMASH region (no marker found)')
    region['excluded'] = len(genes) - len(region['inside'])
    return region


# ─── Validation against known clusters ───────────────────────────────────────

def validate(spec, regions=None):
    """Check the file is coherent, and that no rule contradicts a known cluster."""
    problems, notes = [], []

    known = {k['region']: k for k in spec['known_clusters']}
    # A product-class rule must rest on examples that cannot share an explanation by
    # descent. Three lipid rules each fit ONE cluster and inverted; the failure mode is
    # fitting a single example, so the guard counts independent lineages. The rule names
    # its supporting clusters explicitly -- inferring them from the rule id was fragile.
    pc = spec.get('product_class_rules', {})
    need = (pc.get('admission_test') or {}).get('min_independent_lineages', 2)
    by_id = {k['id']: k for k in spec['known_clusters']}
    for r in spec['rules']:
        if r.get('tier') != 'product_class':
            continue
        sup = r.get('supported_by')
        if not sup:
            problems.append(
                f'rule "{r["id"]}" is a product-class rule and must list the '
                f'known_clusters it rests on in "supported_by".')
            continue
        missing = [c for c in sup if c not in by_id]
        if missing:
            problems.append(f'rule "{r["id"]}" cites unknown cluster(s) {missing}')
            continue
        lineages = {by_id[c]['organism'].split()[0] for c in sup}
        if len(lineages) < need:
            problems.append(
                f'rule "{r["id"]}" rests on {len(lineages)} independent lineage(s) '
                f'({", ".join(sorted(lineages))}); {need} required. Three prior rules '
                f'fit a single example and then inverted.')

    if regions:
        for label, reg in regions.items():
            if label not in known:
                continue
            k = known[label]
            res = evaluate(reg, spec, compile_sets(spec))
            for r in spec['rules']:
                if r['verdict'] != 'required':
                    continue
                if res.get(r['id']) is False and k['truth'].get('is_phosphonate_bgc'):
                    problems.append(
                        f'rule "{r["id"]}" fails on {label} ({k["organism"]}), which is '
                        f'a confirmed phosphonate BGC. Either the rule or the record '
                        f'is wrong.')
            notes.append(f'  {label}  {k["organism"]}: ' +
                         ', '.join(f'{i}={"pass" if v else "FAIL" if v is False else "n/a"}'
                                   for i, v in res.items() if v is not None))
    return problems, notes


# ─── Region loading ──────────────────────────────────────────────────────────

def load_regions(antismash_dir, db_path, cutoff=0.3):
    from Bio import SeqIO
    db = sqlite3.connect(db_path)
    members = []
    for path, fid in db.execute(
            """SELECT g.path, bf.family_id FROM bgc_record_family bf
               JOIN bgc_record r ON r.id=bf.record_id JOIN gbk g ON g.id=r.gbk_id
               JOIN family f ON f.id=bf.family_id WHERE f.cutoff=?""", (cutoff,)):
        p = Path(path)
        members.append((p.parent.name, p.name, fid))
    out = {}
    for genome, fname, fid in members:
        path = Path(antismash_dir) / genome / fname
        if not path.exists():
            continue
        label = fname[:-4]
        doms, smcog = collections.defaultdict(list), collections.defaultdict(list)
        edge = None
        for rec in SeqIO.parse(str(path), 'genbank'):
            for f in rec.features:
                if f.type == 'PFAM_domain':
                    tag = (f.qualifiers.get('locus_tag') or ['?'])[0]
                    for x in f.qualifiers.get('db_xref', []):
                        if x.startswith('PF'):
                            doms[tag].append(x.split('.')[0])
                elif f.type == 'region':
                    ce = (f.qualifiers.get('contig_edge') or [''])[0]
                    edge = ce.strip().lower() == 'true'
        genes = []
        for rec in SeqIO.parse(str(path), 'genbank'):
            for f in rec.features:
                if f.type != 'CDS':
                    continue
                q = f.qualifiers
                tag = (q.get('locus_tag') or q.get('gene') or ['?'])[0]
                accs = doms.get(tag, [])
                sm = re.findall(r'(SMCOG\d+)', ' '.join(q.get('gene_functions', [])))
                roles = {category(a) for a in accs} - {'other'}
                genes.append(dict(tag=tag, product=(q.get('product') or [''])[0],
                                  accessions=accs, smcog=sm,
                                  start=int(f.location.start), end=int(f.location.end),
                                  strand=1 if f.location.strand in (1, None) else -1,
                                  role=sorted(roles)[0] if roles else 'other'))
        genes.sort(key=lambda g: g['start'])
        out[label] = dict(label=label, genome=genome, family=fid,
                          genes=genes, contig_edge=edge)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rules', type=Path, default=DEFAULT_RULES)
    ap.add_argument('--antismash', type=Path)
    ap.add_argument('--db', type=Path)
    ap.add_argument('--cutoff', type=float, default=0.3)
    ap.add_argument('--out', type=Path)
    ap.add_argument('--validate', action='store_true',
                    help='check the rule file, and the rules against known clusters')
    a = ap.parse_args()

    spec = load_rules(a.rules)
    sets = compile_sets(spec)
    print(f'{a.rules.name}: v{spec["version"]}, {len(spec["rules"])} rules, '
          f'{len(spec["known_clusters"])} known clusters, '
          f'{len(spec.get("product_class_rules", {}).get("rules", []))} product-class rules')

    regions = None
    if a.antismash and a.db:
        regions = load_regions(a.antismash, a.db, a.cutoff)
        for reg in regions.values():
            apply_boundary(reg, sets)
        print(f'loaded {len(regions)} regions')

    if a.validate:
        problems, notes = validate(spec, regions)
        if notes:
            print('\nknown clusters:')
            print('\n'.join(notes))
        if problems:
            print('\nPROBLEMS:')
            for p in problems:
                print(f'  - {p}')
            sys.exit(1)
        print('\nvalidation passed — no rule contradicts a known cluster')
        if not regions:
            print('  (rule file only; pass --antismash/--db to test against real regions)')
        return

    if not regions:
        ap.error('--antismash and --db are required unless --validate is used alone')

    rows = []
    for label, reg in sorted(regions.items()):
        res = evaluate(reg, spec, sets)
        required = [r['id'] for r in spec['rules'] if r['verdict'] == 'required']
        excl = [r['id'] for r in spec['rules'] if r['verdict'] == 'excludes']
        passes = all(res.get(i) for i in required)
        blocked = any(res.get(i) is False for i in excl)
        rows.append(dict(
            region=label, genome=reg['genome'], family=reg['family'],
            verdict=('phosphonate BGC' if passes and not blocked else
                     'excluded' if blocked else 'incomplete'),
            boundary=reg['boundary'], genes=len(reg['genes']),
            genes_inside=len(reg['inside']), excluded=reg.get('excluded', 0),
            contig_edge=reg['contig_edge'],
            **{i: ('' if res[i] is None else ('pass' if res[i] else 'fail')) for i in res}))

    cols = (['region', 'genome', 'family', 'verdict', 'boundary', 'genes',
             'genes_inside', 'excluded', 'contig_edge'] +
            [r['id'] for r in spec['rules']])
    if a.out:
        with a.out.open('w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        print(f'wrote {a.out}')
    v = collections.Counter(r['verdict'] for r in rows)
    b = collections.Counter(r['boundary'] for r in rows)
    print('\nverdicts : ' + ', '.join(f'{k} {n}' for k, n in v.most_common()))
    print('boundary : ' + ', '.join(f'{k} {n}' for k, n in b.most_common()))


if __name__ == '__main__':
    main()
