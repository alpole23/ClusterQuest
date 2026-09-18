#!/usr/bin/env python3
"""Compare two pipeline runs BGC by BGC, so a change's effect is data rather than memory.

Each change to detection or clustering moves the numbers the report rests on, and the
runs that show it are large, gitignored and eventually deleted. This extracts the small
comparable core of two result directories — what each BGC was called, how many genes it
had, which family it landed in — into files small enough to keep in the repository.

    python scripts/compare_runs.py --results results \\
        --before Erwiniaceae_pre_recovery --after Erwiniaceae \\
        --outdir docs/comparisons/orf_recovery

BGCs are matched on genome plus region filename, which is stable across runs; antiSMASH
stamps a run date into every region GenBank, so file hashes are not (see CLAUDE.md).

Writes per_bgc.tsv, per_gcf.tsv and summary.json.
"""
import argparse
import json
import sqlite3
from collections import defaultdict
from pathlib import Path


def region_files(antismash_dir):
    """{(genome, region_file): path} for every region GenBank in a run."""
    out = {}
    for gbk in Path(antismash_dir).glob('*/*.region*.gbk'):
        out[(gbk.parent.name, gbk.name)] = gbk
    return out


def gene_and_product(gbk):
    """(CDS count, product) read straight from the GenBank text.

    Counting features by line prefix rather than parsing with Biopython keeps this
    dependency-free and is exact for antiSMASH output, where every feature key starts
    in column 6.
    """
    cds, product = 0, ''
    for line in gbk.read_text().splitlines():
        if line.startswith('     CDS '):
            cds += 1
        elif not product and line.strip().startswith('/product='):
            product = line.split('=', 1)[1].strip().strip('"')
        elif line.startswith('ORIGIN'):
            break
    return cds, product


def families(bigscape_db, cutoff):
    """{(genome, region_file): family_id} from a run's clustering database."""
    if not Path(bigscape_db).is_file():
        return {}
    con = sqlite3.connect(bigscape_db)
    out = {}
    for path, fam in con.execute(
            """SELECT g.path, bf.family_id FROM gbk g
               JOIN bgc_record br ON br.gbk_id = g.id
               JOIN bgc_record_family bf ON bf.record_id = br.id
               JOIN family f ON f.id = bf.family_id
               WHERE br.record_type = 'region' AND f.cutoff = ?""", (cutoff,)):
        parts = Path(path).parts
        out[(parts[-2], parts[-1])] = fam
    return out


def rand_index(pairs_a, pairs_b, keys):
    """Adjusted Rand index over co-membership, the measure the validation matrix uses."""
    keys = sorted(keys)
    n = len(keys)
    if n < 2:
        return None
    contingency = defaultdict(int)
    row, col = defaultdict(int), defaultdict(int)
    for k in keys:
        contingency[(pairs_a[k], pairs_b[k])] += 1
        row[pairs_a[k]] += 1
        col[pairs_b[k]] += 1
    comb2 = lambda x: x * (x - 1) / 2
    sum_ij = sum(comb2(v) for v in contingency.values())
    sum_i = sum(comb2(v) for v in row.values())
    sum_j = sum(comb2(v) for v in col.values())
    expected = sum_i * sum_j / comb2(n)
    maximum = (sum_i + sum_j) / 2
    return None if maximum == expected else (sum_ij - expected) / (maximum - expected)


def read_tsv(path):
    if not Path(path).is_file():
        return []
    rows = [l.rstrip('\n').split('\t') for l in open(path) if l.strip()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', type=Path, default=Path('results'))
    ap.add_argument('--before', required=True, help='run directory name, e.g. Erwiniaceae_pre_recovery')
    ap.add_argument('--after', required=True)
    ap.add_argument('--taxon-db', default=None,
                    help='basename of the BiG-SCAPE db inside each run (default: <after>.db, '
                         'falling back to the only .db present)')
    ap.add_argument('--cutoff', type=float, default=0.3)
    ap.add_argument('--outdir', type=Path, required=True)
    a = ap.parse_args()

    runs = {}
    for label, name in (('before', a.before), ('after', a.after)):
        antismash = a.results / 'antismash_results' / name
        bs_dir = a.results / 'bigscape_results' / name
        dbs = sorted(bs_dir.glob('*.db'))
        runs[label] = {
            'name': name,
            'regions': region_files(antismash),
            'families': families(dbs[0], a.cutoff) if dbs else {},
            'novelty': read_tsv(a.results / 'main_analysis_results' / name / 'novelty_ranking.tsv'),
            'genomes_screened': len(list(antismash.glob('*/'))),
        }

    keys = sorted(set(runs['before']['regions']) | set(runs['after']['regions']))
    a.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for key in keys:
        genome, region = key
        row = {'genome': genome, 'region': region}
        for label in ('before', 'after'):
            gbk = runs[label]['regions'].get(key)
            cds, product = gene_and_product(gbk) if gbk else ('', '')
            row[f'cds_{label}'] = cds
            row[f'product_{label}'] = product
            row[f'family_{label}'] = runs[label]['families'].get(key, '')
        row['cds_delta'] = (row['cds_after'] - row['cds_before']
                            if row['cds_after'] != '' and row['cds_before'] != '' else '')
        rows.append(row)

    cols = ['genome', 'region', 'cds_before', 'cds_after', 'cds_delta',
            'product_before', 'product_after', 'family_before', 'family_after']
    with open(a.outdir / 'per_bgc.tsv', 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # Per-GCF view, from the novelty ranking each run already publishes.
    gcf_cols = ['gcf', 'members', 'genomes', 'coupling_class', 'rank', 'intact']
    with open(a.outdir / 'per_gcf.tsv', 'w') as fh:
        fh.write('run\t' + '\t'.join(gcf_cols) + '\n')
        for label in ('before', 'after'):
            for r in runs[label]['novelty']:
                fh.write(label + '\t' + '\t'.join(str(r.get(c, '')) for c in gcf_cols) + '\n')

    shared = [r for r in rows if r['family_before'] != '' and r['family_after'] != '']
    ari = rand_index({(r['genome'], r['region']): r['family_before'] for r in shared},
                     {(r['genome'], r['region']): r['family_after'] for r in shared},
                     [(r['genome'], r['region']) for r in shared])
    gained = [r for r in rows if r['cds_delta'] != '' and r['cds_delta'] > 0]
    summary = {
        'before': runs['before']['name'],
        'after': runs['after']['name'],
        'cutoff': a.cutoff,
        'genomes_with_antismash_results': {k: runs[k]['genomes_screened'] for k in runs},
        'bgcs': {k: len(runs[k]['regions']) for k in runs},
        'bgcs_only_in_before': [f'{g}/{r}' for g, r in
                                sorted(set(runs['before']['regions']) - set(runs['after']['regions']))],
        'bgcs_only_in_after': [f'{g}/{r}' for g, r in
                               sorted(set(runs['after']['regions']) - set(runs['before']['regions']))],
        'families': {k: len(set(runs[k]['families'].values())) for k in runs},
        'gcf_membership_ari': ari,
        'bgcs_compared': len(shared),
        'genes': {
            'bgcs_that_gained_genes': len(gained),
            'total_genes_gained': sum(r['cds_delta'] for r in gained),
            'largest_gain': max((r['cds_delta'] for r in gained), default=0),
            'cds_total_before': sum(r['cds_before'] for r in rows if r['cds_before'] != ''),
            'cds_total_after': sum(r['cds_after'] for r in rows if r['cds_after'] != ''),
        },
        'product_changes': sorted({f"{r['product_before']} -> {r['product_after']}"
                                   for r in rows
                                   if r['product_before'] and r['product_after']
                                   and r['product_before'] != r['product_after']}),
    }
    (a.outdir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
