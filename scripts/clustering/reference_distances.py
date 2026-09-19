#!/usr/bin/env python3
"""Distance from every BGC to the characterised reference clusters.

Reads the database produced by the reference pass (a copy of the finished clustering
database, with the references added) and joins each measured pair back to the GCF the
main run assigned — so the answer is phrased in the families the report already shows,
and the published clustering database never has to contain a reference.

References are identified by content hash, not by path. BiG-SCAPE deduplicates input on
the sha256 of the file (read as text, so a CRLF file does not hash as its raw bytes),
and when a reference is byte-identical to one of the query BGCs it silently drops the
reference and keeps the query — logged at INFO and easy to miss. Hashing finds the
reference either way, and that case is reported as an exact match rather than lost.

Emits:
  --distances  one row per reference x BGC pair, with the BiG-SCAPE components
  --summary    per-GCF nearest reference, and per-reference nearest GCF
"""
import argparse
import hashlib
import json
import os
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path


def bigscape_hash(path):
    """The content hash BiG-SCAPE stores in gbk.hash (genbank/gbk.py)."""
    with open(path, 'r') as fh:
        return hashlib.sha256(fh.read().encode('utf-8')).hexdigest()


def genome_of(gbk_path):
    """Genome name from a staged region path, as the other DB readers derive it."""
    parts = Path(gbk_path).parts
    for i, part in enumerate(parts):
        if part == 'antismash_input' and i + 1 < len(parts):
            return parts[i + 1]
    return parts[-2] if len(parts) > 1 else '?'


def region_of(gbk_path):
    return os.path.basename(gbk_path).replace('.gbk', '')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pass-db', type=Path, required=True,
                    help='database from the reference pass (queries + references)')
    ap.add_argument('--main-db', type=Path, required=True,
                    help='the run\'s own clustering database, for GCF assignments')
    ap.add_argument('--reference-dir', type=Path, required=True)
    ap.add_argument('--cutoff', type=float, default=0.3)
    ap.add_argument('--distances', type=Path, required=True)
    ap.add_argument('--summary', type=Path, required=True)
    a = ap.parse_args()

    refs = {}   # content hash -> reference name
    for gbk in sorted(a.reference_dir.glob('*.gbk')):
        refs[bigscape_hash(gbk)] = gbk.name.replace('.region001.gbk', '').replace('.gbk', '')
    if not refs:
        sys.exit(f'{a.reference_dir}: no .gbk references')

    con = sqlite3.connect(a.pass_db)
    con.row_factory = sqlite3.Row

    # gbk id -> reference name, for the references present in this database
    ref_gbk, dropped = {}, {}
    for row in con.execute('SELECT id, path, hash FROM gbk'):
        name = refs.get(row['hash'])
        if name is None:
            continue
        ref_gbk[row['id']] = name
        # A reference whose path is not in the reference directory is one BiG-SCAPE
        # deduplicated against an identical query: the query record IS that cluster.
        if Path(row['path']).parent.resolve() != a.reference_dir.resolve():
            dropped[name] = row['path']
    missing = sorted(set(refs.values()) - set(ref_gbk.values()))

    # The main run's family per region record. The pass database is a copy, so record
    # ids are the same on both sides and can be joined directly.
    main = sqlite3.connect(a.main_db)
    family = {}
    for rec_id, fam in main.execute(
            """SELECT brf.record_id, brf.family_id FROM bgc_record_family brf
               JOIN family f ON f.id = brf.family_id
               WHERE f.cutoff = ?""", (a.cutoff,)):
        family[rec_id] = fam
    family_size = defaultdict(int)
    for fam in family.values():
        family_size[fam] += 1

    rows = []
    for r in con.execute(
            """SELECT d.record_a_id AS a_id, d.record_b_id AS b_id, d.distance, d.jaccard,
                      d.adjacency, d.dss, ga.id AS ga_id, gb.id AS gb_id,
                      ga.path AS a_path, gb.path AS b_path
               FROM distance d
               JOIN bgc_record ra ON ra.id = d.record_a_id JOIN gbk ga ON ga.id = ra.gbk_id
               JOIN bgc_record rb ON rb.id = d.record_b_id JOIN gbk gb ON gb.id = rb.gbk_id
               WHERE ra.record_type = 'region' AND rb.record_type = 'region'"""):
        a_ref, b_ref = r['ga_id'] in ref_gbk, r['gb_id'] in ref_gbk
        if a_ref == b_ref:
            continue  # query-query, or reference-reference: neither is asked for here
        ref_name = ref_gbk[r['ga_id'] if a_ref else r['gb_id']]
        q_path = r['b_path'] if a_ref else r['a_path']
        q_id = r['b_id'] if a_ref else r['a_id']
        rows.append({
            'reference': ref_name,
            'genome': genome_of(q_path),
            'region': region_of(q_path),
            'family_id': family.get(q_id, ''),
            'distance': round(r['distance'], 4),
            'jaccard': round(r['jaccard'], 4),
            'adjacency': round(r['adjacency'], 4),
            'dss': round(r['dss'], 4),
        })

    # A deduplicated reference has no pairs of its own: the query it matched carries
    # them. Record the identity so the summary can still name the family it belongs to.
    for name, path in dropped.items():
        rec = con.execute("SELECT br.id FROM bgc_record br JOIN gbk g ON g.id = br.gbk_id "
                          "WHERE g.path = ? AND br.record_type = 'region'",
                          (path,)).fetchone()
        rows.append({
            'reference': name, 'genome': genome_of(path), 'region': region_of(path),
            'family_id': family.get(rec['id'], '') if rec else '',
            'distance': 0.0, 'jaccard': 1.0, 'adjacency': 1.0, 'dss': 1.0,
        })

    rows.sort(key=lambda r: (r['reference'], r['distance']))
    cols = ['reference', 'genome', 'region', 'family_id', 'distance', 'jaccard',
            'adjacency', 'dss']
    with open(a.distances, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # Per GCF: its nearest reference, and how many members are within the cutoff.
    per_family = defaultdict(lambda: {'nearest': None, 'distance': None, 'within_cutoff': {}})
    for r in rows:
        fam = r['family_id']
        if fam == '':
            continue
        entry = per_family[fam]
        if entry['distance'] is None or r['distance'] < entry['distance']:
            entry['nearest'], entry['distance'] = r['reference'], r['distance']
        if r['distance'] <= a.cutoff:
            entry['within_cutoff'][r['reference']] = \
                entry['within_cutoff'].get(r['reference'], 0) + 1

    per_reference = {}
    for name in sorted(set(refs.values())):
        mine = [r for r in rows if r['reference'] == name]
        best = min(mine, key=lambda r: r['distance']) if mine else None
        per_reference[name] = {
            'loaded': name in ref_gbk.values(),
            'identical_to_query': dropped.get(name),
            'nearest_distance': best['distance'] if best else None,
            'nearest_genome': best['genome'] if best else None,
            'nearest_family': best['family_id'] if best else None,
            'bgcs_within_cutoff': sum(1 for r in mine if r['distance'] <= a.cutoff),
        }

    summary = {
        'cutoff': a.cutoff,
        'references': per_reference,
        'references_not_loaded': missing,
        'families': {str(f): {
            'nearest_reference': v['nearest'],
            'distance': v['distance'],
            'members': family_size.get(f, 0),
            'members_within_cutoff': v['within_cutoff'],
        } for f, v in sorted(per_family.items(), key=lambda kv: str(kv[0]))},
    }
    a.summary.write_text(json.dumps(summary, indent=2) + '\n')

    print(f'{len(rows)} reference-BGC pairs over {len(ref_gbk)} loaded references')
    if dropped:
        print(f'identical to a query BGC (BiG-SCAPE dropped the file): '
              f'{", ".join(sorted(dropped))}')
    if missing:
        # Silent rejection is this directory's characteristic failure: a file without
        # an antiSMASH region feature never appears among the loaded records.
        print(f'WARNING: {len(missing)} reference(s) never loaded: {", ".join(missing)}',
              file=sys.stderr)


if __name__ == '__main__':
    main()
