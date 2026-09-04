#!/usr/bin/env python3
"""Merge per-partition BiG-SCAPE databases into one.

Partitioned clustering produces one SQLite database per partition, but every
downstream consumer — clustering stats, GCF representatives, the biosynthetic
tree, rarefaction, the report — reads a single `{taxon}.db`. This reassembles
them so partitioning is invisible above this line.

**Every id is a per-partition autoincrement**, so they collide. Each table is
copied with its primary key offset into a fresh range and every foreign key
rewritten through the same mapping. The table order below is a topological sort:
a table is only copied once everything it points at has been mapped.

**`distance` is deliberately incomplete after merging.** It holds only
within-partition comparisons, which is the entire point — the cross-partition
pairs are the ones pepM identity established cannot cluster together
(ARI 1.0000 against unpartitioned runs on three sets). Anything reading
`distance` must therefore treat a missing pair as "not compared", not as
"distance zero". `analysis/pepm_all_by_all.py` already does, by inner-joining.

**`run` and `edge_params` are deduplicated rather than offset.** Every partition
ran with identical parameters, so keeping one row each preserves meaning;
offsetting would invent parameter sets that never existed.

    python scripts/clustering/merge_bigscape_dbs.py \\
        --inputs part_000.db part_001.db ... --out merged.db
"""
import argparse
import sqlite3
import sys
from pathlib import Path

# (table, primary key or None, {column: table it references})
# Ordered so a table's referents are always mapped before it is copied.
SPEC = [
    ('gbk',                 'id',  {}),
    ('bgc_record',          'id',  {'gbk_id': 'gbk', 'parent_id': 'bgc_record'}),
    ('cds',                 'id',  {'gbk_id': 'gbk'}),
    ('hsp',                 'id',  {'cds_id': 'cds'}),
    ('hsp_alignment',       None,  {'hsp_id': 'hsp'}),
    ('scanned_cds',         None,  {'cds_id': 'cds'}),
    ('family',              'id',  {'center_id': 'bgc_record', 'run_id': 'run'}),
    ('bgc_record_family',   None,  {'record_id': 'bgc_record', 'family_id': 'family'}),
    ('connected_component', None,  {'record_id': 'bgc_record', 'run_id': 'run'}),
    ('distance',            None,  {'record_a_id': 'bgc_record',
                                    'record_b_id': 'bgc_record',
                                    'edge_param_id': 'edge_params'}),
]
SHARED = ('run', 'edge_params')   # identical across partitions; keep one row each


def columns(con, table):
    return [r[1] for r in con.execute(f'PRAGMA table_info({table})')]


def copy_schema(src, dest):
    """Recreate tables, indexes and triggers exactly as BiG-SCAPE built them."""
    for (sql,) in src.execute(
            "SELECT sql FROM sqlite_master WHERE sql IS NOT NULL "
            "AND name NOT LIKE 'sqlite_%'"):
        dest.execute(sql)
    dest.commit()


def merge(inputs, out):
    out_path = Path(out)
    if out_path.exists():
        out_path.unlink()
    dest = sqlite3.connect(out_path)
    first = sqlite3.connect(f'file:{inputs[0]}?mode=ro', uri=True)
    copy_schema(first, dest)

    present = {r[0] for r in dest.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}

    # The shared tables come from the first partition and stand for all of them.
    shared_map = {}
    for table in SHARED:
        if table not in present:
            continue
        cols = columns(first, table)
        rows = first.execute(f'SELECT * FROM {table}').fetchall()
        if not rows:
            continue
        dest.execute(f"INSERT INTO {table} ({','.join(cols)}) "
                     f"VALUES ({','.join('?' * len(cols))})", rows[0])
        keep = rows[0][cols.index('id')] if 'id' in cols else None
        shared_map[table] = keep
    first.close()
    dest.commit()

    stats = {t: 0 for t, _, _ in SPEC}
    for part, db in enumerate(inputs):
        src = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
        # idmap[table][old id] -> new id, built as each table is copied
        idmap = {t: {} for t, _, _ in SPEC}
        for table in SHARED:
            idmap[table] = {}

        for table, pk, refs in SPEC:
            if table not in present:
                continue
            cols = columns(src, table)
            rows = src.execute(f'SELECT * FROM {table}').fetchall()
            if not rows:
                continue
            base = dest.execute(
                f'SELECT COALESCE(MAX({pk}), 0) FROM {table}').fetchone()[0] if pk else 0

            out_rows = []
            for i, row in enumerate(rows):
                r = list(row)
                if pk:
                    old = r[cols.index(pk)]
                    new = base + i + 1
                    idmap[table][old] = new
                    r[cols.index(pk)] = new
                for col, target in refs.items():
                    if col not in cols:
                        continue
                    j = cols.index(col)
                    if r[j] is None:
                        continue
                    if target in SHARED:
                        r[j] = shared_map.get(target, r[j])
                    else:
                        # A self-reference (bgc_record.parent_id) resolves within
                        # the same table, which is already mapped by this point
                        # because rows are numbered before the column is rewritten.
                        r[j] = idmap[target].get(r[j], r[j])
                out_rows.append(r)

            dest.executemany(
                f"INSERT INTO {table} ({','.join(cols)}) "
                f"VALUES ({','.join('?' * len(cols))})", out_rows)
            stats[table] += len(out_rows)
        src.close()
        dest.commit()
        print(f'  partition {part:>3}: merged {Path(db).name}')

    # sqlite_sequence drives AUTOINCREMENT; leaving it stale would let a later
    # writer reuse an id we just assigned.
    for table, pk, _ in SPEC:
        if pk and table in present:
            hi = dest.execute(f'SELECT COALESCE(MAX({pk}), 0) FROM {table}').fetchone()[0]
            dest.execute('INSERT OR REPLACE INTO sqlite_sequence(name, seq) VALUES (?, ?)',
                         (table, hi))
    dest.commit()
    return dest, stats


def verify(dest, inputs):
    """Referential integrity, plus the counts a caller would want to see."""
    problems = []
    checks = [
        ('bgc_record', 'gbk_id', 'gbk', 'id'),
        ('cds', 'gbk_id', 'gbk', 'id'),
        ('hsp', 'cds_id', 'cds', 'id'),
        ('bgc_record_family', 'record_id', 'bgc_record', 'id'),
        ('bgc_record_family', 'family_id', 'family', 'id'),
        ('family', 'center_id', 'bgc_record', 'id'),
        ('distance', 'record_a_id', 'bgc_record', 'id'),
        ('distance', 'record_b_id', 'bgc_record', 'id'),
    ]
    for child, col, parent, pcol in checks:
        n = dest.execute(
            f'SELECT COUNT(*) FROM {child} c LEFT JOIN {parent} p ON c.{col} = p.{pcol} '
            f'WHERE c.{col} IS NOT NULL AND p.{pcol} IS NULL').fetchone()[0]
        if n:
            problems.append(f'{child}.{col} has {n} rows pointing nowhere')

    dupes = dest.execute(
        'SELECT COUNT(*) FROM (SELECT record_id FROM bgc_record_family '
        'GROUP BY record_id HAVING COUNT(*) > 1)').fetchone()[0]
    if dupes:
        problems.append(f'{dupes} BGC records assigned to more than one family')
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs', nargs='+', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--expect_regions', type=int, default=None,
                    help='fail unless the merged database holds exactly this many '
                         'region records; pass the partitioner\'s BGC count')
    ap.add_argument('--partitions', type=Path, default=None,
                    help="the partitioner's TSV; derives --expect_regions from it")
    args = ap.parse_args()

    # A missing partition database is invisible without this: the merge succeeds,
    # referential integrity holds, and the result is simply short some BGCs.
    # Caught in testing when 14 of 19 partition databases merged cleanly to 179
    # regions where the reference had 185.
    expect = args.expect_regions
    if expect is None and args.partitions:
        with args.partitions.open() as fh:
            next(fh, None)
            expect = sum(1 for line in fh if line.strip())

    inputs = sorted(args.inputs)
    print(f'merging {len(inputs)} partition databases')
    dest, stats = merge(inputs, args.out)

    regions = dest.execute(
        "SELECT COUNT(*) FROM bgc_record WHERE record_type='region'").fetchone()[0]
    fams = dest.execute('SELECT COUNT(DISTINCT id) FROM family').fetchone()[0]
    print(f'\n  {regions} regions, {fams} families, '
          f'{stats["distance"]:,} within-partition comparisons')
    print(f'  (all-pairs would be {regions * (regions - 1) // 2:,}; the difference is '
          f'the cross-partition pairs pepM identity ruled out)')

    problems = verify(dest, inputs)
    if expect is not None and regions != expect:
        problems.append(
            f'expected {expect} regions but merged {regions} — '
            f'{expect - regions} lost, most likely a partition database that was '
            f'never produced')
    dest.close()
    if problems:
        for p in problems:
            print(f'  ERROR: {p}', file=sys.stderr)
        return 1
    print('  referential integrity OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
