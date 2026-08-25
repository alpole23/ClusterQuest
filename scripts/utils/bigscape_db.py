"""Shared read access to the BiG-SCAPE SQLite database.

The BGC tree scripts (bgc_pfam_tree, bgc_synteny_tree, bgc_architecture_tree,
bgc_all_bgcs_tree, bgc_gcf_tree) all pull the same records out of the same tables;
the queries live here so the schema is described in exactly one place.

Schema notes:
  bgc_record          one row per BGC record (region / cand_cluster / protocluster)
  gbk                 the source GenBank file (path, organism)
  cds / scanned_cds   coding sequences, with genomic coordinates
  hsp                 Pfam domain hits per CDS (accession + bit_score)
  family / bgc_record_family   GCF assignments per cutoff
"""

import os
import sqlite3

# Minimum Pfam hit bit score. antiSMASH/BiG-SCAPE report weak hits too; 20 is the
# threshold every tree script has used.
MIN_BIT_SCORE = 20

# Default GCF distance cutoff (matches the pipeline's bigscape_cutoffs default)
DEFAULT_CUTOFF = 0.3

_BGC_SELECT = """
    SELECT
        br.id       AS bgc_id,
        br.product,
        g.id        AS gbk_id,
        g.path      AS gbk_path,
        g.organism
    FROM bgc_record br
    JOIN gbk g ON br.gbk_id = g.id
"""

BGC_QUERY = _BGC_SELECT + """
    WHERE LOWER(br.product) LIKE ?
    ORDER BY g.path
"""

BGC_QUERY_FAMILY = _BGC_SELECT + """
    JOIN bgc_record_family brf ON br.id = brf.record_id
    WHERE LOWER(br.product) LIKE ?
      AND brf.family_id = ?
    ORDER BY g.path
"""

FAMILY_QUERY = """
    SELECT brf.family_id, f.cutoff
    FROM bgc_record_family brf
    JOIN family f ON brf.family_id = f.id
    WHERE brf.record_id = ?
    ORDER BY f.cutoff
"""

# All domain hits for a BGC, ordered by genomic position then descending score,
# so the first row seen for a CDS is its best hit.
CDS_DOMAIN_QUERY = f"""
    SELECT
        c.id        AS cds_id,
        c.nt_start,
        c.strand,
        h.accession,
        h.bit_score
    FROM cds c
    JOIN scanned_cds sc ON sc.cds_id = c.id
    JOIN hsp h ON h.cds_id = c.id
    WHERE c.gbk_id = ?
      AND h.accession != ''
      AND h.bit_score >= {MIN_BIT_SCORE}
    ORDER BY c.nt_start ASC, h.bit_score DESC
"""

# Every distinct domain in a BGC (not just the best hit per CDS)
DOMAIN_SET_QUERY = f"""
    SELECT DISTINCT h.accession
    FROM cds c
    JOIN scanned_cds sc ON sc.cds_id = c.id
    JOIN hsp h ON h.cds_id = c.id
    WHERE c.gbk_id = ?
      AND h.accession != ''
      AND h.bit_score >= {MIN_BIT_SCORE}
"""


def connect(db_path):
    """Open the BiG-SCAPE DB with row access by column name."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def fetch_bgc_records(cur, bgc_type_filter, family_id=None):
    """All BGC records whose product matches bgc_type_filter (substring, case-insensitive).

    Restricts to one GCF family when family_id is given.
    """
    pattern = f'%{bgc_type_filter.lower()}%'
    if family_id is not None:
        cur.execute(BGC_QUERY_FAMILY, (pattern, family_id))
    else:
        cur.execute(BGC_QUERY, (pattern,))
    return cur.fetchall()


def fetch_families(cur, bgc_id):
    """GCF assignments for one BGC record → [{'family_id': int, 'cutoff': float}]."""
    cur.execute(FAMILY_QUERY, (bgc_id,))
    return [{'family_id': r['family_id'], 'cutoff': r['cutoff']} for r in cur.fetchall()]


def fetch_domain_set(cur, gbk_id):
    """Set of every Pfam accession annotated in a BGC."""
    cur.execute(DOMAIN_SET_QUERY, (gbk_id,))
    return {r[0] for r in cur.fetchall()}


def fetch_best_domain_per_cds(cur, gbk_id, strip_version=True):
    """Best-scoring Pfam accession per CDS, ordered by genomic position.

    This is the "gene order" view of a BGC: one domain per gene, left to right.
    """
    cur.execute(CDS_DOMAIN_QUERY, (gbk_id,))

    seen_cds = set()
    ordered = []
    for row in cur.fetchall():
        if row['cds_id'] in seen_cds:
            continue
        seen_cds.add(row['cds_id'])
        accession = row['accession']
        ordered.append(accession.split('.')[0] if strip_version else accession)
    return ordered


def record_metadata(row, label):
    """Metadata fields shared by every tree script's per-BGC record."""
    return {
        'label':    label,
        'product':  row['product'],
        'organism': row['organism'] or os.path.basename(os.path.dirname(row['gbk_path'])),
        'gbk_path': row['gbk_path'],
    }
