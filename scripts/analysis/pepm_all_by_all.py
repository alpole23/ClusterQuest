#!/usr/bin/env python3
"""All-by-all pepM comparison, to test whether pepM identity can partition BiG-SCAPE.

BiG-SCAPE compares every BGC against every other, so its memory grows
quadratically and that is the constraint that stops a large run clustering in
one pass. If pepM identity tracked gene-neighbourhood similarity closely enough,
it would be a cheap key to split the input on. This measures, for a given run,
how far the input could be split and whether doing so would separate BGCs that
BiG-SCAPE placed in one family.

Two things make it cheap rather than a second clustering run.

**BiG-SCAPE already computed the neighbourhood axis.** Its `distance` table holds
every pair -- exactly n(n-1)/2 rows, verified -- so nothing is recompared here;
the pairs are joined on their record ids.

**Identity comes from a profile alignment, not BLAST.** `hmmalign` against
PF13714 with pairwise deletion of gapped sites is linear in sequence count where
all-by-all alignment is quadratic, which is what makes this tractable at scales
where BiG-SCAPE itself is not. Only match columns count, so a fusion protein's
extra residues fall out as insertions rather than dragging identity down.

This used to also plot pepM identity against neighbourhood similarity, as a
replication of Yu et al. PNAS 2013;110(51):20759 Fig. 2B. That figure was
removed because on a taxonomically broad run it mostly measures taxonomy.
Splitting the Enterobacterales pairs by GTDB relationship, the pooled r of 0.789
over the paper's 0.6-1.0 window becomes 0.353 between genera, 0.194 within a
genus and 0.385 within a species -- 92.4% of the covariance in that window is
BETWEEN those strata, not within them. The three groups are near-disjoint clouds
on both axes (identity 0.725 / 0.967 / 0.994, similarity 0.20 / 0.63 / 0.87), so
the regression is a line through three points. Within a species, where pepM
spans 0.991-1.000, neighbourhood similarity still has sd 0.225: pepM says
nothing there. Over the full identity range the confound is weaker (63.4%
between strata) but the fitted window is what the figure reported.

None of that touches the partitioning question below, which asks only whether a
pepM cut separates same-family pairs -- an operational claim about where the
input can be split, not a claim that identity predicts chemistry.

    python scripts/analysis/pepm_all_by_all.py \\
        --db results/bigscape_results/Erwiniaceae/Erwiniaceae.db \\
        --pfam results/databases/pfam38.2/Pfam-A.hmm \\
        --outdir results/main_analysis_results/Erwiniaceae/pepm_all_by_all
"""
import argparse
import collections
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

PEPM_ACCESSION = 'PF13714'   # PEP_mutase — the phosphonate detection rule's hallmark
GAP = 0                      # sentinel for a gap/deleted site in the encoded alignment


# ---------------------------------------------------------------- extraction

def fetch_pepm_sequences(con, accession, organism_prefix=None):
    """One pepM protein per region record: the highest-scoring PF13714 CDS.

    A handful of regions carry two copies (2 of 333 in Erwiniaceae), so the hit
    is chosen by bit score rather than assuming uniqueness. Region records map
    1:1 onto gbk files for antiSMASH region output, which is why the CDS join
    goes through gbk_id.
    """
    # --organism lets one run be sliced by taxon without re-clustering. The
    # `distance` values are pairwise comparisons of two BGCs' domain content, so
    # they do not depend on which other genomes were in the run; only family
    # assignment does. Slicing therefore gives an honest pepM-vs-neighbourhood
    # relationship for the subset, but its GCF labels are the parent run's.
    sql = """
        SELECT br.id AS record_id, cds.aa_seq AS seq, hsp.bit_score AS score
        FROM bgc_record br
        JOIN cds ON cds.gbk_id = br.gbk_id
        JOIN hsp ON hsp.cds_id = cds.id
        JOIN gbk g ON g.id = br.gbk_id
        WHERE br.record_type = 'region' AND hsp.accession LIKE ?
    """
    args = [f'{accession}%']
    if organism_prefix:
        sql += ' AND g.organism LIKE ?'
        args.append(f'{organism_prefix}%')
    sql += ' ORDER BY br.id, hsp.bit_score DESC'
    rows = con.execute(sql, args).fetchall()
    best = {}
    for r in rows:
        if r['record_id'] not in best and r['seq']:
            best[r['record_id']] = r['seq'].replace('*', '')
    return best


def extract_profile(pfam_hmm, accession, dest, hmmfetch):
    """Pull one profile out of Pfam-A.hmm.

    hmmfetch needs an .ssi index that hmmpress does not create, and the Pfam
    copy may sit on a read-only store, so the index is built beside a symlink
    in a writable directory rather than next to the original.
    """
    work = dest.parent / '_pfam'
    work.mkdir(parents=True, exist_ok=True)
    link = work / 'Pfam-A.hmm'
    if not link.exists():
        os.symlink(Path(pfam_hmm).resolve(), link)
    if not (work / 'Pfam-A.hmm.ssi').exists():
        subprocess.run([hmmfetch, '--index', str(link)], check=True,
                       capture_output=True, text=True)

    # Pfam stores versioned accessions (ACC PF13714.13) and hmmfetch's index keys
    # on the exact string, so a bare PF13714 misses. Resolve the version from the
    # file rather than hardcoding it — params.pfam_release is a pin that will be
    # bumped, and the accession version moves with it.
    key = None
    with link.open() as fh:
        for line in fh:
            if line.startswith('ACC ') and line.split()[1].split('.')[0] == accession:
                key = line.split()[1]
                break
    if key is None:
        raise SystemExit(f'{accession} not found in {pfam_hmm}')

    with dest.open('w') as fh:
        subprocess.run([hmmfetch, str(link), key], check=True,
                       stdout=fh, stderr=subprocess.PIPE, text=True)
    return dest


def align(seqs, profile, outdir, hmmalign):
    """hmmalign the pepMs to the profile; return (ids, match-column matrix).

    Only match columns are kept — uppercase residues and '-' deletions in
    Stockholm. Insert columns (lowercase, '.') are per-sequence excursions from
    the profile and are not homologous across sequences, so counting them would
    penalise fusion proteins for residues nobody else was aligned against.
    """
    fasta = outdir / 'pepm.faa'
    with fasta.open('w') as fh:
        for rid, seq in sorted(seqs.items()):
            fh.write(f'>{rid}\n{seq}\n')

    sto = outdir / 'pepm.sto'
    with sto.open('w') as fh:
        subprocess.run([hmmalign, '--amino', '--trim', str(profile), str(fasta)],
                       check=True, stdout=fh, stderr=subprocess.PIPE, text=True)

    aln = {}
    for line in sto.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line == '//':
            continue
        name, _, chunk = line.partition(' ')
        aln[name] = aln.get(name, '') + chunk.strip()

    ids = sorted(aln, key=int)
    rows = []
    for i in ids:
        # keep match columns only: uppercase residue, or '-' meaning deleted
        rows.append([ord(ch) if ch.isupper() else (GAP if ch == '-' else None)
                     for ch in aln[i]])
    keep = [k for k in range(len(rows[0])) if rows[0][k] is not None]
    mat = np.zeros((len(ids), len(keep)), dtype=np.uint8)
    for r, row in enumerate(rows):
        mat[r] = [row[k] if row[k] is not None else GAP for k in keep]
    return ids, mat


# ---------------------------------------------------------------- comparison

def pairwise_identity(ids, mat):
    """Every pair's identity, deleting sites missing from either sequence.

    Vectorised one row against all later rows: 55k pairs over ~300 columns is
    milliseconds, where a Python-level double loop is minutes.
    """
    n = len(ids)
    out = {}
    for i in range(n - 1):
        a = mat[i]
        rest = mat[i + 1:]
        both = (a != GAP) & (rest != GAP)          # pairwise deletion
        compared = both.sum(axis=1)
        matches = (both & (a == rest)).sum(axis=1)
        with np.errstate(invalid='ignore', divide='ignore'):
            ident = np.where(compared > 0, matches / np.maximum(compared, 1), np.nan)
        for k, j in enumerate(range(i + 1, n)):
            if compared[k]:
                out[(ids[i], ids[j])] = (float(ident[k]), int(compared[k]))
    return out


def join_neighbourhood(con, identity):
    """Attach BiG-SCAPE's similarity to each pepM pair.

    `jaccard` is shared domain content, the analogue of the paper's fraction of
    homologous genes shared; `1 - distance` is BiG-SCAPE's composite similarity,
    kept alongside because it is what GCF membership is actually cut on.
    """
    rows = []
    for (a, b), (ident, ncols) in identity.items():
        r = con.execute(
            """SELECT distance, jaccard FROM distance
               WHERE (record_a_id = ? AND record_b_id = ?)
                  OR (record_a_id = ? AND record_b_id = ?) LIMIT 1""",
            (a, b, b, a),
        ).fetchone()
        if r is not None:
            rows.append((a, b, ident, ncols, r['jaccard'], 1.0 - r['distance']))
    return rows


def partition_analysis(rows, thresholds, gcf_similarity_cut=0.70):
    """Could pepM identity partition BiG-SCAPE's all-pairs problem?

    Two questions, and the second is the one that decides it.

    *Is a cut lossless?* Every pair BiG-SCAPE puts in one GCF must survive the
    cut, or partitioning would split real families. `lost` counts pairs above
    the GCF similarity cutoff that the threshold would have separated.

    **`lost` is a lower bound, not a count.** It uses pairwise similarity as a
    proxy for family membership, but BiG-SCAPE families are transitively closed
    clusters: 9% of same-family pairs in Streptomyces sit *below* the 0.70
    cutoff, joined through a third BGC rather than directly. Measured against
    an actual partitioned re-run at 0.90, this predicted 7 splits where 23
    occurred (26 same-family pairs were separable). Zero remains a reliable
    all-clear — nothing separable means nothing splits — but a non-zero value
    understates the damage, so treat any non-zero as disqualifying rather than
    as a budget. `scripts/bench_bigscape_partitioned.py` measures the truth.

    *Does the cut actually divide anything?* Single-linkage components at the
    threshold are the jobs BiG-SCAPE would then run. Since its cost is
    quadratic, what matters is the largest component's share of the whole:
    work falls as that fraction squared, so a threshold that leaves one
    dominant component buys almost nothing however lossless it is.
    """
    nodes = {r[0] for r in rows} | {r[1] for r in rows}
    out = []
    for t in thresholds:
        parent = {n: n for n in nodes}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        lost = 0
        for a, b, ident, _, _, sim in rows:
            if ident >= t:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
            elif sim >= gcf_similarity_cut:
                lost += 1

        sizes = sorted(collections.Counter(find(n) for n in nodes).values(), reverse=True)
        frac = sizes[0] / len(nodes)
        out.append({
            'threshold': t, 'components': len(sizes), 'largest': sizes[0],
            'largest_fraction': frac, 'same_gcf_pairs_lost': lost,
            'lossless': lost == 0,
            # BiG-SCAPE is O(n^2), so the surviving work is the sum of squared shares
            'relative_work': sum((s / len(nodes)) ** 2 for s in sizes),
        })
    return out


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--db', type=Path, required=True, help='BiG-SCAPE SQLite database')
    ap.add_argument('--pfam', type=Path, required=True, help='Pfam-A.hmm')
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--accession', default=PEPM_ACCESSION)
    ap.add_argument('--organism', default=None,
                    help='restrict to organisms with this prefix, e.g. "Pantoea"')
    ap.add_argument('--hmmfetch', default=shutil.which('hmmfetch') or 'hmmfetch')
    ap.add_argument('--hmmalign', default=shutil.which('hmmalign') or 'hmmalign')
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f'file:{args.db}?mode=ro', uri=True)
    con.row_factory = sqlite3.Row

    seqs = fetch_pepm_sequences(con, args.accession, args.organism)
    n_regions = con.execute(
        "SELECT COUNT(*) FROM bgc_record WHERE record_type='region'").fetchone()[0]
    if args.organism:
        print(f'restricted to organism prefix {args.organism!r}')
    print(f'pepM sequences: {len(seqs)} of {n_regions} regions '
          f'({100 * len(seqs) / n_regions:.1f}%)')
    if len(seqs) < 3:
        # Not an error: a taxon can legitimately carry one or two phosphonate
        # BGCs (S. hygroscopicus yields 2 across 39 genomes), and a correlation
        # over fewer than three points is undefined rather than wrong. Failing
        # here took down an otherwise complete run AFTER all detection and
        # clustering had finished, which is the worst possible moment.
        print(f'only {len(seqs)} pepM sequences; a pairwise correlation needs '
              f'at least 3. Skipping this analysis, the run is unaffected.',
              file=sys.stderr)
        args.outdir.mkdir(parents=True, exist_ok=True)
        (args.outdir / 'pepm_all_by_all.json').write_text(json.dumps({
            'skipped': True,
            'reason': f'only {len(seqs)} pepM sequences of {n_regions} regions',
            'n_pepm': len(seqs),
            'n_regions': n_regions,
        }, indent=2) + '\n')
        return 0

    profile = extract_profile(args.pfam, args.accession,
                              args.outdir / f'{args.accession}.hmm', args.hmmfetch)
    ids, mat = align(seqs, profile, args.outdir, args.hmmalign)
    print(f'alignment: {len(ids)} sequences x {mat.shape[1]} match columns')

    identity = pairwise_identity(ids, mat)
    rows = join_neighbourhood(con, identity)
    expected = len(ids) * (len(ids) - 1) // 2
    print(f'pairs: {len(identity):,} compared, {len(rows):,} joined to BiG-SCAPE '
          f'(all-pairs would be {expected:,})')

    with (args.outdir / 'pepm_vs_neighbourhood.tsv').open('w') as fh:
        fh.write('record_a\trecord_b\tpepm_identity\taligned_sites\t'
                 'jaccard\tbigscape_similarity\n')
        for r in rows:
            fh.write(f'{r[0]}\t{r[1]}\t{r[2]:.6f}\t{r[3]}\t{r[4]:.6f}\t{r[5]:.6f}\n')

    summary = {'accession': args.accession, 'organism': args.organism,
               'sequences': len(ids),
               'match_columns': int(mat.shape[1]), 'pairs': len(rows)}

    parts = partition_analysis(rows, [0.5, 0.6, 0.7, 0.8, 0.9])
    summary['partitioning'] = parts
    print('\n  partitioning BiG-SCAPE by pepM identity:')
    print(f"  {'cut':>6}{'components':>12}{'largest':>9}{'share':>8}"
          f"{'same-GCF lost':>15}{'work vs one job':>17}")
    for p_ in parts:
        print(f"  {p_['threshold']:>6.2f}{p_['components']:>12}{p_['largest']:>9}"
              f"{p_['largest_fraction']:>7.0%}{p_['same_gcf_pairs_lost']:>15,}"
              f"{p_['relative_work']:>16.0%}")

    (args.outdir / 'pepm_all_by_all.json').write_text(json.dumps(summary, indent=2))
    print(f'wrote {args.outdir}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
