#!/usr/bin/env python3
"""Measure how BiG-SCAPE runtime and memory grow with BGC count.

The 1M-genome projection in the benchmark report rests on BiG-SCAPE being
all-pairs O(n^2). That is true of the algorithm, but it is not something the
pipeline runs could show: they produced 320 and 333 BGCs, far too close together
to fit a curve. This measures the curve directly.

Getting more BGCs without re-running antiSMASH means replicating the real ones.
BiG-SCAPE deduplicates on the sha256 of the raw file bytes (GBK.__hash__ over
GBK.hash, set in gbk.py), so each copy is rewritten with a unique LOCUS,
ACCESSION and VERSION and survives as a distinct BGC.

What that does and does not measure honestly:

  measured faithfully   all-pairs distance computation over real domain content,
                        which is where the O(n^2) lives, and the fixed cost of
                        loading Pfam and building the database
  distorted             family assignment sees many near-identical BGCs, so the
                        number of families is meaningless here; only timing and
                        memory should be read off this

All BGCs are phosphonate, so `--classify category` puts them in a single bin and
compares every pair — the worst case, which is the one worth sizing for.

    python scripts/bench_bigscape_scaling.py --sizes 333,500,750,1000,1500,2000
"""
import argparse
import csv
import json
import math
import os
import random
import re
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path


def unique_copy(src: Path, dest: Path, tag: str) -> None:
    """Copy a region GBK, stamping a unique identity so it is not deduplicated."""
    text = src.read_text()
    # LOCUS is fixed-width-ish but BiG-SCAPE parses via Bio.SeqIO, which tolerates
    # a longer name; keep the rest of the line intact.
    text = re.sub(r'^LOCUS       (\S+)', lambda m: f'LOCUS       {m.group(1)}_{tag}',
                  text, count=1, flags=re.M)
    text = re.sub(r'^ACCESSION   (\S+)', lambda m: f'ACCESSION   {m.group(1)}_{tag}',
                  text, count=1, flags=re.M)
    text = re.sub(r'^VERSION     (\S+)', lambda m: f'VERSION     {m.group(1)}_{tag}',
                  text, count=1, flags=re.M)
    dest.write_text(text)


def build_set(pool: list[Path], n: int, outdir: Path, seed: int) -> int:
    """Materialise n region GBKs into outdir/genome_XXXX/, returning the count."""
    if outdir.exists():
        shutil.rmtree(outdir)
    rng = random.Random(seed)
    # BiG-SCAPE walks **/*.gbk; mirror the antiSMASH layout of one dir per genome
    # so file discovery costs resemble a real run.
    for i in range(n):
        src = pool[i] if i < len(pool) else rng.choice(pool)
        gdir = outdir / f'g{i // 4:05d}'
        gdir.mkdir(parents=True, exist_ok=True)
        unique_copy(src, gdir / f'{src.stem}_c{i:06d}.region001.gbk', f'c{i:06d}')
    return n


def run_bigscape(bigscape: str, indir: Path, outdir: Path, pfam: Path,
                 cores: int, cutoff: str) -> dict:
    if outdir.exists():
        shutil.rmtree(outdir)
    cmd = [bigscape, 'cluster', '-i', str(indir), '-o', str(outdir),
           '--pfam-path', str(pfam), '--alignment-mode', 'auto',
           '--gcf-cutoffs', cutoff, '--include-singletons',
           '--cores', str(cores), '--classify', 'category']
    # BiG-SCAPE shells out to fasttree for its GCF trees. Nextflow gets that for
    # free by activating the conda env; invoking the binary by path does not, so
    # put its own bin directory on PATH or the run dies after the distances are
    # already computed.
    env = dict(os.environ)
    binf = Path(bigscape).resolve().parent
    env['PATH'] = f"{binf}{os.pathsep}{env.get('PATH', '')}"

    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    wall = time.monotonic() - t0
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        'wall_s': round(wall, 1),
        'cpu_s': round((after.ru_utime - before.ru_utime) +
                       (after.ru_stime - before.ru_stime), 1),
        # ru_maxrss is the largest single child, in KB on Linux — not a sum, so
        # it does not have the double-counting problem Nextflow's peak_rss has.
        'max_rss_gb': round(after.ru_maxrss / 1048576, 2),
        'exit': proc.returncode,
        'stderr_tail': proc.stderr.strip().splitlines()[-1][:200] if proc.stderr.strip() else '',
    }


def count_results(outdir: Path) -> dict:
    db = next(outdir.rglob('*.db'), None)
    if db is None:
        return {'bgcs': 0, 'families': 0, 'comparisons': 0}
    import sqlite3
    con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
    q = lambda s: con.execute(s).fetchone()[0]
    try:
        return {'bgcs': q("SELECT COUNT(*) FROM bgc_record WHERE record_type='region'"),
                'families': q("SELECT COUNT(DISTINCT id) FROM family"),
                # The validity check for the whole benchmark: if this is not
                # n(n-1)/2 then BiG-SCAPE pruned pairs and the curve means
                # something else entirely.
                'comparisons': q("SELECT COUNT(*) FROM distance")}
    except sqlite3.Error:
        return {'bgcs': 0, 'families': 0, 'comparisons': 0}
    finally:
        con.close()


def fit_models(rows: list[dict]) -> dict:
    """Fit both a power law and fixed+quadratic; the second is the one to use.

    Comparisons are all-pairs (verified against the distance table), but cost
    *per* comparison falls as n grows, so over a small range the runtime looks
    sub-quadratic and a log-log fit returns an exponent well under 2. Measured
    at 333-4000 BGCs it reads 1.245, which understates a 121k-BGC run tenfold.
    Extrapolate the fixed+quadratic model, and fit it on the largest points
    available, where the quadratic term actually dominates.
    """
    pts = [(r['n'], r['cpu_s']) for r in rows if r['n'] > 0 and r['cpu_s'] > 0]
    if len(pts) < 3:
        return {}

    def lsq(xs, ys):
        m = len(xs); mx = sum(xs) / m; my = sum(ys) / m
        den = sum((x - mx) ** 2 for x in xs)
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
        return my - b * mx, b

    def r2(ys, pred):
        my = sum(ys) / len(ys)
        tot = sum((y - my) ** 2 for y in ys)
        return 1 - sum((y - p) ** 2 for y, p in zip(ys, pred)) / tot if tot else None

    n = [p[0] for p in pts]; y = [p[1] for p in pts]
    la, lb = lsq([math.log(v) for v in n], [math.log(v) for v in y])
    power = {'exponent': round(lb, 3), 'coefficient': round(math.exp(la), 8),
             'r2': round(r2(y, [math.exp(la) * v ** lb for v in n]), 4)}

    # Fit the quadratic on the upper half, where fixed cost no longer dominates.
    upper = sorted(pts)[len(pts) // 2:] if len(pts) >= 6 else sorted(pts)
    ux = [v * v for v, _ in upper]; uy = [c for _, c in upper]
    a, b = lsq(ux, uy)
    worst = max(abs(c - (a + b * v * v)) / c for v, c in upper)
    quad = {'fixed_cpu_s': round(a, 1), 'per_pair_coeff': b,
            'fitted_on': [v for v, _ in upper], 'worst_residual': round(worst, 4)}
    return {'power_law': power, 'fixed_quadratic': quad}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--source', type=Path,
                    default=Path('results/antismash_results/Erwiniaceae'),
                    help='directory tree holding real *.region*.gbk files')
    ap.add_argument('--pfam', type=Path, help='Pfam-A.hmm (default: found under results/databases)')
    ap.add_argument('--bigscape', default='bigscape', help='bigscape executable')
    ap.add_argument('--sizes', default='333,500,750,1000,1500,2000',
                    help='comma-separated BGC counts to measure')
    ap.add_argument('--cores', type=int, default=8)
    ap.add_argument('--cutoff', default='0.30')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--workdir', type=Path,
                    default=Path('results/bench_bigscape'), help='scratch + results')
    ap.add_argument('--out', type=Path, help='TSV path (default: <workdir>/scaling.tsv)')
    args = ap.parse_args()

    pool = sorted(args.source.rglob('*.region*.gbk'))
    if not pool:
        print(f'no region GBKs under {args.source}', file=sys.stderr)
        return 1

    pfam = args.pfam
    if pfam is None:
        pfam = next(Path('results/databases').rglob('Pfam-A.hmm'), None)
    if pfam is None or not pfam.exists():
        print('Pfam-A.hmm not found; pass --pfam', file=sys.stderr)
        return 1

    sizes = [int(s) for s in args.sizes.split(',') if s.strip()]
    args.workdir.mkdir(parents=True, exist_ok=True)
    out = args.out or args.workdir / 'scaling.tsv'

    print(f'pool: {len(pool)} real BGCs from {args.source}')
    print(f'pfam: {pfam}')
    print(f'sizes: {sizes}  cores: {args.cores}\n')

    rows = []
    for n in sizes:
        indir = args.workdir / f'in_{n}'
        odir = args.workdir / f'out_{n}'
        print(f'[{n:>6} BGCs] building...', flush=True)
        build_set(pool, n, indir, args.seed)
        print(f'[{n:>6} BGCs] running bigscape...', flush=True)
        res = run_bigscape(args.bigscape, indir, odir, pfam, args.cores, args.cutoff)
        res.update(count_results(odir))
        res['n'] = n
        res['replication'] = round(n / len(pool), 2)
        rows.append(res)
        status = 'ok' if res['exit'] == 0 else f"EXIT {res['exit']} {res['stderr_tail']}"
        print(f"[{n:>6} BGCs] {res['wall_s']:>8.1f}s wall  {res['cpu_s']:>9.1f}s cpu  "
              f"{res['max_rss_gb']:>6.2f} GB  loaded={res['bgcs']}  {status}\n", flush=True)
        # Reclaim the staged copies immediately; at 2000 BGCs these are ~120 MB each.
        shutil.rmtree(indir, ignore_errors=True)

        # Merge rather than overwrite. A later --sizes run must not destroy the
        # earlier points; the fit needs every size that was ever measured.
        cols = ['n', 'replication', 'bgcs', 'families', 'comparisons', 'wall_s',
                'cpu_s', 'max_rss_gb', 'exit', 'stderr_tail']
        merged = {}
        if out.exists():
            with out.open(newline='') as fh:
                for old in csv.DictReader(fh, delimiter='\t'):
                    try:
                        merged[int(old['n'])] = {k: old.get(k, '') for k in cols}
                    except (TypeError, ValueError):
                        continue
        for r in rows:
            merged[r['n']] = r
        with out.open('w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            w.writerows(merged[k] for k in sorted(merged))
        all_rows = [{'n': int(v['n']), 'cpu_s': float(v['cpu_s'] or 0),
                     'wall_s': float(v['wall_s'] or 0), 'exit': int(v['exit'] or 1)}
                    for v in (merged[k] for k in sorted(merged))]

    fit = fit_models([r for r in all_rows if r['exit'] == 0])
    if fit:
        (args.workdir / 'fit.json').write_text(json.dumps(fit, indent=2, default=float))
        pw, qd = fit['power_law'], fit['fixed_quadratic']
        print(f"\nfitted over {len(all_rows)} sizes: {[r['n'] for r in all_rows]}")
        print(f"  power law       cpu_s = {pw['coefficient']} * n^{pw['exponent']} "
              f"(R2={pw['r2']})  <- DO NOT extrapolate, see fit_models()")
        print(f"  fixed+quadratic cpu_s = {qd['fixed_cpu_s']} + "
              f"{qd['per_pair_coeff']:.4e} * n^2   fitted on {qd['fitted_on']}, "
              f"worst residual {qd['worst_residual']:.1%}")
        for target in (121_000, 184_000):
            h = (qd['fixed_cpu_s'] + qd['per_pair_coeff'] * target ** 2) / 3600
            print(f"  projected {target:>7,} BGCs: {h:>10,.0f} CPU-h"
                  f"   ({h / args.cores:>7,.0f} wall-h at {args.cores} cores)")
    print(f'\nwrote {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
