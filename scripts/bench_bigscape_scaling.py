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
        return {'bgcs': 0, 'families': 0}
    import sqlite3
    con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
    q = lambda s: con.execute(s).fetchone()[0]
    try:
        return {'bgcs': q("SELECT COUNT(*) FROM bgc_record WHERE record_type='region'"),
                'families': q("SELECT COUNT(DISTINCT id) FROM family")}
    except sqlite3.Error:
        return {'bgcs': 0, 'families': 0}
    finally:
        con.close()


def fit_exponent(rows: list[dict]) -> dict:
    """Least-squares slope of log(time) against log(n) — the empirical exponent."""
    pts = [(math.log(r['n']), math.log(r['wall_s']))
           for r in rows if r['n'] > 0 and r['wall_s'] > 0]
    if len(pts) < 3:
        return {}
    n = len(pts)
    mx = sum(x for x, _ in pts) / n
    my = sum(y for _, y in pts) / n
    denom = sum((x - mx) ** 2 for x, _ in pts)
    slope = sum((x - mx) * (y - my) for x, y in pts) / denom
    intercept = my - slope * mx
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in pts)
    ss_tot = sum((y - my) ** 2 for _, y in pts)
    return {'exponent': round(slope, 3),
            'r2': round(1 - ss_res / ss_tot, 4) if ss_tot else None,
            'coefficient': round(math.exp(intercept), 8)}


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

        cols = ['n', 'replication', 'bgcs', 'families', 'wall_s', 'cpu_s',
                'max_rss_gb', 'exit', 'stderr_tail']
        with out.open('w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)

    fit = fit_exponent([r for r in rows if r['exit'] == 0])
    if fit:
        (args.workdir / 'fit.json').write_text(json.dumps(fit, indent=2))
        print(f"log-log fit: wall_s = {fit['coefficient']} * n^{fit['exponent']}  "
              f"(R2={fit['r2']})")
        e = fit['exponent']
        print(f"  1.0 = linear, 2.0 = all-pairs. Measured {e}.")
        for target in (121_000, 184_000):
            hours = fit['coefficient'] * target ** e / 3600
            print(f"  projected {target:>7,} BGCs: {hours:>10,.0f} wall-hours "
                  f"at {args.cores} cores")
    print(f'\nwrote {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
