#!/usr/bin/env python3
"""Reclaim disk from published antiSMASH results without losing analysable data.

At 1,735 genomes the Pantoea run left 43 GB under ``antismash_results/``, and a
directory for a genome with no BGC is the same size as one with BGCs — ~23 MB either
way. The bulk is not regions:

    {genome}.gbk    8.6 MB   annotated genome (only if --antismash_summary_gbk)
    {genome}.json   6.7 MB   antiSMASH JSON — required by CHECK_ANTISMASH_REUSE
    {genome}.zip    5.6 MB   archive of the very same directory
    js/images/css   708 KB   byte-identical in every genome's directory

So the largest safe win is not "delete BGC-negative genomes" — it is deleting the
redundant ``.zip`` from *every* directory, which costs nothing at all.

Three tiers, increasingly aggressive:

  archives   drop {genome}.zip everywhere. Lossless: it archives the loose files
             sitting beside it. ~24% of the tree. Retrospective only — ANTISMASH now
             passes --no-zip-output, so runs from that change onward never write one
             and this tier finds nothing to reclaim in their output.
  strip      for BGC-negative genomes, also drop the .gbk and the HTML report assets.
             The .gbk half is retrospective: ANTISMASH now defaults to
             --no-summary-gbk, so new runs only have one unless the user opted in.
             but KEEP {genome}.json and .antismash_meta so --reuse_antismash_from still
             recognises the genome as analysed and does not re-run antiSMASH on it.
  purge      remove BGC-negative directories outright. Frees the most, and breaks
             reuse for those genomes: a later run with --reuse_antismash_from will
             find nothing and re-run antiSMASH, which for a low-prevalence taxon is
             most of the compute bill.

Dry-run by default. Nothing is deleted without --apply.

    python scripts/prune_antismash_results.py \
        --antismash_dir results/antismash_results/Pantoea \
        --counts results/main_analysis_results/Pantoea/region_counts.tsv

Note this prunes only the published copy. Nextflow's work/ holds another copy of every
result (publishDir mode 'copy'), and accumulates one per run — reclaim those with
`nextflow clean -f` once you no longer need to -resume.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Kept for every genome, whatever the tier: without the JSON, CHECK_ANTISMASH_REUSE
# reports "no existing results" and antiSMASH re-runs from scratch.
# {genome}.json and the dot-file .antismash_meta. The `strip` tier never targets
# either, so they survive by omission rather than by an explicit keep-list.
REUSE_ESSENTIAL = ('{genome}.json', '.antismash_meta')
HTML_ASSET_DIRS = ('css', 'js', 'images')
HTML_ASSET_FILES = ('index.html', 'regions.js')


def bgc_status(counts_file):
    """{genome_name: total_bgcs} from region_counts.tsv."""
    with open(counts_file) as f:
        rows = [ln for ln in f if not ln.startswith('#')]
    out = {}
    for r in csv.DictReader(rows, delimiter='\t'):
        name = (r.get('record') or '').split('|')[0]
        for suffix in ('.gbff', '.gbk', '.gb'):
            if name.lower().endswith(suffix):
                name = name[:-len(suffix)]
                break
        try:
            out[name] = int(r.get('total_count') or 0)
        except ValueError:
            out[name] = 0
    return out


def dir_size(path):
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


def pipeline_running():
    """True if a Nextflow run looks active — pruning published output under a live run
    races with publishDir and can delete files it is still writing."""
    try:
        out = subprocess.run(['pgrep', '-af', 'nextflow run'],
                             capture_output=True, text=True, timeout=10)
        return bool(out.stdout.strip())
    except Exception:
        return False


def plan(antismash_dir, status, tier):
    """[(path, bytes, reason)] of everything the chosen tier would remove."""
    targets = []
    for entry in sorted(Path(antismash_dir).iterdir()):
        if not entry.is_dir():
            continue
        genome = entry.name
        negative = status.get(genome, 0) == 0
        known = genome in status

        # purge removes the whole directory, so its .zip must not also be counted
        # separately — that double-counts and understates what purge actually frees
        if tier == 'purge' and negative and known:
            targets.append((entry, dir_size(entry), 'BGC-negative, whole directory'))
            continue

        zip_file = entry / f'{genome}.zip'
        if zip_file.exists():
            targets.append((zip_file, zip_file.stat().st_size, 'redundant archive'))

        if tier == 'archives' or not negative or not known:
            continue

        # tier == 'strip': drop the bulky annotated genome and the report chrome,
        # keep the JSON and meta so reuse still recognises this genome
        gbk = entry / f'{genome}.gbk'
        if gbk.exists():
            targets.append((gbk, gbk.stat().st_size, 'BGC-negative, annotated genome'))
        for d in HTML_ASSET_DIRS:
            p = entry / d
            if p.is_dir():
                targets.append((p, dir_size(p), 'BGC-negative, report assets'))
        for f in HTML_ASSET_FILES:
            p = entry / f
            if p.exists():
                targets.append((p, p.stat().st_size, 'BGC-negative, report page'))
    return targets


def human(n):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if abs(n) < 1024:
            return f'{n:.1f} {unit}'
        n /= 1024
    return f'{n:.1f} PB'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--antismash_dir', required=True, help='results/antismash_results/<taxon>')
    ap.add_argument('--counts', required=True, help='region_counts.tsv for the same taxon')
    ap.add_argument('--tier', choices=('archives', 'strip', 'purge'), default='strip',
                    help='archives: .zip only (lossless). strip: also thin BGC-negative '
                         'genomes but keep reuse working (default). purge: delete '
                         'BGC-negative directories, breaking reuse for them.')
    ap.add_argument('--apply', action='store_true',
                    help='actually delete; without it nothing is removed')
    ap.add_argument('--force', action='store_true',
                    help='proceed even if a Nextflow run appears to be active')
    args = ap.parse_args()

    ad = Path(args.antismash_dir)
    if not ad.is_dir():
        print(f'No such directory: {ad}', file=sys.stderr)
        return 1
    if not Path(args.counts).exists():
        print(f'No such counts file: {args.counts}', file=sys.stderr)
        return 1

    if args.apply and pipeline_running() and not args.force:
        print('A Nextflow run appears to be active. Pruning published results while\n'
              'publishDir is writing can delete files mid-copy. Wait for it to finish,\n'
              'or pass --force if you are certain.', file=sys.stderr)
        return 2

    status = bgc_status(args.counts)
    negatives = sum(1 for v in status.values() if v == 0)
    before = dir_size(ad)
    targets = plan(ad, status, args.tier)
    freed = sum(sz for _, sz, _ in targets)

    by_reason = {}
    for _, sz, reason in targets:
        by_reason[reason] = by_reason.get(reason, [0, 0])
        by_reason[reason][0] += 1
        by_reason[reason][1] += sz

    print(f'  genomes in counts   : {len(status)}  ({negatives} BGC-negative)')
    print(f'  antismash_dir size  : {human(before)}')
    print(f'  tier                : {args.tier}')
    print()
    for reason, (n, sz) in sorted(by_reason.items(), key=lambda x: -x[1][1]):
        print(f'    {reason:34} {n:>6} items  {human(sz):>10}')
    print(f'\n  would free          : {human(freed)}  ({100 * freed / before:.0f}% of the tree)')
    print(f'  remaining           : {human(before - freed)}')

    if not args.apply:
        print('\n  DRY RUN — nothing deleted. Re-run with --apply to remove.')
        return 0

    removed = 0
    for path, sz, _ in targets:
        try:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
            removed += sz
        except OSError as e:
            print(f'  could not remove {path}: {e}', file=sys.stderr)
    print(f'\n  freed {human(removed)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
