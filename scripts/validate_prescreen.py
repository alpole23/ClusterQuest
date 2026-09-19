#!/usr/bin/env python3
"""Validate the pepM pre-screen against an unscreened arm of the same genomes.

The screen decides, per genome, whether antiSMASH is worth running. The only
honest way to grade that decision is to run antiSMASH on everything and see
what it would have thrown away, which is what the unscreened arm is for.

    python scripts/validate_prescreen.py \
        --ground-truth results_heldout_off --screened results_heldout_on \
        --taxon "Bacteroides fragilis" \
        --outdir docs/comparisons/pepm_prescreen/heldout_bacteroides

Ground truth is "antiSMASH called at least one region in this genome", read
from the unscreened run's region_tabulation.tsv. Verdicts come from the
screened run's prescreen_*.tsv. A genome present in one arm and not the other
is reported rather than dropped -- the arms are supposed to see the same set,
and a silent set difference would corrupt every rate below.
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def sanitize_taxon(name):
    """Mirror Utils.sanitizeTaxon: the on-disk spelling of a taxon."""
    return ''.join(c if c.isalnum() or c in '_-' else '_' for c in name).strip('_')


def read_verdicts(results, taxon):
    """Every prescreen_*.tsv shard of one run, keyed by genome."""
    root = Path(results) / 'prescreen_results' / sanitize_taxon(taxon)
    shards = sorted(root.glob('prescreen_*.tsv'))
    if not shards:
        sys.exit(f"no prescreen output under {root} -- was the screen on for that run?")
    verdicts = {}
    for shard in shards:
        with open(shard) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                verdicts[row['genome']] = {
                    'mode': row['mode'],
                    'cds_per_mb': int(float(row['cds_per_mb'])),
                    'bitscore': float(row['best_bitscore']),
                    'pass': row['pass'].strip().lower() == 'yes',
                }
    return verdicts


def read_bgc_positive(results, taxon):
    """Genomes antiSMASH called at least one region in, and the region count."""
    tab = Path(results) / 'main_analysis_results' / sanitize_taxon(taxon) / 'region_tabulation.tsv'
    if not tab.exists():
        sys.exit(f"no region tabulation at {tab}")
    regions = {}
    with open(tab) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            regions[row['file']] = regions.get(row['file'], 0) + 1
    return regions


def read_analysed(results, taxon):
    """Every genome the run analysed, from region_counts.tsv (BGC-negative included)."""
    counts = Path(results) / 'main_analysis_results' / sanitize_taxon(taxon) / 'region_counts.tsv'
    if not counts.exists():
        return set()
    names = set()
    with open(counts) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('record\t'):
                continue
            record = line.split('\t')[0].strip()
            if record:
                names.add(strip_genome_suffix(record))
    return names


def strip_genome_suffix(name):
    """region_counts.tsv carries the extension; nothing else does."""
    for suffix in ('.gbff', '.gbk', '.gb', '.genbank'):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def best_cuts(scored):
    """Where the cutoff could sit: the highest lossless one, the cleanest one.

    `scored` is (bitscore, is_true_positive). Returned cuts are expressed as
    the bitscore a genome must reach, so a cut of c keeps score >= c.
    """
    positives = sorted(s for s, tp in scored if tp)
    negatives = sorted(s for s, tp in scored if not tp)
    if not positives:
        return {'note': 'no true positives; sensitivity is undefined on this set'}
    weakest_tp = positives[0]
    lossless = {
        'highest_lossless_cut': weakest_tp,
        'false_positives_there': sum(1 for s in negatives if s >= weakest_tp),
    }
    clean = None
    top_negative = negatives[-1] if negatives else 0.0
    kept = sum(1 for s in positives if s > top_negative)
    clean = {
        'zero_fp_cut': top_negative + 0.1 if negatives else 0.0,
        'true_positives_kept_there': f'{kept}/{len(positives)}',
    }
    return {**lossless, **clean,
            'true_positive_range': [positives[0], positives[-1]],
            'negative_range': [negatives[0], negatives[-1]] if negatives else None}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ground-truth', required=True,
                    help='results dir of the UNSCREENED arm (pepm_prescreen=false)')
    ap.add_argument('--screened', required=True,
                    help='results dir of the SCREENED arm (pepm_prescreen=true)')
    ap.add_argument('--taxon', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--cutoff', type=float, default=100.0,
                    help='the cutoff the screened run used (default: 100)')
    ap.add_argument('--question', default='Does the pepM pre-screen hold on this clade?')
    ap.add_argument('--dataset', default='')
    args = ap.parse_args()

    verdicts = read_verdicts(args.screened, args.taxon)
    regions = read_bgc_positive(args.ground_truth, args.taxon)
    analysed = read_analysed(args.ground_truth, args.taxon) or set(verdicts)

    # A set difference between the arms invalidates every rate below, so say so.
    only_truth = sorted(analysed - set(verdicts))
    only_screen = sorted(set(verdicts) - analysed)

    scored, missed, false_positives = [], [], []
    for genome, v in sorted(verdicts.items()):
        is_tp = genome in regions
        scored.append((v['bitscore'], is_tp))
        if is_tp and not v['pass']:
            missed.append({'genome': genome, 'bitscore': v['bitscore'],
                           'mode': v['mode'], 'regions': regions[genome]})
        if not is_tp and v['pass']:
            false_positives.append({'genome': genome, 'bitscore': v['bitscore'],
                                    'mode': v['mode']})

    n_positive = sum(1 for _, tp in scored if tp)
    n_negative = len(scored) - n_positive
    passed = sum(1 for v in verdicts.values() if v['pass'])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(outdir / 'prescreen_verdicts.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['genome', 'mode', 'cds_per_mb', 'best_bitscore', 'pass',
                    'bgc_positive', 'regions'])
        for genome, v in sorted(verdicts.items()):
            w.writerow([genome, v['mode'], v['cds_per_mb'], v['bitscore'],
                        'yes' if v['pass'] else 'no',
                        'yes' if genome in regions else 'no',
                        regions.get(genome, 0)])

    summary = {
        'question': args.question,
        'date': __import__('datetime').date.today().isoformat(),
        'taxon': args.taxon,
        'dataset': args.dataset,
        'method': ('Ground truth from an unscreened arm - every genome through antiSMASH - '
                   'against the screen per-genome verdicts from the screened arm.'),
        'arms_agree_on_genome_set': not (only_truth or only_screen),
        'genomes': {
            'screened_arm': len(verdicts),
            'unscreened_arm': len(analysed),
            'only_in_unscreened': only_truth[:20],
            'only_in_screened': only_screen[:20],
        },
        'ground_truth': {
            'bgc_positive_genomes': n_positive,
            'regions': sum(regions.values()),
            'of_genomes': len(analysed),
        },
        'result': {
            'sensitivity': f'{n_positive - len(missed)}/{n_positive}',
            'missed': missed,
            'passed_the_screen': passed,
            'false_positives': len(false_positives),
            'negatives': n_negative,
            'false_positive_rate_pct': round(100.0 * len(false_positives) / n_negative, 2)
                                       if n_negative else None,
            'cutoff': args.cutoff,
        },
        'cutoff_analysis': best_cuts(scored),
    }

    with open(outdir / 'summary.json', 'w') as fh:
        json.dump(summary, fh, indent=2)
        fh.write('\n')

    r = summary['result']
    print(f"genomes            {len(verdicts)} screened, {n_positive} BGC-positive "
          f"({sum(regions.values())} regions)")
    print(f"sensitivity        {r['sensitivity']}")
    print(f"passed the screen  {passed}  ({len(false_positives)} false positives "
          f"of {n_negative} negatives"
          + (f", {r['false_positive_rate_pct']}%)" if n_negative else ")"))
    ca = summary['cutoff_analysis']
    if 'highest_lossless_cut' in ca:
        print(f"true positives     {ca['true_positive_range'][0]}-{ca['true_positive_range'][1]}")
        if ca['negative_range']:
            print(f"negatives          {ca['negative_range'][0]}-{ca['negative_range'][1]}")
        print(f"highest lossless   {ca['highest_lossless_cut']} "
              f"({ca['false_positives_there']} false positives there)")
        print(f"zero-FP cut        {ca['zero_fp_cut']} "
              f"(keeps {ca['true_positives_kept_there']})")
    if not summary['arms_agree_on_genome_set']:
        print(f"\nWARNING: the two arms do not cover the same genomes "
              f"({len(only_truth)} only unscreened, {len(only_screen)} only screened). "
              f"Every rate above is computed over the screened arm regardless.")
    print(f"\nwrote {outdir}/summary.json and prescreen_verdicts.tsv")


if __name__ == '__main__':
    main()
