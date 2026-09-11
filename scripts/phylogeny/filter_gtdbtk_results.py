#!/usr/bin/env python3
"""
Filter GTDB-Tk results to only include genomes from a subset run.

This script filters the GTDB-Tk summary TSV to only include genomes specified in
the genome list. Used for cross-taxon result reuse, where a larger taxon run's
results are filtered for a subset taxon.

Tree pruning was removed with the GTDB-Tk tree itself: GTDBTK_CLASSIFY is sharded
for wall time and per-shard trees span disjoint genome sets, so no run produces a
whole-set tree to prune.

Usage:
    filter_gtdbtk_results.py <genome_list> <reuse_summary> <output_dir>

Arguments:
    genome_list   - File with list of genome paths (one per line)
    reuse_summary - Path to source GTDB-Tk summary TSV
    output_dir    - Output directory for filtered results
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter GTDB-Tk results for a subset of genomes'
    )
    parser.add_argument('genome_list', type=Path,
                        help='File with list of genome paths')
    parser.add_argument('reuse_summary', type=Path,
                        help='Source GTDB-Tk summary TSV')
    parser.add_argument('output_dir', type=Path,
                        help='Output directory')
    return parser.parse_args()


def read_genome_list(genome_list_path: Path) -> set:
    """
    Read genome list and create set of valid genome names.
    Includes both raw names and usr_ prefixed versions (GTDB-Tk adds this prefix).
    """
    current_genomes = set()
    with open(genome_list_path) as f:
        for line in f:
            # Extract basename without .fna extension
            genome = Path(line.strip()).stem
            current_genomes.add(genome)
            # Also add with usr_ prefix (GTDB-Tk adds this)
            current_genomes.add(f"usr_{genome}")
    return current_genomes


def filter_summary(reuse_summary: Path, current_genomes: set, output_path: Path):
    """Filter GTDB-Tk summary TSV to only include specified genomes."""
    import pandas as pd

    summary_df = pd.read_csv(reuse_summary, sep='\t')
    filtered_df = summary_df[summary_df['user_genome'].isin(current_genomes)]

    print(f"Filtered summary: {len(filtered_df)} rows (from {len(summary_df)})")
    filtered_df.to_csv(output_path, sep='\t', index=False)


def main():
    args = parse_args()

    # Setup output directories
    outdir = args.output_dir
    outdir.mkdir(exist_ok=True)

    # Read genome list
    current_genomes = read_genome_list(args.genome_list)
    print(f"Filtering results for {len(current_genomes)//2} genomes")

    # Filter summary
    summary_output = outdir / "gtdbtk.bac120.summary.tsv"
    filter_summary(args.reuse_summary, current_genomes, summary_output)

    print("GTDB-Tk result filtering complete")


if __name__ == '__main__':
    main()
