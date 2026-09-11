#!/usr/bin/env python3
"""Convert GenBank files to FASTA format (one .fna per input)."""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO


def convert(input_file, output_file):
    """Convert one GenBank file. Returns True on success."""
    records = list(SeqIO.parse(input_file, "genbank"))

    if not records:
        print(f"ERROR: No sequences found in {input_file}")
        return False

    with open(output_file, 'w') as f:
        SeqIO.write(records, f, "fasta")

    total_bp = sum(len(r.seq) for r in records)
    print(f"Converted {input_file} -> {output_file}")
    print(f"  Sequences: {len(records)}")
    print(f"  Total length: {total_bp:,} bp")
    return True


def main():
    parser = argparse.ArgumentParser(description="Convert GenBank files to FASTA format")
    parser.add_argument("input_files", nargs='+', help="Input GenBank file(s)")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    args = parser.parse_args()

    failed = []
    for input_file in args.input_files:
        output_file = Path(args.output_dir) / (Path(input_file).stem + ".fna")
        try:
            if not convert(input_file, output_file):
                failed.append(input_file)
        except Exception as e:
            print(f"ERROR: Failed to convert {input_file}: {e}")
            failed.append(input_file)

    print(f"Converted {len(args.input_files) - len(failed)}/{len(args.input_files)} genomes")

    # Only fail the batch if every genome failed; individual failures are tolerated
    # (matches the per-genome 'tolerant' behaviour this process had before batching).
    if failed and len(failed) == len(args.input_files):
        sys.exit(1)


if __name__ == "__main__":
    main()
