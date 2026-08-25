#!/usr/bin/env python3
"""Rename a batch of genome files based on the name map.

Reads a manifest of "assembly_id<TAB>staged_file" lines and copies each genome to
<mapped_name>.gbff in the output directory.
"""

import argparse
import json
import sys
from pathlib import Path
import shutil


def load_manifest(manifest_file):
    pairs = []
    with open(manifest_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            assembly_id, genome_file = line.split('\t')
            pairs.append((assembly_id, genome_file))
    return pairs


def rename_genome(assembly_id, genome_file, name_map, output_dir):
    """Copy one genome to its mapped name. Returns True on success."""
    if assembly_id not in name_map:
        print(f"ERROR: Assembly ID {assembly_id} not found in name map")
        return False

    new_name = name_map[assembly_id].replace(' ', '_')
    output_file = Path(output_dir) / f"{new_name}.gbff"

    shutil.copy2(genome_file, output_file)

    if not output_file.exists():
        print(f"ERROR: Output file was not created for {assembly_id}")
        return False

    file_size = output_file.stat().st_size
    if file_size == 0:
        print(f"ERROR: Output file is empty for {assembly_id}")
        return False

    print(f"Successfully created {output_file} ({file_size:,} bytes)")
    return True


def main():
    parser = argparse.ArgumentParser(description="Rename genome files based on name map")
    parser.add_argument("manifest", help="TSV of assembly_id<TAB>genome_file")
    parser.add_argument("name_map", help="Path to name_map.json")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    args = parser.parse_args()

    with open(args.name_map) as f:
        name_map = json.load(f)

    pairs = load_manifest(args.manifest)
    failed = [a for a, g in pairs if not rename_genome(a, g, name_map, args.output_dir)]

    print(f"Renamed {len(pairs) - len(failed)}/{len(pairs)} genomes")

    # Fail the batch only if nothing could be renamed; individual misses are
    # reported above and simply drop out of the run.
    if failed and len(failed) == len(pairs):
        sys.exit(1)


if __name__ == "__main__":
    main()
