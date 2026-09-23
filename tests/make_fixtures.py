#!/usr/bin/env python3
"""Create the synthetic genome fixture used by the batching test.

Writes <outdir>/data/<assembly_id>/genomic.gbff for a handful of tiny genomes plus a
name_map.json, mirroring the layout NCBI_DATASETS_DOWNLOAD produces. One genome is
deliberately corrupt so the test can check that a batch tolerates a bad member.
"""

import argparse
import json
import os

GENBANK_TEMPLATE = """LOCUS       CONTIG_{n}                60 bp    DNA     linear   BCT 01-JAN-2020
DEFINITION  test genome {n}.
ACCESSION   CONTIG_{n}
VERSION     CONTIG_{n}.1
FEATURES             Location/Qualifiers
     source          1..60
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--n', type=int, default=7, help='Number of valid genomes')
    parser.add_argument('--with-corrupt', action='store_true',
                        help='Add one unparseable genome to test batch tolerance')
    args = parser.parse_args()

    name_map = {}
    for i in range(1, args.n + 1):
        assembly_id = f"GCF_{i:09d}.1"
        genome_dir = os.path.join(args.outdir, 'data', assembly_id)
        os.makedirs(genome_dir, exist_ok=True)
        with open(os.path.join(genome_dir, 'genomic.gbff'), 'w') as f:
            f.write(GENBANK_TEMPLATE.format(n=i))
        name_map[assembly_id] = f"Test organism {i}"

    if args.with_corrupt:
        assembly_id = f"GCF_{args.n + 1:09d}.1"
        genome_dir = os.path.join(args.outdir, 'data', assembly_id)
        os.makedirs(genome_dir, exist_ok=True)
        with open(os.path.join(genome_dir, 'genomic.gbff'), 'w') as f:
            f.write("this is not a GenBank file\n")
        name_map[assembly_id] = "Broken organism"

    # Uniquely-named copies, as the pipeline's rename step would produce. NCBI
    # names every genome `genomic.gbff`, so anything staging several of them in one
    # task needs distinct names; GENBANK_TO_FASTA used to get these from
    # RENAME_GENOMES, which no longer exists as a stage.
    renamed_dir = os.path.join(args.outdir, 'renamed')
    os.makedirs(renamed_dir, exist_ok=True)
    for assembly_id, organism in name_map.items():
        src = os.path.join(args.outdir, 'data', assembly_id, 'genomic.gbff')
        with open(src) as fh, open(os.path.join(
                renamed_dir, organism.replace(' ', '_') + '.gbff'), 'w') as out:
            out.write(fh.read())

    # Fake antiSMASH results for the reuse-copy path, including a hidden metadata
    # file and a nested directory so the test can check copy fidelity.
    for genome in ('Genome_A', 'Genome_B'):
        d = os.path.join(args.outdir, 'reuse', genome)
        os.makedirs(os.path.join(d, 'nested'), exist_ok=True)
        with open(os.path.join(d, f'{genome}.json'), 'w') as f:
            json.dump({'records': []}, f)
        with open(os.path.join(d, '.antismash_meta'), 'w') as f:
            f.write("version=7.1.0\nparams_hash=deadbeef\n")
        with open(os.path.join(d, 'nested', 'deep.txt'), 'w') as f:
            f.write("nested payload\n")

    with open(os.path.join(args.outdir, 'name_map.json'), 'w') as f:
        json.dump(name_map, f)

    print(f"Wrote {len(name_map)} genomes to {args.outdir}/data")


if __name__ == '__main__':
    main()
