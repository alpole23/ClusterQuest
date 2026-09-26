#!/usr/bin/env python3
"""Fetch the pepM-bearing MIBiG clusters and make BiG-SCAPE able to read them.

MIBiG carries 13 clusters with a pepM homologue, found by searching MIBiG 4.0's
proteins with the 7 curated pepM references (five more than a keyword search
finds). Bringing them in with `--bigscape_mibig_version 4.0` also brings ~2,400
unrelated clusters into the same distance matrix, which shifts family
composition — that is why that option is off by default. Supplying just these
through `--reference-dir` instead measures the same distances without touching
the published clustering.

**BiG-SCAPE 2.0.1 cannot read a MIBiG 4.0 GenBank as published**, and patching
the symptoms does not work. MIBiG writes `Version :: False` in the
antiSMASH-Data header and a `region` feature with no `candidate_cluster_numbers`.
Fixing those two gets you a third failure, because BiG-SCAPE's AS5 reader walks
a four-level hierarchy -- region, cand_cluster, protocluster, proto_core -- and
MIBiG supplies only the first. `make_reference_bgc.py` already builds the whole
chain, and its own comment warns against discovering the levels one rejection at
a time, which is precisely what repairing in place turns into.

So MIBiG's partial region is STRIPPED and the hierarchy rebuilt by that tool.
The MIBiG boundary is preserved: it declared the whole record to be the cluster,
and so does the rebuild.

The quieter trap is the filename: BiG-SCAPE only ingests `.gbk` files whose names
contain "cluster" or "region", so `BGC0000897.gbk` would be skipped in silence.
Everything written here is named `<accession>_<product>.region001.gbk`.

    python scripts/genome/fetch_mibig_references.py \\
        --outdir assets/phosphonate_reference_bgcs
"""
import argparse
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ARCHIVE = 'https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz'

# The 13 pepM-bearing MIBiG 4.0 clusters, minus one. BGC0000383 is EXCLUDED and
# the reason is worth recording: it is deposited as the luminmycin/glidobactin
# NRPS/PKS cluster of Photorhabdus, and the pepM is there because a pantaphos-like
# BGC sits adjacent in the deposit, unnoticed by its authors. Including it would
# label a phosphonate distance with an unrelated product. It also explains a
# result this repository already measured and could not account for: in the
# Erwiniaceae KnownClusterBlast run, luminmycin/glidobactin was the best hit on
# 236 of 334 regions.
CLUSTERS = {
    'BGC0000806': 'phosphonoglycans_Glycomyces',
    'BGC0000807': 'phosphonoglycans_Stackebrandtia',
    'BGC0000897': 'dehydrophos',
    'BGC0000904': 'FR900098',
    'BGC0000926': 'rhizocticin_A',
    'BGC0000937': 'fosfazinomycin',
    'BGC0000938': 'fosfomycin_Sfradiae',
    'BGC0001411': 'polysaccharideB_Bfragilis_2AEP',
    'BGC0001739': 'phosphonoacetic_acid',
    'BGC0001859': 'fosfomycin_Psyringae',
    'BGC0002036': 'dehydrofosmidomycin',
    'BGC0002670': 'fosfonochlorin',
}

STRIP = {'region', 'cand_cluster', 'protocluster', 'proto_core'}


def strip_region_features(raw, tmpdir, accession):
    """MIBiG's partial antiSMASH hierarchy removed, so the tool will rebuild it."""
    from Bio import SeqIO
    src = Path(tmpdir) / f'{accession}.raw.gbk'
    src.write_text(raw)
    recs = list(SeqIO.parse(str(src), 'genbank'))
    removed = 0
    for rec in recs:
        keep = [f for f in rec.features if f.type not in STRIP]
        removed += len(rec.features) - len(keep)
        rec.features = keep
    out = Path(tmpdir) / f'{accession}.stripped.gbk'
    SeqIO.write(recs, str(out), 'genbank')
    return out, removed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--archive', help='local mibig_gbk_4.0.tar.gz; downloaded if absent')
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tar_path = Path(args.archive) if args.archive else Path(tmp) / 'mibig.tar.gz'
        if not args.archive:
            print(f'downloading {ARCHIVE}')
            subprocess.run(['curl', '-sSL', '-o', str(tar_path), ARCHIVE], check=True)
        with tarfile.open(tar_path) as tf:
            written = 0
            for acc, product in sorted(CLUSTERS.items()):
                member = f'mibig_gbk_4.0/{acc}.gbk'
                try:
                    raw = tf.extractfile(member).read().decode()
                except (KeyError, AttributeError):
                    print(f'  {acc}: NOT in the archive — has MIBiG renumbered it?')
                    continue
                stripped, removed = strip_region_features(raw, tmp, acc)
                dest = args.outdir / f'{acc}_{product}.region001.gbk'
                note = (f'MIBiG 4.0 {acc}. Its own region feature was stripped and the '
                        f'antiSMASH hierarchy rebuilt: MIBiG supplies only a region, '
                        f'while BiG-SCAPE 2.0.1 walks region/cand_cluster/protocluster/'
                        f'proto_core. Boundary unchanged -- MIBiG declares the whole '
                        f'record to be the cluster.')
                r = subprocess.run(
                    [sys.executable, str(Path(__file__).with_name('make_reference_bgc.py')),
                     '--in', str(stripped), '--out', str(dest), '--note', note],
                    capture_output=True, text=True)
                if r.returncode != 0:
                    print(f'  {acc}: make_reference_bgc.py failed — '
                          f'{r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "?"}')
                    continue
                written += 1
                print(f'  {acc} -> {dest.name}  ({removed} MIBiG features stripped, '
                      f'hierarchy rebuilt)')

    print(f'\n{written} of {len(CLUSTERS)} written to {args.outdir}')
    return 0 if written == len(CLUSTERS) else 1


if __name__ == '__main__':
    sys.exit(main())
