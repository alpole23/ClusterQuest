#!/usr/bin/env python3
"""Re-cut antiSMASH region GenBanks to a narrower neighbourhood, without re-running
antiSMASH.

Why this exists: `antismash_phosphonate_neighbourhood` sets how much flank a region
carries, and the flank is 79% of what BiG-SCAPE compares (median 4 Pfam domains inside
the rule core against 19 across the whole region). Asking "what would the families look
like at 8 kb" therefore needs a different *comparison*, not a different *detection* --
and detection genuinely does not move with the neighbourhood: the same 334 regions in
the same 299 genomes were found at 5 kb and at 10 kb on Erwiniaceae. So the honest and
cheap experiment is to take the regions antiSMASH already found at 10 kb and re-cut
them to core +/- N, then re-run BiG-SCAPE over the result.

**This does not re-use the stored Pfam scan, and must not.** The obvious shortcut --
load the records from a finished database, narrow each one's nt_start/nt_stop, and call
BiG-SCAPE's comparison code on them -- was tried and abandoned. That comparison code
reproduces stored distances exactly (60 of 60 verified), but narrowing a record in place
desynchronises the CDS indices its LCS extension walks, and it raised `KeyError` inside
`score_extend_rev`. A crash is the good outcome there; the same inconsistency could as
easily have returned plausible wrong distances. Write real files and let BiG-SCAPE scan
them.

Usage:
    python trim_regions.py --in  results/antismash_results/Enterobacterales \\
                           --out trimmed_8kb --kb 8
    bigscape cluster -i trimmed_8kb -o out8 --pfam-path .../Pfam-A.hmm \\
        --alignment-mode auto --gcf-cutoffs 0.2,0.3,0.4 --include-singletons \\
        --classify category --cores 16
"""
import argparse
import os
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation

# The features BiG-SCAPE's parser walks. A sliced record keeps only features wholly
# inside the slice, and every one of these spans the old record, so all four would be
# dropped and the file would be rejected as non-antiSMASH. They are rebuilt by hand.
STRUCTURAL = ('region', 'cand_cluster', 'protocluster', 'proto_core')


def core_bounds(record):
    """(start, stop) of the rule core, in record coordinates."""
    for feat in record.features:
        if feat.type == 'proto_core':
            return int(feat.location.start), int(feat.location.end)
    return None


def trim(record, kb):
    """Return the record re-cut to core +/- kb, or None if it has no core."""
    bounds = core_bounds(record)
    if bounds is None:
        return None
    cs, ce = bounds
    lo = max(0, cs - kb * 1000)
    hi = min(len(record), ce + kb * 1000)

    structural = {f.type: f for f in record.features if f.type in STRUCTURAL}
    out = record[lo:hi]                     # Biopython shifts the CDS coordinates
    out.id, out.name, out.description = record.id, record.name, record.description
    out.annotations = dict(record.annotations)

    # Rebuild the four structural features against the new coordinates, clamped. The
    # region and cand_cluster span the whole trimmed record by definition; the
    # protocluster and proto_core keep their real extents, shifted.
    rebuilt = []
    for kind in STRUCTURAL:
        feat = structural.get(kind)
        if feat is None:
            continue
        if kind in ('region', 'cand_cluster'):
            start, end = 0, hi - lo
        else:
            start = max(0, int(feat.location.start) - lo)
            end = min(hi - lo, int(feat.location.end) - lo)
        quals = dict(feat.qualifiers)
        if kind == 'protocluster':
            # Keep the strand suffix antiSMASH writes, e.g. "[5000:8337](+)". BiG-SCAPE
            # reads the proto_core FEATURE rather than this qualifier, but the report's
            # own scripts parse it, so the file stays faithful rather than merely
            # parseable.
            old = (feat.qualifiers.get('core_location') or [''])[0]
            strand = old[old.index(']') + 1:] if ']' in old else ''
            quals['core_location'] = [f'[{cs - lo}:{ce - lo}]{strand}']
        rebuilt.append(SeqFeature(SimpleLocation(start, end, feat.location.strand),
                                  type=kind, qualifiers=quals))
    out.features = rebuilt + [f for f in out.features if f.type not in STRUCTURAL]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in', dest='src', required=True,
                    help='antiSMASH results directory (searched recursively)')
    ap.add_argument('--out', required=True, help='directory to write trimmed GBKs into')
    ap.add_argument('--kb', type=int, required=True,
                    help='neighbourhood in kb, each side of the rule core')
    args = ap.parse_args()

    src, dst = Path(args.src), Path(args.out)
    gbks = sorted(p for p in src.rglob('*.region*.gbk'))
    if not gbks:
        sys.exit(f'no *.region*.gbk under {src}')

    written = skipped = clipped = 0
    for path in gbks:
        record = SeqIO.read(path, 'genbank')
        before = len([f for f in record.features if f.type == 'CDS'])
        out = trim(record, args.kb)
        if out is None:
            skipped += 1
            continue
        after = len([f for f in out.features if f.type == 'CDS'])
        clipped += before - after
        # one directory per genome, mirroring the antiSMASH layout BiG-SCAPE expects
        target = dst / path.parent.name
        target.mkdir(parents=True, exist_ok=True)
        SeqIO.write(out, target / path.name, 'genbank')
        written += 1

    print(f'{written} regions re-cut to core +/- {args.kb} kb, {skipped} skipped '
          f'(no proto_core), {clipped} CDS dropped in total '
          f'({clipped / max(1, written):.2f} per region)')


if __name__ == '__main__':
    main()
