#!/usr/bin/env python3
"""Set a detection rule's NEIGHBOURHOOD in the installed antiSMASH rule file.

antiSMASH sizes a region as its rule core plus a fixed neighbourhood, and the strict
`phosphonate` rule declares NEIGHBOURHOOD 5 (kb). That is too small for this chemistry:
on *Pantoea ananatis* LMG 5342 the region stops 4.2 kb short of the HiVir cluster, so
the MFS transporter, the hypothetical, the FMN reductase and the second ATP-grasp fall
outside it — in all 215 members of that family. Measured on one genome:

    NEIGHBOURHOOD  5   796,909-810,246   13,338 bp    8,337 of 12,526 bp of cluster
    NEIGHBOURHOOD 10   791,909-815,246   23,338 bp    ALL of it, +0.8 kb past its end
    NEIGHBOURHOOD 20   781,909-825,246   43,338 bp    ALL of it, +10.8 kb past its end

`--hmmdetection-strictness relaxed` does NOT help: `phosphonate-like`, the relaxed rule
that carries NEIGHBOURHOOD 20, never fires on this cluster — verified, no trace of it in
the output — so the strict rule and its 5 kb are what apply either way. antiSMASH has no
command-line option for a bacterial neighbourhood (only fungal multipliers), and no
option to point at a different rule file, so the installed file is the only lever.

The neighbourhood is symmetric, so raising it also pulls in the same amount of upstream
DNA the cluster does not contain. 10 is the tuned value for HiVir, not a general one.

This edits a file inside the conda environment Nextflow manages, which is shared by
every concurrent ANTISMASH task. It is therefore idempotent (a file already at the
requested value is left alone) and atomic (written to a temporary file in the same
directory, then os.replace), so a task reading the rules always sees one whole version.
Each task patches before invoking antiSMASH, so a task's own run always sees the value
it asked for. The environment is recreated whenever the conda spec changes or the cache
is cleared, which is why this runs per task rather than once.
"""
import argparse
import os
import re
import sys
import tempfile
from pathlib import Path


def default_rule_file():
    import antismash
    return (Path(antismash.__file__).parent / 'detection' / 'hmm_detection' /
            'cluster_rules' / 'strict.txt')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--neighbourhood', type=int, required=True,
                    help='new NEIGHBOURHOOD value, in kb')
    ap.add_argument('--rule', default='phosphonate', help='rule name to edit')
    ap.add_argument('--file', type=Path, default=None,
                    help='rule file (default: strict.txt inside the installed antiSMASH)')
    a = ap.parse_args()

    path = a.file or default_rule_file()
    if not path.is_file():
        sys.exit(f'{path}: antiSMASH rule file not found')

    text = path.read_text()
    block = re.search(rf'^RULE {re.escape(a.rule)}\n.*?(?=^RULE |\Z)', text, re.S | re.M)
    if not block:
        sys.exit(f'{path}: no "RULE {a.rule}"')
    current = re.search(r'^(\s*)NEIGHBOURHOOD (\d+)', block.group(0), re.M)
    if not current:
        sys.exit(f'{path}: rule "{a.rule}" declares no NEIGHBOURHOOD')

    if int(current.group(2)) == a.neighbourhood:
        print(f'{a.rule}: NEIGHBOURHOOD already {a.neighbourhood}')
        return 0

    patched = re.sub(r'^(\s*)NEIGHBOURHOOD \d+', rf'\g<1>NEIGHBOURHOOD {a.neighbourhood}',
                     block.group(0), count=1, flags=re.M)
    new_text = text[:block.start()] + patched + text[block.end():]

    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix='.rules-')
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write(new_text)
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)          # atomic for concurrent readers
    except BaseException:
        os.path.exists(tmp) and os.unlink(tmp)
        raise

    check = re.search(rf'^RULE {re.escape(a.rule)}\n.*?(?=^RULE |\Z)',
                      path.read_text(), re.S | re.M)
    got = re.search(r'^\s*NEIGHBOURHOOD (\d+)', check.group(0), re.M)
    if not got or int(got.group(1)) != a.neighbourhood:
        sys.exit(f'{path}: NEIGHBOURHOOD did not take; regions would be sized wrongly')
    print(f'{a.rule}: NEIGHBOURHOOD {current.group(2)} -> {a.neighbourhood} in {path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
