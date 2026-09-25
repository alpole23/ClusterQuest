#!/usr/bin/env python3
"""The pepM screen is invoked from two places; assert they pass the same flags.

`PEPM_PRESCREEN` screens genomes supplied through `--input_genomes`, and
`FETCH_RENAME_SCREEN` screens genomes it has just downloaded. Both shell out to
the same `scripts/analysis/pepm_prescreen.py`, so the screening *algorithm*
cannot drift between them. The *arguments* can, and did.

The fused call site was written without `--threads` and `--diamond`. Nothing
failed, because the two happened to agree: `process_medium` is 4 CPUs and the
script's `--threads` default is 4, and diamond was on PATH from the conda
environment. Raise that label to 8 and one path would silently use half the
threads of the other, with no error and no visible difference in output --
exactly the failure mode this repository keeps finding.

Two invocations of one decision need to be checked BY CONSTRUCTION rather than
by remembering. This compares the flag sets and fails on any difference.

Values are deliberately NOT compared: the two legitimately differ in how they
name their inputs (`*.gbff` after a download versus a staged `${genomes}`), and
in `task.index` for the report name. What must match is WHICH knobs are set.

    python tests/check_screen_flags.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CALL_SITES = [
    ROOT / 'modules' / 'analysis' / 'pepm_prescreen.nf',
    ROOT / 'modules' / 'genome' / 'fetch_rename_screen.nf',
]
SCRIPT = 'analysis/pepm_prescreen.py'


def flags_of(path):
    """The set of --flags passed to pepm_prescreen.py in this module, or None."""
    text = path.read_text()
    idx = text.find(SCRIPT)
    if idx < 0:
        return None
    # The invocation runs until a line that does not end in a continuation.
    tail, out = text[idx:], []
    for line in tail.splitlines():
        out.append(line)
        if not line.rstrip().endswith('\\'):
            break
    return set(re.findall(r'--([a-z_]+)', '\n'.join(out)))


def main():
    found = {p.name: flags_of(p) for p in CALL_SITES}
    missing = [n for n, f in found.items() if f is None]
    if missing:
        print(f'pepm_prescreen.py is not invoked in: {", ".join(missing)}. '
              f'If a call site moved, update CALL_SITES in this test.')
        return 1

    names = list(found)
    a, b = found[names[0]], found[names[1]]
    if a == b:
        print(f'pepM screen flags agree across {len(names)} call sites: '
              f'{" ".join("--" + f for f in sorted(a))}')
        return 0

    print('pepM screen call sites pass different flags:')
    for name, only in ((names[0], a - b), (names[1], b - a)):
        if only:
            print(f'  only in {name}: {" ".join("--" + f for f in sorted(only))}')
    print('\nBoth run the same script, so this does not change the algorithm — it '
          'changes how it is configured, silently. Add the missing flag, or if the '
          'difference is deliberate, record why here.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
