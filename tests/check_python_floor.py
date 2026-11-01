#!/usr/bin/env python3
"""Compile every pipeline script with the OLDEST Python any process runs.

The suite used to compile with whatever `python3` the developer has. That is the
newest interpreter on the box, and it accepts syntax the pipeline's conda
environments reject. A backslash inside an f-string expression is legal from
Python 3.12 (PEP 701) and a SyntaxError before it; the antiSMASH environment runs
3.11, so `viz/gcf_network.py` compiled clean locally, passed the suite, merged, and
killed VISUALIZE_RESULTS on a real Erwiniaceae run after every expensive task had
already succeeded. `viz/cutoff_comparison.py` carried the same error and had never
been executed at all.

The floor is discovered from the installed conda environments rather than pinned, so
it tracks reality instead of a number someone has to remember to update. With no
environments present the check SKIPS rather than passing: a silent pass here is what
let the bug through.
"""
import pathlib
import py_compile
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parent.parent
CONDA = [pathlib.Path.home() / '.clusterquest' / 'conda',
         REPO / 'work' / 'conda']


def interpreters():
    """{(major, minor): path} for every python in an installed env."""
    found = {}
    for root in CONDA:
        if not root.is_dir():
            continue
        for py in sorted(root.glob('env-*/bin/python3.*')):
            if py.is_symlink() or not py.is_file():
                continue
            try:
                v = subprocess.run([str(py), '-c',
                                    'import sys;print(sys.version_info[:2])'],
                                   capture_output=True, text=True, timeout=30)
                ver = eval(v.stdout.strip())          # "(3, 11)"
            except Exception:
                continue
            found.setdefault(ver, py)
    return found


def main():
    found = interpreters()
    if not found:
        print('SKIP  no conda environments found; cannot determine the Python floor')
        return 0
    floor = min(found)
    py = found[floor]
    print(f'      floor is Python {floor[0]}.{floor[1]} ({py})')

    scripts = sorted(REPO.glob('scripts/**/*.py'))
    bad = []
    for s in scripts:
        r = subprocess.run([str(py), '-c',
                            f'import py_compile,sys;'
                            f'py_compile.compile({str(s)!r}, cfile=sys.argv[1],'
                            f' doraise=True)',
                            tempfile.mkstemp(suffix='.pyc')[1]],
                           capture_output=True, text=True)
        if r.returncode:
            msg = (r.stderr.strip().splitlines() or ['?'])[-1]
            bad.append((s.relative_to(REPO), msg))

    if bad:
        for path, msg in bad:
            print(f'FAIL  {path}: {msg}')
        return 1
    print(f'      {len(scripts)} scripts compile at the floor')
    return 0


if __name__ == '__main__':
    sys.exit(main())
