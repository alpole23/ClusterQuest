#!/usr/bin/env python3
"""Fail if any process script calls `sync`.

WHY THIS EXISTS
---------------
`sync` has killed two runs of this pipeline.

A bare `sync` flushes EVERY mounted filesystem, so one sick mount anywhere on the
host hangs a task that never touched it. On WSL2 the Windows-served 9p mounts do
exactly that: the call parks in uninterruptible D state, no signal reaches it, and
Nextflow is left with `IllegalThreadStateException: process hasn't exited`, which
aborts the RUN rather than the task. That killed a Bacteria GenBank run ~8,800 fetch
batches in, after the downloads had already succeeded.

`sync -f .` was adopted as the fix, documented as flushing only the filesystem
holding the work directory. It does not. On 2026-11-03 all four concurrent fetch
tasks wedged in `sync -f .` for 55 minutes with the work filesystem healthy and
writable, `Dirty` at 136 kB, `Writeback` at 0 and zero sectors moving, waiting on
`request_wait_answer` -- the FUSE/9p request path. SIGTERM moved them from S to D
state and they became unkillable.

Neither form is usable, and `timeout` cannot bound either one, because SIGTERM does
not reach a D-state process and `timeout` waits for a child that never exits. So the
rule is simply that no process script calls it, and this test is what keeps the
sentence in CLAUDE.md from drifting back into code.

To make a write durable, declare it as a process output and let Nextflow stage it.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# A shell command, not the word inside a comment or a longer identifier such as
# `rsync`, `syncfs` or `--sync`. Matches `sync` at the start of a command position:
# line start, or after a pipe, &&, ||, ; or `(`.
CALL = re.compile(r'(?:^|[|&;(]|&&|\|\|)[ \t]*sync\b(?![-\w])', re.M)
COMMENT = re.compile(r'^\s*(#|//)')

FAIL = []


def check(ok, label, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {label}' + (f' ({detail})' if detail else ''))
    if not ok:
        FAIL.append(label)


def offenders():
    hits = []
    for path in sorted(list(ROOT.glob('modules/**/*.nf'))
                       + list(ROOT.glob('subworkflows/**/*.nf'))
                       + [ROOT / 'main.nf']):
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if COMMENT.match(line):
                continue
            if CALL.search(line):
                hits.append(f'{path.relative_to(ROOT)}:{n}: {line.strip()}')
    return hits


def main():
    hits = offenders()
    check(not hits, 'no process script calls sync', '; '.join(hits[:4]))

    # The pattern has to actually fire, or this file is decoration. Two forms that
    # must be caught and three that must not.
    must_catch = ['    sync -f .', '    sync', '    cd x && sync -f .']
    must_pass = ['    rsync -a a b', '    # sync -f . was removed', '    syncfs_probe']
    caught = [s for s in must_catch if CALL.search(s)]
    passed = [s for s in must_pass if not CALL.search(s)]
    check(len(caught) == len(must_catch), 'the matcher catches every form of the call',
          f'missed {[s.strip() for s in must_catch if s not in caught]}')
    check(len(passed) == len(must_pass), 'and does not fire on rsync, syncfs or a comment',
          f'false positives {[s.strip() for s in must_pass if s not in passed]}')

    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
