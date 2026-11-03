#!/usr/bin/env python3
"""Check that run reports are per-run, and that everything reading them agrees.

WHY THIS EXISTS
---------------
trace/report/timeline are written with `overwrite = true`. The path was keyed by
TAXON, so any later run of the SAME taxon against the same outdir erased the
previous run's provenance -- and a `-preview`, which executes nothing, erased it
too. Four 114-byte header-only traces sat in results/pipeline_info/ as a result,
and on 2026-11-03 two previews destroyed the Erwiniaceae regression gate's trace,
report and timeline an hour after it finished. Its per-task figures were only
recoverable because each work dir keeps its own .command.trace.

`overwrite = false` is NOT the alternative. Measured on 26.10.5: TraceHelper
refuses an existing file, onFlowCreate catches it and logs a WARN, and the run
writes NO trace at all while still reporting [SUCCESS]. It changes which run you
lose, silently. So `overwrite = true` plus a per-run directory is the fix, and both
halves are checked here -- flipping either one reintroduces the loss.

The path is also read from four places outside the config, and a reader that
disagrees with the writer fails quietly: VISUALIZE_RESULTS simply omits its
resources panel. So the consumers are checked too.
"""
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FAIL = []


def check(ok, label, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {label}' + (f' ({detail})' if detail else ''))
    if not ok:
        FAIL.append(label)


def main():
    cfg = (ROOT / 'nextflow.config').read_text()

    # ── the writer ───────────────────────────────────────────────────────────
    stamp = re.search(r'params\.pipeline_info_dir\s*=(.{0,600}?)\n\n', cfg, re.S)
    decl = stamp.group(1) if stamp else ''
    check('pipeline_info' in decl and 'yyyyMMdd_HHmmss' in decl,
          'pipeline_info_dir carries a per-run stamp')
    check('ZoneOffset.UTC' in decl,
          'the stamp is UTC, so a DST fall-back cannot repeat an hour')
    # Date.format(String) lives in groovy-dateutil, not groovy-core, since 2.5.
    check('new java.util.Date()' not in decl,
          'the stamp does not depend on the groovy-dateutil extension')
    # Two evaluations would put trace and report in two different directories.
    check(decl.count('Instant.now()') == 1,
          'the stamp is evaluated exactly once', f"{decl.count('Instant.now()')} times")
    # The strict config DSL rejects a `def` here; nothing may reintroduce one.
    check(not re.search(r'^\s*def\s', decl, re.M),
          'no variable declaration, which 26.x strict config DSL rejects')

    for scope in ('trace', 'report', 'timeline'):
        block = re.search(scope + r'\s*\{(.*?)\n\}', cfg, re.S)
        body = block.group(1) if block else ''
        check('params.pipeline_info_dir' in body,
              f'{scope} writes under pipeline_info_dir, so it inherits the stamp')
        check('overwrite = true' in body,
              f'{scope} keeps overwrite = true '
              f'(false writes nothing and still reports SUCCESS)')

    # ── the readers ──────────────────────────────────────────────────────────
    bgc = (ROOT / 'subworkflows' / 'bgc_analysis.nf').read_text()
    check('${params.pipeline_info_dir}/pipeline_trace.tsv' in bgc,
          'bgc_analysis.nf reads the same value the config wrote',
          'a drift here silently drops the report resources panel')

    schema = json.loads((ROOT / 'nextflow_schema.json').read_text())
    prop = schema['properties'].get('pipeline_info_dir')
    # additionalProperties is false, so removing the property would make
    # validateParameters() reject the param on every run.
    check(prop is not None, 'pipeline_info_dir is still declared in the schema')
    check(prop is not None and 'default' not in prop,
          'the schema carries no default for it',
          'a stamped default would differ on every --write')
    check("VOLATILE_DEFAULT = {'pipeline_info_dir'}"
          in (ROOT / 'tests' / 'check_schema.py').read_text(),
          'check_schema.py --write will not bake a wall clock into the schema')

    # Both shell readers must match the pre-2026-11-03 layout as well: every outdir
    # written before then still holds it, and they have to read those too.
    status = (ROOT / 'scripts' / 'run_status.sh').read_text()
    check(status.count('pipeline_info/$TAXON/*/pipeline_trace.tsv') >= 1
          and status.count('pipeline_info/$TAXON/pipeline_trace.tsv') >= 1,
          'run_status.sh resolves both layouts')
    check('REL=${T#' in status and 'TAXON=${REL%%/*}' in status,
          'run_status.sh takes the taxon as the first segment under pipeline_info',
          'one dirname would return the run stamp instead')

    slurm = (ROOT / 'submit_slurm.sh').read_text()
    check('${OUTDIR}/pipeline_info/${TAXON_DIR}' in slurm,
          'submit_slurm.sh scopes its trace lookup to its own taxon',
          'an unscoped glob reports another taxon under this job header')
    check(slurm.count('TAXON_DIR=$(printf') == 1,
          'submit_slurm.sh derives the taxon directory once, not twice')

    arch = (ROOT / 'scripts' / 'bench_trace_archive.sh').read_text()
    check('GLOB=${1:?trace glob}' in arch,
          'the trace archiver takes a glob, not a path it cannot know in advance')
    check(re.search(r'while true; do.*?ls -1t \$GLOB', arch, re.S) is not None,
          'the archiver re-resolves the glob inside its loop')

    # ── the one piece of real logic, exercised ───────────────────────────────
    # A preview or an aborted run leaves a NEWER directory holding a header-only
    # trace. Selecting purely by mtime then reports "0 tasks" for a healthy run.
    fn = re.search(r'^newest_trace\(\) \{.*?^\}', status, re.S | re.M)
    if not fn:
        check(False, 'run_status.sh defines newest_trace()')
    else:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp) / 'pipeline_info' / 'Tax'
            real = base / '20261103_120000'
            empty = base / '20261103_130000'          # newer, header only
            real.mkdir(parents=True)
            empty.mkdir(parents=True)
            (real / 'pipeline_trace.tsv').write_text('hdr\nrow\n')
            (empty / 'pipeline_trace.tsv').write_text('hdr\n')
            import os
            import time
            now = time.time()
            os.utime(real / 'pipeline_trace.tsv', (now - 60, now - 60))
            os.utime(empty / 'pipeline_trace.tsv', (now, now))
            out = subprocess.run(
                ['bash', '-c', fn.group(0) + f'\nnewest_trace "{base}/*/pipeline_trace.tsv"'],
                capture_output=True, text=True).stdout.strip()
            check(out.endswith('20261103_120000/pipeline_trace.tsv'),
                  'newest_trace skips a newer header-only trace',
                  out or 'nothing returned')

    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
