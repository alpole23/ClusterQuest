#!/usr/bin/env python3
"""Check that machine sizing cannot produce a configuration the machine refuses.

WHY THIS NEEDS A TEST AT ALL
----------------------------
Nextflow's local executor does not run an oversized task slowly -- it refuses to
admit it, and the run dies:

    Process requirement exceeds available memory -- req: 48 GB; avail: 12 GB
    [FAILED] completed=0 failed=0 cached=0

So an arithmetic slip in a budget does not degrade performance, it kills a run
that may be two days in. Every check here exists because the failure it guards is
silent until the task is submitted.

Three things can drift independently and all three are checked:

  * a budget larger than the machine it was computed for -- the failure above;
  * nextflow.config's derived `executor.cpus` / `maxForks` disagreeing with what
    check_machine.py computes for the same box, so the floor and the override
    describe different machines;
  * the BiG-SCAPE partition cap drifting from the memory fit it is supposed to
    invert. That one was a REAL bug: nextflow.config documented the tree-PSS fit
    GB = 1.80 + 2.04e-7*n^2 while partition_bgcs.py still inverted the superseded
    ru_maxrss fit, returning a cap ~27% too high.

Run by tests/run_tests.sh.
"""
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))

FAIL = []

# The grid of machines sizing has to survive. The first is this workstation, the
# last two are below the declared floor and must be refused rather than clamped.
MACHINES = [
    (16, 56.9), (16, 64.0), (8, 32.0), (8, 16.0), (4, 16.0),
    (4, 8.0), (2, 16.0),
]


def check(ok, label, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {label}' + (f' ({detail})' if detail else ''))
    if not ok:
        FAIL.append(label)


def main():
    import check_machine as cm
    import configure_run as cr

    schema = cr.load_schema()
    answers = {k: v.get('default') for k, v in schema.items()}

    # 1. No budget may exceed the machine. antiSMASH is checked at DOUBLE its
    #    request, because it escalates once on an OOM kill -- `memory = 6.GB *
    #    task.attempt` -- and the retry is admitted against the same machine.
    bad = []
    for cores, ram in MACHINES:
        b = cm.budgets(cores, ram)
        for name in ('gtdbtk_gb', 'bigscape_gb', 'medium_gb'):
            if b[name] > b['usable_gb'] + 1e-9:
                bad.append(f'{cores}c/{ram}G {name}={b[name]:.1f} > usable {b["usable_gb"]:.1f}')
        if b['antismash_gb'] * 2 > b['usable_gb'] + 1e-9:
            bad.append(f'{cores}c/{ram}G antismash retry '
                       f'{b["antismash_gb"] * 2:.1f} > usable {b["usable_gb"]:.1f}')
        if b['usable_gb'] >= ram:
            bad.append(f'{cores}c/{ram}G leaves the OS nothing')
    check(not bad, 'no budget exceeds the machine it was sized for', '; '.join(bad[:3]))

    # 2. The widest cpus request must be admissible. `process_high` and
    #    `process_high_memory` ask for 8, and the profile oversubscribes CPU 2x,
    #    which is exactly why the core floor is 4 and not 2.
    widest = 8
    bad = [f'{c}c offers {c * 2}' for c, _r in MACHINES
           if c >= cm.MIN_CORES and c * 2 < widest]
    check(not bad, f'every supported core count admits a {widest}-cpu task',
          '; '.join(bad))

    # 3. This workstation's budgets must reproduce the absolutes they replaced, or
    #    the Enterobacterales re-run is not comparable to the committed trace.
    #    GTDB-Tk is the one DELIBERATE change: 48 -> 16 GB, because 48 was never a
    #    measurement and the sampled peak across its whole 3.6 h window was 4 GB.
    b = cm.budgets(16, 56.9)
    expected = {'bigscape_gb': 16.0, 'medium_gb': 8.0, 'antismash_gb': 6.0,
                'antismash_forks': 7, 'executor_cpus': 32, 'gtdbtk_gb': 16.0}
    bad = [f'{k}={b[k]} want {v}' for k, v in expected.items() if b[k] != v]
    check(not bad, 'this box reproduces the pre-sizing values (GTDB-Tk 48->16 aside)',
          '; '.join(bad))

    # 4. The floor refuses rather than clamps. Below it the arithmetic still works
    #    -- an 8 GB box gets antiSMASH 3 GB, 6 on retry -- so nothing but this
    #    check stops an untested configuration from launching.
    for cores, ram, why in ((4, 8.0, 'RAM'), (2, 16.0, 'cores')):
        machine = dict(cores=cores, ram_gb=ram, disk_free_gb=5000.0,
                       disk_probe='t', platform='t', notes=[])
        got = cm.assess(machine, answers, schema)['verdict']
        check(got == 'refused', f'{cores}c/{ram:g}G is refused ({why} below floor)',
              '' if got == 'refused' else f'verdict {got!r}')

    # 5. A machine at the floor with plenty of disk must be usable, not refused.
    machine = dict(cores=cm.MIN_CORES, ram_gb=cm.MIN_RAM_GB, disk_free_gb=5000.0,
                   disk_probe='t', platform='t', notes=[])
    got = cm.assess(machine, answers, schema)
    check(got['verdict'] != 'refused',
          f'{cm.MIN_CORES}c/{cm.MIN_RAM_GB:g}G at the floor is supported',
          f'verdict {got["verdict"]!r}: {"; ".join(got["blockers"])}')

    # 6. Disk, not RAM, is what gates GTDB-Tk -- so a floor machine with too little
    #    disk must degrade run_gtdbtk rather than be refused outright.
    machine = dict(cores=4, ram_gb=16.0, disk_free_gb=60.0,
                   disk_probe='t', platform='t', notes=[])
    got = cm.assess(machine, answers, schema)
    degraded = [d['param'] for d in got['degradations'] if not d.get('advisory')]
    check(got['verdict'] == 'degraded' and 'run_gtdbtk' in degraded,
          'too little disk degrades run_gtdbtk instead of refusing',
          f'verdict {got["verdict"]!r}, degradations {degraded}')

    # 7. partition_bgcs.py must invert the SAME fit nextflow.config documents. The
    #    two disagreed until 2026-11-02 and the cap came out ~27% high; that is
    #    survivable while the cap is set by hand and not once memory is sized from
    #    the machine, because then this is the only thing between a small box and
    #    an OOM that discards the clustering.
    part = (ROOT / 'scripts' / 'clustering' / 'partition_bgcs.py').read_text()
    m = re.search(r'cap = int\(\(\(budget - ([\d.]+)\) / ([\d.e-]+)\) \*\* 0\.5\)', part)
    if not m:
        check(False, 'partition_bgcs.py cap derivation found')
    else:
        got_fit = (float(m.group(1)), float(m.group(2)))
        want = (cm.BIGSCAPE_INTERCEPT_GB, cm.BIGSCAPE_COEFF)
        check(got_fit == want, 'partition cap inverts the tree-PSS fit',
              '' if got_fit == want else f'{got_fit} != {want}')

    # 8. And that fit must still reproduce CLAUDE.md's published capacity table,
    #    which is what the manuscript and fig4 quote.
    published = {16: 8350, 48: 15060, 64: 17470, 128: 24880}
    bad = []
    for gb, want in published.items():
        got_cap = cm.bigscape_capacity(gb)
        if abs(got_cap - want) > max(20, want * 0.01):
            bad.append(f'{gb} GB -> {got_cap:,} want ~{want:,}')
    check(not bad, 'capacity matches the published table', '; '.join(bad))

    # 9. The floor baked into nextflow.config's `local` profile must itself fit the
    #    declared floor machine, since that profile is what a bare `nextflow run`
    #    gets with no machine.config at all.
    cfg = (ROOT / 'nextflow.config').read_text()
    floor_usable = cm.budgets(cm.MIN_CORES, cm.MIN_RAM_GB)['usable_gb']
    bad = []
    for label in ('process_high_memory', 'process_high', 'process_medium'):
        m = re.search(rf"withLabel: '{label}' \{{\s*\n?\s*memory = '(\d+) GB'", cfg)
        if not m:
            bad.append(f'{label}: not found in the local profile')
        elif int(m.group(1)) > floor_usable:
            bad.append(f'{label}={m.group(1)} GB > {floor_usable:.1f} usable at the floor')
    check(not bad, 'the local profile floor fits the floor machine', '; '.join(bad))

    # 10. nextflow.config derives executor.cpus and ANTISMASH maxForks in place,
    #     because the strict config DSL allows bare expressions. check_machine.py
    #     computes the same two numbers. They must agree, or the portable floor and
    #     the generated override describe different machines.
    import os
    here = os.cpu_count() or 1
    b_here = cm.budgets(here, 999.0)
    for pat, key, what in (
            (r'executor\.cpus = Runtime\.runtime\.availableProcessors\(\) \* (\d+)',
             'executor_cpus', 'executor.cpus'),
            (r'maxForks = Math\.max\(1, \(Runtime\.runtime\.availableProcessors\(\) - (\d+)\)\.intdiv\((\d+)\)\)',
             'antismash_forks', 'ANTISMASH maxForks')):
        m = re.search(pat, cfg)
        if not m:
            check(False, f'{what} is derived in nextflow.config')
            continue
        if key == 'executor_cpus':
            got_val = here * int(m.group(1))
        else:
            got_val = max(1, (here - int(m.group(1))) // int(m.group(2)))
        check(got_val == b_here[key], f'{what} agrees with check_machine.py',
              '' if got_val == b_here[key] else f'config {got_val} vs script {b_here[key]}')

    # 10b. The cgroup walk must find this process's OWN cgroup, not just the root.
    #      Reading /sys/fs/cgroup/memory.max alone happens to work under Docker,
    #      where the container's cgroup is its namespace root, and silently fails
    #      for a systemd scope on the host -- measured 2026-11-02, a
    #      `systemd-run --scope -p MemoryMax=16G` was detected as the host's
    #      56.9 GB. Under a real limit the walk is also exercised end to end by
    #      the systemd-run case in the sizing validation; here, assert the shape.
    paths = cm._cgroup_paths()
    root = Path('/sys/fs/cgroup')
    if not root.is_dir():
        print('  SKIP  cgroup walk reaches this process\'s own cgroup (no /sys/fs/cgroup)')
    else:
        ok = bool(paths) and paths[-1] == root and all(
            p == root or root in p.parents for p in paths)
        check(ok, "cgroup walk runs from this process's cgroup up to the root",
              '' if ok else f'{[str(p) for p in paths[:4]]}')
        # A limit found anywhere in the chain must be the minimum, not the first hit.
        lim = cm._cgroup_ram_gb()
        check(lim is None or lim > 0, 'cgroup memory limit is positive when present',
              '' if (lim is None or lim > 0) else f'{lim}')

    # 11. The generated config must actually parse and resolve to the budgets it
    #     claims. Everything above is arithmetic about Nextflow; this is Nextflow.
    #     Skipped without a JVM, like the rest of the suite's nextflow checks.
    nf = None
    for cand in ('nextflow',):
        from shutil import which
        nf = which(cand)
    if not nf:
        print('  SKIP  generated config resolves under Nextflow (nextflow not on PATH)')
    else:
        import tempfile
        machine = dict(cores=4, ram_gb=16.0, disk_free_gb=5000.0,
                       disk_probe='t', platform='t', notes=[])
        res = cm.assess(machine, answers, schema)
        with tempfile.TemporaryDirectory() as td:
            cfg_path = Path(td) / 'machine.config'
            cfg_path.write_text(cm.render_config(res))
            out = subprocess.run([nf, '-c', str(cfg_path), 'config',
                                  '-profile', 'local', str(ROOT)],
                                 capture_output=True, text=True, timeout=600)
            if out.returncode != 0:
                check(False, 'generated config resolves under Nextflow',
                      (out.stderr or out.stdout).strip().splitlines()[0][:120]
                      if (out.stderr or out.stdout).strip() else 'non-zero exit')
            else:
                want_gtdb = f"memory = '{res['budgets']['gtdbtk_gb']:.0f} GB'"
                want_cpus = f"cpus = {res['budgets']['executor_cpus']}"
                missing = [w for w in (want_gtdb, want_cpus) if w not in out.stdout]
                check(not missing, 'generated config resolves to its own budgets',
                      f'missing {missing}' if missing else '')

    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
