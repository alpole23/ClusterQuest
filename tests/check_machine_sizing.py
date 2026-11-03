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
    #
    #    gtdbtk_gb is DELIBERATELY EXEMPT, and test 1b below is why. It is a
    #    measured requirement (~37 GiB, 2026-11-03), not a budget to be fitted:
    #    clamping it to a small machine is what let a 16 GB declaration tell
    #    Nextflow it could admit three 36.8 GiB shards at once. A machine that
    #    cannot afford it must be told to turn GTDB-Tk OFF, which is test 1b.
    bad = []
    for cores, ram in MACHINES:
        b = cm.budgets(cores, ram)
        for name in ('bigscape_gb', 'medium_gb'):
            if b[name] > b['usable_gb'] + 1e-9:
                bad.append(f'{cores}c/{ram}G {name}={b[name]:.1f} > usable {b["usable_gb"]:.1f}')
        if b['antismash_gb'] * 2 > b['usable_gb'] + 1e-9:
            bad.append(f'{cores}c/{ram}G antismash retry '
                       f'{b["antismash_gb"] * 2:.1f} > usable {b["usable_gb"]:.1f}')
        if b['usable_gb'] >= ram:
            bad.append(f'{cores}c/{ram}G leaves the OS nothing')
    check(not bad, 'no budget exceeds the machine it was sized for', '; '.join(bad[:3]))

    # 1b. A machine that cannot hold GTDB-Tk must be told to disable it, rather
    #     than handed a budget that fits and a task that will not be admitted.
    small = cm.assess(cm.detect(force_cores=8, force_ram_gb=32.0,
                                force_disk_gb=20000, outdir='.'),
                      dict(answers, run_gtdbtk=True), schema)
    prop = [d for d in small['degradations'] if d['param'] == 'run_gtdbtk']
    check(bool(prop), 'a 32 GB box is told to turn GTDB-Tk off',
          f'degradations: {[d["param"] for d in small["degradations"]]}')
    check(any('not be admitted' in d['why'] for d in prop),
          'and told why -- the task would not be admitted at all')

    big = cm.assess(cm.detect(force_cores=16, force_ram_gb=128.0,
                              force_disk_gb=20000, outdir='.'),
                    dict(answers, run_gtdbtk=True), schema)
    check(not [d for d in big['degradations'] if d['param'] == 'run_gtdbtk'],
          'a 128 GB box is not')

    # 1c. The allocation has to be declared where it cannot be clamped away. Its
    #     only label user is GTDBTK_CLASSIFY, but machine.config overrides labels,
    #     so the authoritative figure is a withName -- which beats any withLabel.
    labels = (ROOT / 'conf' / 'labels.config').read_text()
    import re
    blk = re.search(r"withName: 'GTDBTK_CLASSIFY' \{(.*?)\}", labels, re.S)
    check(blk is not None and "'40 GB'" in blk.group(1),
          'conf/labels.config declares GTDBTK_CLASSIFY memory by name',
          blk.group(1).strip() if blk else 'no withName block')
    check(f"{cm.GTDBTK_ALLOC_GB:.0f} GB" == '40 GB',
          'and check_machine.py agrees with it',
          f'{cm.GTDBTK_ALLOC_GB} GB')

    # 1d. Nothing may still be telling a user that GTDB-Tk needs 4 GB. It needs
    #     ~37 GiB; the 4 GB came from 60 s sampling that stepped over a 30 s spike.
    main_nf_txt = (ROOT / 'main.nf').read_text()
    check('sampled peak was 4 GB' not in main_nf_txt,
          'main.nf no longer claims a 4 GB GTDB-Tk peak')
    check(cm.GTDBTK_PEAK_GB > 30,
          'and the recorded peak is the measured one', f'{cm.GTDBTK_PEAK_GB} GiB')

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
    #    GTDB-Tk went 48 -> 16 on 2026-11-02 on the strength of 60 s sampling, and
    #    back to 40 on 2026-11-03 when 2-5 s sampling found the real ~37 GiB peak.
    #    So the original 48 was nearer right than the correction, and 40 is measured.
    b = cm.budgets(16, 56.9)
    expected = {'bigscape_gb': 16.0, 'medium_gb': 8.0, 'antismash_gb': 6.0,
                'antismash_forks': 7, 'executor_cpus': 32, 'gtdbtk_gb': 40.0}
    bad = [f'{k}={b[k]} want {v}' for k, v in expected.items() if b[k] != v]
    check(not bad, 'this box reproduces the pre-sizing values (GTDB-Tk 48->40 aside)',
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

    # ── the work directory and publish_mode = 'link' ────────────────────────
    #
    # `link` hard-links published output into outdir, and a hard link cannot cross
    # a filesystem. An outdir on a second drive with the work directory left at
    # Nextflow's default ./work therefore fails EVERY publish with `Invalid
    # cross-device link` -- AFTER the task has run, so the compute is paid for and
    # thrown away. Hit for real on 2026-11-02 with outdir on /mnt/bigdisk.
    #
    # Devices are stubbed rather than probed: a test that needed two real
    # filesystems would skip on any single-drive machine, which is most of them,
    # and this is the logic that has to hold everywhere.
    real_device = cm._device
    try:
        DEVS = {'/work/here': 1, '/drive2': 2}

        def fake_device(path):
            s = str(Path(path).absolute())
            for prefix, dev in sorted(DEVS.items(), key=lambda kv: -len(kv[0])):
                if s == prefix or s.startswith(prefix + '/'):
                    return dev
            return 1          # anything else shares the cwd's device

        cm._device = fake_device

        # Same filesystem: nothing to declare. Declaring a work directory that
        # Nextflow would have chosen anyway is noise in every single-drive config.
        got, why = cm.work_dir_for('/work/here/results')
        check(got is None and why is None,
              'same filesystem declares no workDir',
              '' if got is None else f'got {got!r}')

        # Different filesystem: a work directory beside outdir, and it must land on
        # OUTDIR's device -- that is the whole point. The sibling is preferred so
        # the work tree stays out of the published output.
        got, why = cm.work_dir_for('/drive2/clusterquest/results')
        ok = got is not None and fake_device(got) == fake_device('/drive2/clusterquest/results')
        check(ok, 'cross-filesystem outdir gets a workDir on its own device',
              '' if ok else f'got {got!r}')
        sibling = got == Path('/drive2/clusterquest/work')
        check(sibling, 'the work directory is a sibling of outdir, not inside it',
              '' if sibling else f'got {got!r}')
        reasoned = bool(why) and 'cross devices' in (why or '')
        check(reasoned, 'the reason names the hard-link limit',
              '' if reasoned else f'why={why!r}')

        # outdir AT the mount point: its parent is on the other device, so the
        # sibling would be wrong and it has to fall back to inside outdir. This is
        # the case a naive `outdir.parent / "work"` gets silently wrong.
        got, _ = cm.work_dir_for('/drive2')
        inside = got == Path('/drive2/work')
        check(inside, 'outdir at a mount point falls back to a work dir inside it',
              '' if inside else f'got {got!r}')

        # And the rendered config has to actually carry it.
        res = cm.assess(cm.detect(force_cores=16, force_ram_gb=56.9,
                                  force_disk_gb=12000,
                                  outdir='/drive2/clusterquest/results'),
                        dict(answers, outdir='/drive2/clusterquest/results'), schema)
        rendered = cm.render_config(res)
        declared = "workDir = '/drive2/clusterquest/work'" in rendered
        check(declared, 'the generated config declares the work directory',
              '' if declared else 'absent from render_config output')

        res_same = cm.assess(cm.detect(force_cores=16, force_ram_gb=56.9,
                                       force_disk_gb=12000,
                                       outdir='/work/here/results'),
                             dict(answers, outdir='/work/here/results'), schema)
        check('workDir' not in cm.render_config(res_same),
              'a single-filesystem machine gets no workDir line')
    finally:
        cm._device = real_device

    # main.nf must refuse this before the first task rather than after it. The
    # check lives in validateParams() and reads params.publish_mode, so a grep is
    # the cheap guard against it being dropped; TEST 1/2 above exercise it live.
    main_nf = (ROOT / 'main.nf').read_text()
    check("sameFileStore(params.outdir, workDir)" in main_nf
          and 'publish_mode' in main_nf,
          'main.nf refuses a cross-device link before running anything')

    # The sizing notice has to fire on BOTH branches. It spoke up only when
    # machine.config was absent until 2026-11-03, which left a sized run silent --
    # so the only way to know sizing had happened was to go back and read the
    # wizard's scrollback. A notice that says nothing on success is not a notice.
    notice = main_nf[main_nf.find('def sized = workflow.configFiles'):]
    notice = notice[:notice.find('// Warnings')]
    check('Resources sized to this machine by' in notice,
          'main.nf says so when the run WAS sized')
    check('portable floor' in notice and 'check_machine.py' in notice,
          'main.nf says so when it was NOT, and how to fix it')

    # It prints the header of the file that did the sizing, which is where the
    # measurements and the generation date live -- so a config carried over from
    # another machine is visible rather than silently in force. Reading the resolved
    # config instead would print an uncalled closure for any dynamic memory
    # directive, and would lean on workflow.session, which is not documented API.
    check('Detected:' in notice and 'Generated by' in notice,
          'the notice prints the sizing header, so a stale config is visible')
    check('workflow.session' not in notice,
          'the notice does not read the resolved config to print a number')

    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
