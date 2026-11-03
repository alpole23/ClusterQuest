#!/usr/bin/env python3
"""Size a run to the machine it will actually run on, and refuse what will not fit.

WHY THIS EXISTS
---------------
Nextflow's local executor does not run an oversized task slowly -- it refuses to
ADMIT it, and the run dies:

    Process requirement exceeds available memory -- req: 48 GB; avail: 12 GB
    [FAILED] completed=0 failed=0 cached=0

So every absolute resource number in a profile is a hidden minimum machine spec.
`-profile local` used to request 48 GB for GTDB-Tk and therefore silently required
a >=48 GB box.

That 48 GB was cut to 4 GB on 2026-11-02, from the RefSeq Enterobacterales run's
60-SECOND sampling (2,608 samples; 206 across the 3.6 h GTDB-Tk window), which
showed mem_used peaking at 4 GB. On 2026-11-03 that was measured again at 5 s and
2 s and the peak is ~37 GiB in a ~30-SECOND spike -- the 60 s sampling had simply
stepped over it. Three independent instruments agree: PSS 36.8 GiB, system
mem_used 38.8 GiB, Nextflow peak_rss 36.9 GB. So the 48 GB this script was written
to mock was nearer the truth than the 4 GB that replaced it, and GTDB-Tk really
does require a big box. It is now allocated 40 GB BY NAME in conf/labels.config.

GTDB-Tk also consumes 150 GB of WORK DISK, which is unchanged. Both gate it now.

This script therefore does two different jobs, and it is worth keeping them apart:

  * RAM and cores decide the RESOURCE REQUESTS, and they are generous because the
    measured peaks are small.
  * DISK decides WHICH STAGES CAN RUN AT ALL, and it is the real constraint.

Sizing lives here rather than in nextflow.config because Nextflow 26's config
parser is a strict DSL. Verified on 26.04.3, 2026-11-02:

  * `def` is rejected anywhere in nextflow.config, including line 1
    ("Variable declarations cannot be mixed with config statements");
  * statements are rejected inside config closures ("Unexpected input: 'if'");
  * a C-style `(int)` cast is rejected ("Unable to parse config file", with no
    line number, which is a long debugging session if you do not know this);
  * lib/*.groovy is NOT on the classpath while the config is parsed
    ("Sizing is not defined"), so a helper class next to lib/Utils.groovy
    cannot be called from there either.

Only bare expressions survive, which is why nextflow.config can derive
`executor.cpus` and `maxForks` in place but not the memory budgets.

WHAT IT WRITES
--------------
A `machine.config` that overrides the `local` profile. A `-c` file beats the
profile regardless of profile order -- verified 2026-11-02 that LOAD order, not
profile order, decides -- so overriding `local` by name is enough and the run
command keeps a single `-profile local`.

    python scripts/check_machine.py
    nextflow run main.nf -profile local -c machine.config -params-file run.json

`scripts/configure_run.py` calls this at the end of the interview, so the normal
path needs none of the above by hand.

SLURM is deliberately untouched: the cluster profile's numbers are declared
capacity, not a property of the submit host, and sizing them from whatever node
you happened to log in to would be wrong. `--profile slurm` short-circuits.
"""
import argparse
import json
import math
import os
import platform
import shutil
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))

GIB = 2 ** 30

# ── The floor ────────────────────────────────────────────────────────────────
# 16 GB / 4 cores. Not arbitrary: 4 cores is where the pipeline's widest request
# still fits. `process_high` and `process_high_memory` ask for 8 cpus, and the
# local profile oversubscribes CPU 2x, so a 4-core box offers 8 slots and admits
# them exactly. At 2 cores the executor offers 4, an 8-cpu task is refused, and
# the run dies the same way the 48 GB request used to.
MIN_CORES = 4
MIN_RAM_GB = 16.0

# ── Measured constants ───────────────────────────────────────────────────────
# BiG-SCAPE peak is quadratic above ~4,000 BGCs. TWO fits exist and you must say
# which: GB = 1.14 + 1.29e-7*n^2 is on ru_maxrss, which reports the max of any ONE
# process, and BiG-SCAPE runs eleven. Tree-wide PSS -- which divides each shared
# page by its sharers, so the shared Pfam data is counted once -- is 1.58x higher
# (stable 1.51-1.65 over five sizes). What the machine has to hold is PSS.
BIGSCAPE_INTERCEPT_GB = 1.80
BIGSCAPE_COEFF = 2.04e-7

# GTDB-Tk. Memory measured 2026-11-03 at 2-5 s; the earlier 4.0 came from 60 s
# sampling that stepped over a ~30 s spike. Standalone at three sizes on the same
# genomes, flags and database: 25 -> 10.67 GiB, 100 -> 34.44, 300 -> 36.22. It
# SATURATES, so a smaller gtdbtk_shard_size buys no memory above ~100 genomes while
# wall time stays fixed per shard at ~27 min (the ANI pre-screen sketches all
# 143,614 reference genomes every time). Caveat: all measured genomes were
# Erwiniaceae, one family; a broader taxon may place higher.
GTDBTK_PEAK_GB = 36.8             # measured PSS peak, 301 genomes
GTDBTK_ALLOC_GB = 40.0            # what conf/labels.config declares, by name
GTDBTK_WORK_SCRATCH_GB = 150.0    # disk, unchanged
# Peak work disk away from the GTDB-Tk window, same run: 22 GB during the
# antiSMASH/ORF phase. That is NOT an input-scaled figure -- fetch batches are
# deleted as they are screened, so what is resident is whatever antiSMASH and
# RECOVER_ORFS hold open at once, which is bounded by maxForks. 22 GB at the 7
# concurrent batches this box allows is ~3.1 GB per fork, so a 4-core laptop
# running 1 fork needs a few GB here, not 22.
WORK_SCRATCH_PER_FORK_GB = 3.1
WORK_SCRATCH_FLOOR_GB = 5.0
# Output directory at order scale: 30 GB for 1,309 retained of 150,690 screened.
# It scales with RETAINED genomes, and retention is clade-specific and not knowable
# before the screen runs -- 0.113% on GenBank Enterobacterales against a measured
# 7.9% false-positive rate on actinomycetes, a 70x spread. So this is the measured
# order-scale figure held as an ASSUMPTION rather than a projection, and the report
# says so. `--outdir-gb` overrides it for a run that is known to be small.
OUTDIR_ORDER_SCALE_GB = 30.0

# The largest BGC count a run of these taxa plausibly produces, used to decide
# whether the BiG-SCAPE partition cap is worth warning about. GenBank
# Enterobacterales yielded 1,476 from 1.35 M genomes; Actinomycetes projects to
# 1,000-4,000 at the measured 7.9% actinomycete screen retention.
PLAUSIBLE_BGC_YIELD = 5000


def bigscape_capacity(mem_gb):
    """Largest single BiG-SCAPE job that fits `mem_gb`, in BGCs.

    The raw inversion of the PSS fit, with no margin: this is the figure CLAUDE.md
    publishes (8,350 at 16 GB). `bigscape_cap` below is deliberately lower.
    """
    budget = mem_gb - BIGSCAPE_INTERCEPT_GB
    return 0 if budget <= 0 else int(math.sqrt(budget / BIGSCAPE_COEFF))


def bigscape_cap(mem_gb, margin=0.85):
    """The PARTITION cap for `mem_gb`, in BGCs -- capacity times the margin.

    Inverts the PSS fit. `margin` is budgeted below 1.0 because the fit is being
    extrapolated past its data (measured to 10,000 BGCs), and the asymmetry is
    stark: the margin costs ~8% of the cap, an OOM kill discards hours of
    clustering. Returns 0 when the budget cannot even cover the intercept.
    """
    budget = mem_gb * margin - BIGSCAPE_INTERCEPT_GB
    return 0 if budget <= 0 else int(math.sqrt(budget / BIGSCAPE_COEFF))


def _read_int(path):
    try:
        text = Path(path).read_text().strip()
    except OSError:
        return None
    return int(text) if text.lstrip('-').isdigit() else None


def _cgroup_paths():
    """Every cgroup directory from this process's own cgroup up to the root.

    The effective limit is the MINIMUM along that chain, and it is not necessarily
    at the root. Reading only /sys/fs/cgroup/<file> happens to work under Docker,
    where the container's cgroup IS the root of its namespace -- but not for a
    systemd scope or slice on the host, where the limit sits at e.g.
    /sys/fs/cgroup/user.slice/user-1000.slice/...scope/memory.max while the root
    reads "max". Measured 2026-11-02: a `systemd-run --scope -p MemoryMax=16G` was
    detected as the host's 56.9 GB until this walked the chain.
    """
    root = Path('/sys/fs/cgroup')
    if not root.is_dir():
        return []
    rel = ''
    try:
        for line in Path('/proc/self/cgroup').read_text().splitlines():
            parts = line.split(':', 2)
            # cgroup v2 is the "0::<path>" line; v1 names a controller.
            if len(parts) == 3 and (parts[0] == '0' or 'memory' in parts[1]):
                rel = parts[2].lstrip('/')
                if parts[0] == '0':
                    break
    except OSError:
        pass

    paths, here = [], (root / rel if rel else root)
    while True:
        paths.append(here)
        if here == root or root not in here.parents:
            break
        here = here.parent
    return paths


def _cgroup_ram_gb():
    """Container or scope memory limit, or None.

    Checked because /proc/meminfo reports the HOST under a container, so a 4 GB
    container on a 64 GB host would otherwise be sized as a 64 GB machine and then
    OOM-killed. cgroup v1 writes a ~8 EB sentinel rather than "max" when unlimited,
    which the upper bound rejects.
    """
    limits = []
    for d in _cgroup_paths():
        for name in ('memory.max', 'memory.limit_in_bytes'):
            raw = _read_int(d / name)
            if raw and 0 < raw < (1 << 50):
                limits.append(raw / GIB)
    return min(limits) if limits else None


def _cgroup_cores():
    """Container or scope CPU quota in whole cores, or None.

    cpu.max is "<quota> <period>". Walked up the chain for the same reason as
    memory: a systemd scope's CPUQuota does not live at the cgroup root.
    """
    quotas = []
    for d in _cgroup_paths():
        try:
            fields = (d / 'cpu.max').read_text().split()
        except OSError:
            continue
        if len(fields) != 2 or fields[0] == 'max':
            continue
        try:
            quota, period = int(fields[0]), int(fields[1])
        except ValueError:
            continue
        if period > 0:
            quotas.append(max(1, quota // period))
    return min(quotas) if quotas else None


def _host_ram_gb():
    kb = None
    try:
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemTotal:'):
                kb = int(''.join(c for c in line if c.isdigit()))
                break
    except OSError:
        pass
    if kb:
        return kb / (1024 * 1024)
    # macOS. Deliberately not a subprocess on Linux, where /proc is definitive.
    if platform.system() == 'Darwin':
        try:
            import subprocess
            out = subprocess.run(['sysctl', '-n', 'hw.memsize'],
                                 capture_output=True, text=True, timeout=5)
            if out.returncode == 0 and out.stdout.strip().isdigit():
                return int(out.stdout.strip()) / GIB
        except Exception:
            pass
    return None


def _device(path):
    """st_dev of a path, walking up to its nearest existing ancestor.

    Neither the outdir nor the work directory need exist yet -- both are created by
    the run being sized -- so the device of the deepest parent that does exist is
    the one they will land on.
    """
    q = Path(path).expanduser().absolute()
    while not q.exists() and q != q.parent:
        q = q.parent
    try:
        return os.stat(q).st_dev
    except OSError:
        return None


def work_dir_for(outdir):
    """Where the work directory has to live for `publish_mode = 'link'` to work.

    `publish_mode = 'link'` hard-links published output into outdir, which is 68%
    of this pipeline's storage saving -- every genome would otherwise be stored
    twice, once in work/ and once published (~52 against ~16 MB/genome, measured
    on Erwiniaceae). A hard link cannot cross a filesystem, so an outdir on a
    second drive with the work directory left at Nextflow's default ./work fails
    EVERY publish with `Invalid cross-device link`, after the task has run.

    Declaring the work directory beside outdir fixes that without giving up the
    hard links, and has a second benefit: the pplacer scratch and the concurrent
    antiSMASH batches land in the work directory, so the free-space figure this
    script probes at outdir is then the one that actually governs them. With them
    on different filesystems a single probe is measuring the wrong disk.

    Returns (Path, why) or (None, None) when Nextflow's default is already on the
    right filesystem and nothing need be declared.
    """
    outdir = Path(outdir).expanduser().absolute()
    dev = _device(outdir)
    if dev is None or dev == _device(Path.cwd()):
        return None, None
    for cand in (outdir.parent / 'work', outdir / 'work'):
        if _device(cand) == dev:
            return cand, (f'outdir is on a different filesystem from the working '
                          f'directory ({Path.cwd()}), and publish_mode is '
                          f"'link' -- a hard link cannot cross devices")
    return None, None


def detect(force_cores=None, force_ram_gb=None, force_disk_gb=None, outdir='.'):
    """Read the machine. `force_*` exist so the tests can size a box this is not."""
    notes = []

    host_cores = os.cpu_count() or 1
    cg_cores = _cgroup_cores()
    cores = force_cores or (min(host_cores, cg_cores) if cg_cores else host_cores)
    if force_cores:
        notes.append(f'cores forced to {force_cores}')
    elif cg_cores and cg_cores < host_cores:
        notes.append(f'cgroup caps CPU at {cg_cores} of {host_cores} host cores')

    host_ram = _host_ram_gb()
    cg_ram = _cgroup_ram_gb()
    if force_ram_gb:
        ram_gb = float(force_ram_gb)
        notes.append(f'RAM forced to {force_ram_gb} GB')
    elif host_ram and cg_ram:
        ram_gb = min(host_ram, cg_ram)
        if cg_ram < host_ram:
            notes.append(f'cgroup caps memory at {cg_ram:.0f} GB '
                         f'of {host_ram:.0f} GB host')
    elif host_ram or cg_ram:
        ram_gb = host_ram or cg_ram
    else:
        # Guess LOW. The admission check is one-sided: under-requesting runs
        # slowly, over-requesting does not run at all.
        ram_gb = MIN_RAM_GB
        notes.append(f'could not read total memory; assuming the {MIN_RAM_GB:.0f} GB floor')

    if force_disk_gb:
        disk_free_gb = float(force_disk_gb)
        notes.append(f'free disk forced to {force_disk_gb} GB')
    else:
        probe = Path(outdir)
        while not probe.exists() and probe != probe.parent:
            probe = probe.parent
        try:
            disk_free_gb = shutil.disk_usage(probe).free / GIB
        except OSError:
            disk_free_gb = 0.0
            notes.append(f'could not read free space at {outdir}')

    return dict(cores=cores, ram_gb=ram_gb, disk_free_gb=disk_free_gb,
                disk_probe=str(outdir), platform=platform.platform(), notes=notes)


def budgets(cores, ram_gb):
    """Per-label memory budgets, in GB, plus the two derived concurrency numbers.

    Each is min(what fits, what the process was measured to need), so a larger
    machine never inflates a request past the measurement and a smaller one is
    clamped rather than refused. The OS keeps 2 GB, or 12% on a large box.
    """
    usable = ram_gb - max(2.0, ram_gb * 0.12)
    return dict(
        usable_gb=usable,
        # NOT clamped to the machine. Declaring less than GTDB-Tk measurably needs
        # is the bug this script exists to prevent in the other direction: the
        # local executor treats `memory` as an ADMISSION constraint and does not
        # enforce it, so a 16 GB declaration let a 36.8 GiB task run while telling
        # Nextflow it could admit three of them at once -- ~110 GiB on a 57 GB box,
        # which any run over gtdbtk_shard_size genomes would have triggered.
        # Declaring the true figure both reserves enough and caps concurrency at
        # one shard for free. A machine that cannot afford it gets a degradation
        # below, proposing run_gtdbtk = false, rather than a budget that lies.
        gtdbtk_gb=GTDBTK_ALLOC_GB,
        # Capped at the 16 GB the local profile already gave BiG-SCAPE, so the
        # published 8,350-BGC ceiling is unchanged on a box that can afford it.
        bigscape_gb=max(4.0, min(usable, 16.0)),
        medium_gb=max(4.0, min(usable, 8.0)),
        # antiSMASH doubles on its one retry, so it has to fit the machine TWICE.
        # Measured peak RSS over 950 genomes: p90 1.2 GB, p99 4.7 GB.
        antismash_gb=max(2.0, min(usable / 2.0, 6.0)),
        # Reproduces the measured 7 on a 16-core box: "7 x 2 = 14 cores, leaving
        # headroom for the fetch tasks running alongside".
        antismash_forks=max(1, (cores - 2) // 2),
        # Oversubscription, because FETCH_RENAME_SCREEN is latency-bound: at 12
        # concurrent fetches nothing was saturated (network ~11 of ~56 MB/s, disk
        # 57.7% busy, CPU ~7 of 16, memory 5 of 56 GB). A blocked network read
        # occupies no core. Any CPU-BOUND process added needs its own maxForks.
        executor_cpus=cores * 2,
    )


def disk_needed(answers, schema, antismash_forks=1, outdir_gb=None):
    """What the planned run needs on disk, split by whether GTDB-Tk runs.

    Database sizes come from configure_run.py's own table rather than a second
    copy, so the two cannot disagree about what a run costs.
    """
    import configure_run as cr

    fmt = {k: answers.get(k, schema.get(k, {}).get('default'))
           for k in ('taxdump_date', 'pfam_release', 'gtdb_release')}
    outdir = Path(answers.get('outdir', 'results'))

    rows = []
    for db in cr.DATABASES:
        if not db['needed'](answers):
            continue
        cached = (outdir / db['path'].format(**fmt)).exists()
        rows.append(dict(key=db['key'].format(**fmt), disk_gb=db['disk'] / GIB,
                         download_gb=db['download'] / GIB, cached=cached,
                         gtdbtk=db['module'] == 'download_gtdbtk_db'))

    def total(include_gtdbtk):
        return sum(r['disk_gb'] for r in rows
                   if not r['cached'] and (include_gtdbtk or not r['gtdbtk']))

    scratch = max(WORK_SCRATCH_FLOOR_GB, WORK_SCRATCH_PER_FORK_GB * antismash_forks)
    out_gb = OUTDIR_ORDER_SCALE_GB if outdir_gb is None else float(outdir_gb)
    gtdbtk_planned = any(r['gtdbtk'] for r in rows)
    return dict(
        rows=rows,
        gtdbtk_planned=gtdbtk_planned,
        scratch_gb=scratch,
        outdir_gb=out_gb,
        outdir_assumed=outdir_gb is None,
        # GTDB-Tk's 150 GB is set by GTDB's reference data, not by how many genomes
        # reach it, so it does NOT scale with forks.
        with_gtdbtk=total(True) + GTDBTK_WORK_SCRATCH_GB + out_gb,
        without_gtdbtk=total(False) + scratch + out_gb,
    )


def assess(machine, answers, schema, outdir_gb=None):
    """Verdict, budgets and the degradations needed to make the run fit."""
    cores, ram_gb, free_gb = machine['cores'], machine['ram_gb'], machine['disk_free_gb']
    b = budgets(cores, ram_gb)
    disk = disk_needed(answers, schema, b['antismash_forks'], outdir_gb)

    blockers, degradations = [], []

    if cores < MIN_CORES:
        blockers.append(
            f'{cores} cores, floor is {MIN_CORES}. The widest process asks for 8 cpus '
            f'and the profile oversubscribes CPU 2x, so {MIN_CORES} cores offer 8 slots '
            f'and admit it exactly; {cores} offer {cores * 2} and the task is refused, '
            f'which kills the run rather than slowing it.')
    if ram_gb < MIN_RAM_GB:
        # Deliberately a declared floor rather than a computed one. The budgets do
        # clamp below this -- an 8 GB box gets antiSMASH 3 GB, 6 GB on retry, which
        # arithmetically fits -- but nothing has been VALIDATED there, and the
        # margins that absorb a missed memory spike are gone. Refusing is honest;
        # running and failing eight hours in is not.
        blockers.append(
            f'{ram_gb:.1f} GB RAM, below the supported floor of {MIN_RAM_GB:.0f} GB. '
            f'The budgets would clamp to fit (antiSMASH {b["antismash_gb"]:.0f} GB, '
            f'{b["antismash_gb"] * 2:.0f} GB on retry) but that configuration is '
            f'untested, and the 60 s memory sampling behind these figures cannot '
            f'rule out a spike the remaining headroom would not absorb.')

    # RAM is now a gate on GTDB-Tk too, not just disk. Checked against USABLE rather
    # than total, because the OS's share is not available to a task.
    if disk['gtdbtk_planned'] and b['gtdbtk_gb'] > b['usable_gb']:
        degradations.append(dict(
            param='run_gtdbtk', value=False,
            why=(f'GTDB-Tk needs {GTDBTK_ALLOC_GB:.0f} GB of RAM (measured peak '
                 f'{GTDBTK_PEAK_GB:.1f} GiB in a ~30 s spike) and only '
                 f'{b["usable_gb"]:.1f} GB is budgetable on this machine, so the '
                 f'task would not be admitted at all'),
            costs=('no GTDB-Tk taxonomy per GCF, so the HGT-candidate analysis '
                   'and the taxonomy columns of the report are unavailable; '
                   'clustering, BGC detection and the branch-point calls are '
                   'unaffected')))

    # Disk is the other gate on GTDB-Tk, and the only degradation that changes the
    # SCIENCE rather than the wall clock. Say so plainly.
    if disk['gtdbtk_planned'] and free_gb < disk['with_gtdbtk']:
        if free_gb >= disk['without_gtdbtk']:
            degradations.append(dict(
                param='run_gtdbtk', value=False,
                why=(f'GTDB-Tk needs about {disk["with_gtdbtk"]:.0f} GB free '
                     f'({GTDBTK_WORK_SCRATCH_GB:.0f} GB of it pplacer work scratch, '
                     f'measured) and {free_gb:.0f} GB is available'),
                costs=('no GTDB-Tk taxonomy per GCF, so the HGT-candidate analysis '
                       'and the taxonomy columns of the report are unavailable; '
                       'clustering, BGC detection and the branch-point calls are '
                       'unaffected')))
        else:
            blockers.append(
                f'{free_gb:.0f} GB free at {machine["disk_probe"]}; even with GTDB-Tk '
                f'off the run needs about {disk["without_gtdbtk"]:.0f} GB')

    # BiG-SCAPE: partitioning is the lever when the cap is the problem, but it is
    # validated only to 1,302 BGCs and has never run on SLURM, so it is proposed
    # rather than assumed -- and only when a PLAUSIBLE yield would exceed the cap.
    #
    # An earlier version compared the cap to bigscape_partition_threshold, which is
    # not a comparison that means anything: the cap is under 10,000 on every machine
    # the local profile supports, so the advisory fired every time. The measured
    # yields are the right yardstick -- Enterobacterales GenBank produced 1,476 BGCs
    # from 1.35 M genomes, and Actinomycetes projects to 1,000-4,000 -- so flag only
    # a cap that those runs could actually hit.
    cap = bigscape_cap(b['bigscape_gb'])
    if (answers.get('clustering', 'bigscape') == 'bigscape'
            and not answers.get('bigscape_partition', False)
            and cap < PLAUSIBLE_BGC_YIELD):
        degradations.append(dict(
            param='bigscape_partition', value=True,
            why=(f'BiG-SCAPE gets {b["bigscape_gb"]:.0f} GB here, a partition cap of '
                 f'~{cap:,} BGCs (GB = {BIGSCAPE_INTERCEPT_GB} + {BIGSCAPE_COEFF:g}*n^2 '
                 f'on the PSS basis, at the 0.85 margin). The measured runs reach '
                 f'{PLAUSIBLE_BGC_YIELD:,} BGCs, so this run could exceed it'),
            costs=('partitioning is validated to 1,302 BGCs and has never run on '
                   'SLURM; a force-chunked component invents families, so treat '
                   "PARTITION_BGCS's \"exceeded --max_partition_size\" warning as "
                   'invalidating the run'),
            advisory=True))

    work_dir, work_why = work_dir_for(answers.get('outdir', 'results'))

    verdict = 'refused' if blockers else ('degraded' if degradations else 'ok')
    return dict(machine=machine, budgets=b, disk=disk, bigscape_cap=cap,
                work_dir=str(work_dir) if work_dir else None, work_why=work_why,
                bigscape_capacity=bigscape_capacity(b['bigscape_gb']),
                blockers=blockers, degradations=degradations, verdict=verdict)


def render_config(result):
    """A machine.config overriding the `local` profile."""
    b, m = result['budgets'], result['machine']
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    # Declared only when Nextflow's default would be on the wrong filesystem, so a
    # single-drive machine's config stays free of a setting it does not need.
    work = ''
    if result.get('work_dir'):
        why = textwrap.fill(
            f"{result['work_why']}, so the work directory is declared beside outdir. "
            f"Without this every publish fails with `Invalid cross-device link` AFTER "
            f"the task has run. -w on the command line still overrides it.",
            width=76, initial_indent='// ', subsequent_indent='// ')
        work = f"//\n{why}\nworkDir = '{result['work_dir']}'\n\n"
    return f"""// Generated by scripts/check_machine.py on {stamp}. Do not hand-edit:
// re-run the script instead, so the numbers keep their provenance.
//
// Detected: {m['cores']} cores, {m['ram_gb']:.1f} GB RAM, \
{m['disk_free_gb']:.0f} GB free at {m['disk_probe']}
// Budgetable after the OS's share: {b['usable_gb']:.1f} GB
// BiG-SCAPE at {b['bigscape_gb']:.0f} GB: holds ~{result['bigscape_capacity']:,} BGCs,
// partition cap ~{result['bigscape_cap']:,} (capacity x the 0.85 margin)
//
// This overrides the `local` profile, which carries the PORTABLE FLOOR (sized for
// 16 GB / 4 cores) rather than this machine. A -c file wins over the profile
// regardless of profile order, so `-profile local -c machine.config` is enough.
{work}profiles {{
    local {{
        executor.cpus = {b['executor_cpus']}
        // Declared, not left to detection. Nextflow reads the HOST's memory, so
        // inside a container or cgroup it would over-admit and the tasks would be
        // OOM-killed instead of queued. This is the figure `detect()` resolved,
        // which already takes the cgroup limit into account.
        executor.memory = '{m['ram_gb']:.0f} GB'

        process {{
            withLabel: 'process_high_memory' {{ memory = '{b['gtdbtk_gb']:.0f} GB' }}
            withLabel: 'process_high'        {{ memory = '{b['bigscape_gb']:.0f} GB' }}
            withLabel: 'process_medium'      {{ memory = '{b['medium_gb']:.0f} GB' }}

            withName: 'ANTISMASH' {{
                memory   = {{ {b['antismash_gb']:.0f}.GB * task.attempt }}
                maxForks = {b['antismash_forks']}
            }}
        }}
    }}
}}
"""


def report(result):
    """Human-readable verdict. Returns a list of lines."""
    m, b, disk = result['machine'], result['budgets'], result['disk']
    out = ['', '=' * 72, '  Machine check', '=' * 72,
           f'  {m["cores"]} cores | {m["ram_gb"]:.1f} GB RAM | '
           f'{m["disk_free_gb"]:.0f} GB free at {m["disk_probe"]}']
    for note in m['notes']:
        out.append(f'    note: {note}')
    out += ['',
            f'  Budgets ({b["usable_gb"]:.1f} GB budgetable after the OS share):',
            f'    GTDB-Tk        {b["gtdbtk_gb"]:5.0f} GB   '
            f'(measured peak {GTDBTK_PEAK_GB:.1f} GiB, ~30 s spike)',
            f'    BiG-SCAPE      {b["bigscape_gb"]:5.0f} GB   '
            f'(holds ~{result["bigscape_capacity"]:,} BGCs; '
            f'partition cap ~{result["bigscape_cap"]:,})',
            f'    antiSMASH      {b["antismash_gb"]:5.0f} GB   '
            f'x2 on retry, {b["antismash_forks"]} concurrent',
            f'    process_medium {b["medium_gb"]:5.0f} GB',
            f'    executor.cpus  {b["executor_cpus"]:5d}      (2x oversubscribed; '
            f'fetch is latency-bound)']

    if result.get('work_dir'):
        out += ['', '  Work directory:', f'    {result["work_dir"]}',
                f'    why:  {result["work_why"]}',
                '    Without it every publish fails with `Invalid cross-device link`,',
                '    after the task has already run. -w overrides it.']

    if disk['rows']:
        out += ['', '  Databases this run needs:']
        for r in disk['rows']:
            state = 'cached' if r['cached'] else f'{r["download_gb"]:.1f} GB to fetch'
            out.append(f'    {r["key"]:<22} {r["disk_gb"]:7.1f} GB on disk   ({state})')
        need = disk['with_gtdbtk'] if disk['gtdbtk_planned'] else disk['without_gtdbtk']
        scratch = (GTDBTK_WORK_SCRATCH_GB if disk['gtdbtk_planned']
                   else disk['scratch_gb'])
        out.append(f'    {"work scratch":<22} {scratch:7.1f} GB'
                   + ('   (GTDB-Tk pplacer; set by GTDB reference data)'
                      if disk['gtdbtk_planned']
                      else f'   ({b["antismash_forks"]} concurrent batches)'))
        out.append(f'    {"output":<22} {disk["outdir_gb"]:7.1f} GB'
                   + ('   (ASSUMED order-scale; --outdir-gb to override)'
                      if disk['outdir_assumed'] else ''))
        out.append(f'    {"total needed":<22} {need:7.0f} GB')

    if result['blockers']:
        out += ['', '  REFUSED:']
        out += [f'    - {b_}' for b_ in result['blockers']]
    for d in result['degradations']:
        label = 'SUGGESTED' if d.get('advisory') else 'REQUIRED TO FIT'
        out += ['', f'  {label}: {d["param"]} = {json.dumps(d["value"])}',
                f'    why:  {d["why"]}', f'    cost: {d["costs"]}']
    out += ['', '=' * 72, '']
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('-p', '--params-file', help='the run that is about to be launched')
    ap.add_argument('-o', '--out', default='machine.config', help='config to write')
    ap.add_argument('--profile', default='local',
                    help="skip entirely for anything but 'local'")
    ap.add_argument('--json', action='store_true', help='machine-readable verdict')
    ap.add_argument('--print-only', action='store_true', help='do not write the config')
    ap.add_argument('--yes', action='store_true', help='accept degradations unprompted')
    ap.add_argument('--force-cores', type=int)
    ap.add_argument('--force-ram-gb', type=float)
    ap.add_argument('--force-disk-gb', type=float)
    ap.add_argument('--outdir-gb', type=float,
                    help='expected output size, overriding the order-scale assumption')
    args = ap.parse_args()

    if args.profile != 'local':
        print(f'profile {args.profile!r}: sizing skipped. A cluster profile declares '
              f'capacity, not a property of the submit host.')
        return 0

    import configure_run as cr
    schema = cr.load_schema()
    answers = {}
    if args.params_file:
        answers = json.loads(Path(args.params_file).read_text())
    for k, spec in schema.items():
        answers.setdefault(k, spec.get('default'))

    machine = detect(args.force_cores, args.force_ram_gb, args.force_disk_gb,
                     answers.get('outdir', 'results'))
    result = assess(machine, answers, schema, args.outdir_gb)

    if args.json:
        print(json.dumps({k: v for k, v in result.items() if k != 'disk'}
                         | {'disk_needed_gb': result['disk']['with_gtdbtk']
                            if result['disk']['gtdbtk_planned']
                            else result['disk']['without_gtdbtk']}, indent=2, default=str))
        return 2 if result['verdict'] == 'refused' else 0

    print('\n'.join(report(result)))

    if result['verdict'] == 'refused':
        print('Nothing written. The run would die on admission, not run slowly.')
        return 2

    if result['degradations'] and not args.yes:
        required = [d for d in result['degradations'] if not d.get('advisory')]
        if required and sys.stdin.isatty():
            reply = input('Apply the changes above and continue? [y/N] ').strip().lower()
            if reply not in ('y', 'yes'):
                print('Nothing written.')
                return 1

    if not args.print_only:
        Path(args.out).write_text(render_config(result))
        print(f'wrote {Path(args.out).resolve()}')
        applied = [f'--{d["param"]} {json.dumps(d["value"])}'
                   for d in result['degradations'] if not d.get('advisory')]
        if applied:
            print('Add to the run (or set them in the params file): '
                  + ' '.join(applied))
    else:
        print(render_config(result))
    return 0


if __name__ == '__main__':
    sys.exit(main())
