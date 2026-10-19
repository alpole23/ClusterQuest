"""Run recover_orfs.py unmodified, timing every external call.

The question is how RECOVER_ORFS spends its time, and the suspicion is that
`diamond makedb` -- which the script runs once per GENOME over the whole protein
pool -- dominates. Patching subprocess.run attributes each call without touching
the script, so this measures exactly what the pipeline runs.
"""
import runpy, subprocess, sys, time
from pathlib import Path

real, log = subprocess.run, []


def timed(cmd, *a, **k):
    t0 = time.perf_counter()
    r = real(cmd, *a, **k)
    dt = time.perf_counter() - t0
    if isinstance(cmd, (list, tuple)):
        label = ' '.join(Path(str(cmd[0])).name.split()[:1] + [str(c) for c in cmd[1:2]])
    else:
        label = str(cmd)[:40]
    log.append((label, dt))
    return r


subprocess.run = timed
sys.argv = ['recover_orfs.py'] + sys.argv[1:]
t0 = time.perf_counter()
try:
    runpy.run_path('/home/alexp/pipeline/scripts/genome/recover_orfs.py', run_name='__main__')
except SystemExit:
    pass
total = time.perf_counter() - t0
acc = sum(d for _, d in log)
print('\n--- wall seconds by phase ---', file=sys.stderr)
for label, dt in log:
    print(f'{dt:9.2f}  {label}  ({100*dt/total:.1f}%)', file=sys.stderr)
print(f'{total-acc:9.2f}  python: parse/gaps/merge/write  ({100*(total-acc)/total:.1f}%)', file=sys.stderr)
print(f'{total:9.2f}  TOTAL', file=sys.stderr)
