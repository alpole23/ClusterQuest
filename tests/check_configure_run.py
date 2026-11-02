#!/usr/bin/env python3
"""Check scripts/configure_run.py against the schema and against main.nf.

The wizard reads types, defaults and help text out of nextflow_schema.json at runtime,
so those cannot drift. Two things it states for itself CAN drift, and both fail in a way
the user sees rather than the author:

  * a parameter name in QUESTIONS that the schema does not declare. The wizard would
    skip it silently -- the question simply never appears -- so the parameter becomes
    unaskable with nothing to say so.
  * a choice list that main.nf's validateParams() does not accept. The wizard would
    offer a value, the user would pick it, and the run would die on validation with the
    wizard's own suggestion in the error.

Run by tests/run_tests.sh.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))

FAIL = []


def check(ok, label, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {label}' + (f' ({detail})' if detail else ''))
    if not ok:
        FAIL.append(label)


def main():
    import configure_run as cr

    schema = json.loads((ROOT / 'nextflow_schema.json').read_text())['properties']

    # 1. Every parameter the wizard asks about exists in the schema.
    asked = [p for p, _t, _m, _pred in cr.QUESTIONS]
    missing = sorted(p for p in asked if p not in schema)
    check(not missing, 'every asked parameter is in the schema', ', '.join(missing))

    # 2. No duplicates -- a repeated entry asks the same question twice.
    dupes = sorted({p for p in asked if asked.count(p) > 1})
    check(not dupes, 'no parameter asked twice', ', '.join(dupes))

    # 3. Choice lists are accepted by main.nf. Its validateParams() builds a
    #    `validX` list per parameter; parse those and compare. A choice the pipeline
    #    rejects is worse than no choice list at all.
    main_nf = (ROOT / 'main.nf').read_text()
    for param, key in (('workflow', 'validWorkflows'),
                       ('clustering', 'validClustering'),
                       ('bigscape_alignment_mode', 'validAlignmentModes'),
                       ('bigscape_classify', 'validClassify')):
        m = re.search(rf'def {key}\s*=\s*\[(.*?)\]', main_nf, re.S)
        if not m:
            check(False, f'{key} found in main.nf')
            continue
        allowed = set(re.findall(r"'([^']*)'", m.group(1)))
        offered = set(cr.CHOICES[param])
        bad = sorted(offered - allowed)
        check(not bad, f'{param} choices accepted by main.nf',
              f'rejected: {", ".join(repr(b) for b in bad)}' if bad else '')

    # 4. Every parameter named in a WARNINGS key is one the wizard actually asks.
    warned = sorted({p for p, _v in cr.WARNINGS} - set(asked))
    check(not warned, 'warnings attach to asked parameters', ', '.join(warned))

    # 5. Defaults round-trip through coerce(): the value shown in brackets must be
    #    something the parser accepts back. This is what catches a nullable parameter
    #    rendered as "null" that coerce() would then reject.
    bad = []
    for param in asked:
        spec = schema.get(param)
        if spec is None:
            continue
        shown = cr._fmt(spec.get('default'))
        if shown == '(empty)':
            continue
        try:
            got = cr.coerce(shown, spec, param)
        except ValueError as exc:
            bad.append(f'{param}: {shown!r} -> {exc}')
            continue
        if got != spec.get('default'):
            bad.append(f'{param}: {shown!r} -> {got!r}, expected {spec.get("default")!r}')
    check(not bad, 'shown defaults parse back to themselves', '; '.join(bad[:3]))

    # 6. Every mode asks at least the parameters its own validation requires.
    #    bgc_analysis dies without input_genomes, so it must be asked there.
    bgc = {p for p, _t, m, _pred in cr.QUESTIONS if 'bgc_analysis' in m}
    check('input_genomes' in bgc, 'bgc_analysis asks for input_genomes')
    dl = {p for p, _t, m, _pred in cr.QUESTIONS if 'download' in m}
    check('taxon' in dl, 'download asks for taxon')
    check('clustering' not in dl, 'download does not ask about clustering')

    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
