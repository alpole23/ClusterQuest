#!/usr/bin/env python3
"""Verify Utils.BOOLEAN_PARAMS still lists every boolean param in nextflow.config.

main.nf rejects a boolean param whose value is not a real boolean, because a
command-line `--run_gtdbtk false` arrives as the STRING "false" and every non-empty
string is true in Groovy — so the flag ENABLES what it appears to disable, silently.
That protection only covers params named in the list, so a boolean added to the config
and not to the list is unprotected and looks protected, which is the worst state.

This is a list-drift check, not a behaviour test: the behaviour is checked by
tests/test_boolean_params.nf, which runs the pipeline with `--run_gtdbtk false` and
asserts it aborts.

Exit 0 if the two agree, 1 otherwise.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def config_booleans():
    """Param names whose default in nextflow.config's top-level params block is true/false."""
    text = (ROOT / 'nextflow.config').read_text()
    start = text.index('params {')
    depth, i = 0, start
    while i < len(text):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                break
        i += 1
    block = text[start:i]
    return {m.group(1) for m in
            re.finditer(r'^\s*([a-z_][a-zA-Z0-9_]*)\s*=\s*(?:true|false)\s*(?://.*)?$',
                        block, re.M)}


def declared_booleans():
    """Param names listed in Utils.BOOLEAN_PARAMS."""
    text = (ROOT / 'lib' / 'Utils.groovy').read_text()
    m = re.search(r'BOOLEAN_PARAMS\s*=\s*\[(.*?)\]', text, re.S)
    if not m:
        print('Utils.BOOLEAN_PARAMS not found in lib/Utils.groovy')
        sys.exit(1)
    return set(re.findall(r"'([^']+)'", m.group(1)))


def main():
    in_config, declared = config_booleans(), declared_booleans()
    missing = sorted(in_config - declared)
    extra = sorted(declared - in_config)

    if missing:
        print('boolean params in nextflow.config but NOT in Utils.BOOLEAN_PARAMS:')
        for name in missing:
            print(f'  - {name}   (`--{name} false` would silently enable it)')
    if extra:
        print('names in Utils.BOOLEAN_PARAMS with no boolean default in nextflow.config:')
        for name in extra:
            print(f'  - {name}')
    if missing or extra:
        print('\nUpdate Utils.BOOLEAN_PARAMS to match nextflow.config.')
        return 1

    print(f'boolean param list current ({len(declared)} params)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
