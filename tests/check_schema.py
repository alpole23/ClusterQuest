#!/usr/bin/env python3
"""Keep nextflow_schema.json in step with the params nextflow.config declares.

The schema is the pipeline's parameter contract: main.nf calls nf-schema's
validateParameters(), and the schema's root `additionalProperties: false` is what makes
`--taxn Pantoea` an error instead of a silent run on the default taxon.

That cuts both ways. A param added to nextflow.config and not to the schema is REJECTED
at runtime, so the pipeline breaks for everyone the moment someone adds a config option
and forgets the schema. This check compares the two and fails on either kind of drift.

    python tests/check_schema.py            # compare, exit 1 on drift
    python tests/check_schema.py --write    # regenerate the schema from the config

Regeneration reads the authoritative param list from `nextflow config -flat` rather than
by parsing the config text: a regex over the config missed one param, and because the
same regex checked its own output, the gap was invisible until a run failed.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = ROOT / 'nextflow_schema.json'


def resolved_params():
    """{name: value-as-written} straight from Nextflow's own config resolution."""
    out = subprocess.run(['nextflow', 'config', '-flat'], cwd=ROOT,
                         capture_output=True, text=True)
    if out.returncode != 0:
        print('nextflow config failed:\n' + out.stderr[:2000])
        sys.exit(1)
    params = {}
    for line in out.stdout.splitlines():
        m = re.match(r'params\.([A-Za-z0-9_]+)\s*=\s*(.*)', line.strip())
        if m:
            params[m.group(1)] = m.group(2).strip()
    return params


def infer(value):
    """JSON-schema type and default for a config value as `nextflow config` prints it."""
    if value in ('true', 'false'):
        return 'boolean', value == 'true'
    if value == 'null':
        return ['string', 'null'], None
    if re.fullmatch(r'-?\d+', value):
        return 'integer', int(value)
    if re.fullmatch(r'-?\d+\.\d+', value):
        return 'number', float(value)
    return 'string', value.strip('\'"')


def descriptions():
    """Param -> the comment beside or above it in nextflow.config, for --help text."""
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
    block, docs, pending, group = text[start:i], {}, [], 'General'
    groups = {}
    for line in block.split('\n'):
        sec = re.match(r'\s*//\s*─+\s*(.+?)\s*─+\s*$', line)
        if sec:
            group, pending = sec.group(1).strip(), []
            continue
        if re.match(r'\s*//', line):
            pending.append(re.sub(r'^\s*//\s?', '', line).strip())
            continue
        m = re.match(r'\s*([a-z_][a-zA-Z0-9_]*)\s*=\s*(.+?)\s*$', line)
        if m:
            inline = m.group(2).split('//', 1)[1].strip() if '//' in m.group(2) else ''
            docs[m.group(1)] = (inline or (pending[0] if pending else ''))[:200]
            groups[m.group(1)] = group
        if not line.strip():
            pending = []
    return docs, groups


def build():
    params = resolved_params()
    docs, groups = descriptions()
    props = {}
    for name in sorted(params):
        typ, default = infer(params[name])
        spec = {'type': typ, 'description': docs.get(name, '')}
        if default is not None:
            spec['default'] = default
        if groups.get(name):
            spec['x-group'] = groups[name]
        props[name] = spec
    return {
        '$schema': 'https://json-schema.org/draft/2020-12/schema',
        '$id': 'https://raw.githubusercontent.com/alpole23/ClusterQuest/main/nextflow_schema.json',
        'title': 'ClusterQuest pipeline parameters',
        'description': ('Phosphonate BGC discovery and characterisation. Generated from '
                        'nextflow.config by tests/check_schema.py --write; regenerate '
                        'rather than hand-edit.'),
        'type': 'object',
        'properties': props,
        # Root level, with every property listed here rather than inside $defs groups:
        # additionalProperties only sees properties declared in the same schema object,
        # so an allOf/$defs layout rejects every grouped param instead of unknown ones.
        'additionalProperties': False,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--write', action='store_true', help='regenerate the schema')
    args = ap.parse_args()

    built = build()
    if args.write:
        SCHEMA.write_text(json.dumps(built, indent=4) + '\n')
        print(f'wrote {SCHEMA.name}: {len(built["properties"])} params')
        return 0

    if not SCHEMA.is_file():
        print(f'{SCHEMA.name} is missing; run tests/check_schema.py --write')
        return 1
    current = json.loads(SCHEMA.read_text())
    have, want = set(current.get('properties', {})), set(built['properties'])
    missing, extra = sorted(want - have), sorted(have - want)
    if missing:
        print('params in nextflow.config but NOT in the schema '
              '(the pipeline would reject them at runtime):')
        for n in missing:
            print(f'  - {n}')
    if extra:
        print('params in the schema that nextflow.config no longer declares:')
        for n in extra:
            print(f'  - {n}')

    wrong = []
    for name in sorted(have & want):
        a, b = current['properties'][name].get('type'), built['properties'][name].get('type')
        if a != b:
            wrong.append(f'  - {name}: schema says {a}, config default implies {b}')
    if wrong:
        print('type mismatches between schema and config:')
        print('\n'.join(wrong))

    if missing or extra or wrong:
        print('\nRegenerate with: python tests/check_schema.py --write')
        return 1
    if not current.get('additionalProperties') is False:
        print('schema root is missing "additionalProperties": false, so an unknown '
              'param would be accepted silently')
        return 1
    print(f'parameter schema current ({len(want)} params)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
