#!/usr/bin/env python3
"""Check docs/parameters.md against nextflow.config.

WHY THIS EXISTS
---------------
The measurements behind each default used to live directly above the value they
explained, which made them impossible to get out of sync. They were moved to
docs/parameters.md on 2026-11-02 for a good reason -- as block comments they became
the wizard's prompts, truncated mid-word at 240 characters, and 19 of the parameters
the wizard asks about were in that state.

The cost of the move is exactly the drift this file exists to catch: a default can
now be changed in the config while the note in the docs goes on describing the old
one. That is worse than a long prompt, because a stale note reads as authoritative.
So the correspondence is checked rather than trusted:

  * every parameter in nextflow.config has a section in the docs;
  * the default the docs print matches the config's value, character for character;
  * no section describes a parameter that no longer exists;
  * no parameter has regained a block comment above it, which would silently become
    its prompt again;
  * no description is long enough to truncate.

Run by tests/run_tests.sh.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))

# The one parameter declared outside the params block, as a derived path
# (`params.pipeline_info_dir = "${params.outdir}/pipeline_info/..."`). It is not a
# knob -- nothing asks for it and setting it by hand would desynchronise the trace,
# report and timeline paths -- so it is exempt rather than documented as one.
DERIVED = {'pipeline_info_dir'}

# Long enough to say what the decision is, short enough that nothing truncates.
# tests/check_schema.py caps descriptions at 240 characters; staying well under that
# is what keeps the prompt a prompt.
MAX_DESCRIPTION = 140

FAIL = []


def check(ok, label, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {label}' + (f' ({detail})' if detail else ''))
    if not ok:
        FAIL.append(label)


def config_params():
    """param -> (value, inline comment, block comment) from the params block."""
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

    found, pending = {}, []
    for line in text[start:i].split('\n'):
        if re.match(r'\s*//\s*─+\s*(.+?)\s*─+\s*$', line):
            pending = []
            continue
        if re.match(r'\s*//', line):
            pending.append(re.sub(r'^\s*//\s?', '', line).rstrip())
            continue
        m = re.match(r'\s*([a-z_][a-zA-Z0-9_]*)\s*=\s*(.+?)\s*$', line)
        if m:
            rhs = m.group(2)
            inline = rhs.split('//', 1)[1].strip() if '//' in rhs else ''
            value = rhs.split('//', 1)[0].strip()
            found[m.group(1)] = (value, inline, '\n'.join(pending).strip())
            pending = []
        elif not line.strip():
            pending = []
    return found


def main():
    cfg = config_params()
    doc_text = (ROOT / 'docs' / 'parameters.md').read_text()

    # Sections look like:  ### `param`   then  **Default** `value`
    doc_defaults = dict(re.findall(
        r'### `([a-z_][a-zA-Z0-9_]*)`\n\n(?:.*?\n\n)?\*\*Default\*\* `([^`]*)`',
        doc_text))
    doc_sections = set(re.findall(r'### `([a-z_][a-zA-Z0-9_]*)`', doc_text))

    # 1. Nothing undocumented.
    missing = sorted(set(cfg) - doc_sections)
    check(not missing, 'every config parameter has a docs section',
          ', '.join(missing[:6]))

    # 2. Nothing documented that no longer exists -- a renamed parameter leaves a
    #    section behind that still reads as current.
    orphans = sorted(doc_sections - set(cfg) - DERIVED)
    check(not orphans, 'no docs section for a parameter that was removed',
          ', '.join(orphans[:6]))

    # 3. The defaults agree. This is the whole point of the file.
    bad = []
    for param, (value, _inline, _note) in sorted(cfg.items()):
        shown = doc_defaults.get(param)
        if shown is None:
            bad.append(f'{param}: no **Default** line')
        elif shown != (value or 'null'):
            bad.append(f'{param}: docs {shown!r} vs config {value!r}')
    check(not bad, 'documented defaults match nextflow.config', '; '.join(bad[:4]))

    # 4. No parameter has regained a block comment. tests/check_schema.py prefers an
    #    inline comment, so a block one only matters when there is no inline -- but
    #    then it becomes the prompt again, which is what the move undid.
    regressed = sorted(p for p, (_v, inline, note) in cfg.items()
                       if note and not inline)
    check(not regressed,
          'no parameter takes its prompt from a block comment',
          ', '.join(regressed[:6]) + ' (put the note in docs/parameters.md '
          'and leave a one-line inline comment)' if regressed else '')

    # 5. Every parameter actually has a prompt, and none is long enough to truncate.
    schema = json.loads((ROOT / 'nextflow_schema.json').read_text())['properties']
    blank = sorted(p for p in cfg if not schema.get(p, {}).get('description'))
    check(not blank, 'every parameter has a non-empty description',
          ', '.join(blank[:6]))
    toolong = sorted((p, len(schema[p]['description'])) for p in cfg
                     if len(schema.get(p, {}).get('description', '')) > MAX_DESCRIPTION)
    check(not toolong, f'no description exceeds {MAX_DESCRIPTION} characters',
          ', '.join(f'{p} ({n})' for p, n in toolong[:4]))

    # 6. The notes that were moved are still there. 17,647 bytes were relocated on
    #    2026-11-02; a docs file that shrinks far below that has lost them, which a
    #    per-parameter check would not notice if whole sections went at once.
    measured = len(doc_text)
    check(measured > 25_000, 'the relocated notes are still in the docs file',
          f'{measured:,} B, expected >25,000')

    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
