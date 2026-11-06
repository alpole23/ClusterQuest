#!/usr/bin/env python3
"""Self-test for utils/report_lint, plus an optional check of a real report.

The linter runs at report-generation time, so its own correctness is what needs
covering here: a checker that silently stops detecting anything is worse than none.

Usage:
    python tests/check_report_js.py [report.html ...]

Exits 0 when the linter behaves and any supplied reports are clean.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'scripts'))
from utils.report_lint import (  # noqa: E402
    undefined_handlers, check_report, unbalanced_braces)

CASES = [
    # (name, html, expected undefined handlers)
    ("missing definition — the filterGenomes bug",
     '<input onkeyup="filterGenomes()"><script>function other() {}</script>',
     ['filterGenomes']),
    ("definition present",
     '<input onkeyup="filterGenomes()"><script>function filterGenomes() {}</script>',
     []),
    ("several attributes, one missing",
     '<a onclick="toggleNode(\'n\')"></a><select onchange="sortGCFs(this.value)"></select>'
     '<script>function toggleNode(x) {}</script>',
     ['sortGCFs']),
    # The 2026-11-04 failure: the handler IS defined, so the undefined-handler check
    # passes -- but a stray brace makes the whole block a SyntaxError and nothing runs.
    ("orphaned function body — defined handler, dead script",
     '<input onkeyup="filterGenomes()"><script>function filterGenomes() {}\n'
     'const x = 1;\n}</script>',
     []),
    ("method calls are not our business",
     '<button onclick="this.form.reset()"></button>',
     []),
    ("browser globals are callable without definition",
     '<a onclick="alert(1)"></a>',
     []),
    ("no handlers at all",
     '<p>nothing here</p>',
     []),
    ("whitespace inside the attribute",
     '<a onclick=" toggleGCF( \'g\' ) "></a><script>function toggleGCF(x) {}</script>',
     []),
]


# A brace that closes nothing is a SyntaxError, and one anywhere in a <script> kills
# every script in the document. Braces are legal inside strings, template literals,
# both comment forms and regex literals, so each of those must NOT trip it -- a checker
# that cries wolf on `const r = /[}]/g` would be turned off within a week.
BRACE_CASES = [
    ("stray brace at top level",        "function a(){}\n}\n",            2),
    ("orphaned body, the real bug",     "function a(){}\n\nconst i = 1;\nfor (;;) {}\n}", 5),
    ("balanced",                        "function a(){ if (x) { y() } }",  None),
    ("brace in a string",               'const s = "}";',                  None),
    ("brace in a single-quoted string", "const s = '}';",                  None),
    ("brace in a template literal",     "const t = `}`;",                  None),
    ("brace in a line comment",         "// }\nfunction a(){}",            None),
    ("brace in a block comment",        "/* } */ function a(){}",          None),
    ("brace in a regex",                "const r = /[}]/g;",               None),
    ("division is not a regex",         "const q = a / b; const w = c / d;", None),
]


def _REPORT_JS():
    """The shipped REPORT_JS string, loaded without importing the viz package."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'report_assets',
        str(Path(__file__).resolve().parent.parent / 'scripts' / 'viz' / 'report_assets.py'))
    ra = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ra)
    return ra.REPORT_JS


# ── The JSON island is not all strings, and the JS reading it must know ──────
#
# `viz/tables.py` ships the Genomes table as array-of-arrays from json.dumps, and one
# column is `int(total_bgcs)` -- a BARE JSON NUMBER. `searchNorm` called
# `text.toLowerCase()` on whichever cell it was handed, so every Genomes query threw
#     TypeError: text.toLowerCase is not a function
# at the first row whose four string columns all missed the query, `filter` abandoned
# the whole scan, and nothing was rendered. Typing killed the box; clearing it brought
# the table back, because the empty query never reaches a numeric cell.
#
# Nothing in this suite could see it. The handler IS defined, the braces ARE balanced,
# and a TypeError thrown inside a 150 ms setTimeout is reported to the console and
# discarded -- the page simply stops changing. Audit 2, A2-7, 2026-11-05.
#
# So this pins the CONTRACT rather than the spelling of the fix: that the producer
# really does emit a non-string column, and that every JS helper fed an island cell
# coerces before treating it as a string. `esc` is included because it was already
# written this way; it is the pattern `searchNorm` should have followed.
ISLAND_CONSUMERS = ('searchNorm', 'esc')

# Called on a String. Any of these applied straight to the parameter is the bug.
STRING_METHODS = ('toLowerCase', 'toUpperCase', 'trim', 'replace', 'indexOf',
                  'split', 'charAt', 'slice', 'startsWith', 'endsWith',
                  'substring', 'includes', 'match', 'padStart')


def _js_function(js, name):
    """The source of `function <name>(param) { ... }`, by brace matching."""
    i = js.find(f'function {name}(')
    if i < 0:
        return None, None
    param = js[js.find('(', i) + 1:js.find(')', i)].strip()
    j = js.find('{', i)
    depth = 0
    for k in range(j, len(js)):
        if js[k] == '{':
            depth += 1
        elif js[k] == '}':
            depth -= 1
            if depth == 0:
                return js[j:k + 1], param
    return None, param


def island_coercion(js):
    """Problems with how the island consumers handle a non-string cell."""
    out = []
    for name in ISLAND_CONSUMERS:
        body, param = _js_function(js, name)
        if body is None:
            out.append(f'{name}() is gone from REPORT_JS; if the island consumers '
                       f'changed, update ISLAND_CONSUMERS here')
            continue
        if 'String(' not in body:
            out.append(f'{name}({param}) does not coerce with String(). The Genomes '
                       f'JSON island carries a bare JSON number (int(total_bgcs) in '
                       f'viz/tables.py), and a string method on it raises TypeError '
                       f'silently inside the search debounce')
        for meth in STRING_METHODS:
            if f'{param}.{meth}(' in body:
                out.append(f'{name}() calls {param}.{meth}() directly on its '
                           f'parameter; wrap it in String(...) first')
    return out


def island_producer_has_a_number():
    """The reason the coercion is required, asserted at the producer.

    The cast is NOT in the island row itself -- that row is
    `[[r['genome_name'], ..., r['total_bgcs'], ...]]`, all dict lookups. The number
    is made one step earlier, where the dict is built (`'total_bgcs': int(...)`), so
    a check that reads only the row expression finds nothing and passes while the bug
    is live. Follow the field names.
    """
    tables = (Path(__file__).resolve().parent.parent
              / 'scripts' / 'viz' / 'tables.py').read_text(encoding='utf-8')
    i = tables.find('data = [[')
    if i < 0:
        return ["viz/tables.py no longer builds the Genomes island with "
                "`data = [[...]]`; re-check which of its columns are non-strings"]
    row = tables[i:tables.find(']]', i)]
    fields = re.findall(r"r\['([^']+)'\]", row)
    if not fields:
        return [f"could not read the island's column fields out of {row!r}"]
    numeric = [f for f in fields
               if re.search(rf"'{re.escape(f)}'\s*:\s*(?:int|float)\s*\(", tables)]
    if not numeric:
        return [f"none of the Genomes island columns {fields} is assigned a numeric "
                f"cast in viz/tables.py. If every column is now a string the String() "
                f"coercion above is merely harmless -- but confirm it before relaxing "
                f"this check, because json.dumps emits whatever type it is handed"]
    return []


def island_in_report_has_a_number(html):
    """Same contract, read off a real rendered report rather than the source."""
    m = re.search(r'<script id="genomeData" type="application/json">(.*?)</script>',
                  html, re.S)
    if not m:
        return []                      # no Genomes island in this report; nothing to say
    try:
        rows = json.loads(m.group(1))
    except ValueError as exc:
        return [f'genomeData island is not valid JSON: {exc}']
    if not rows:
        return []
    kinds = {i: type(v).__name__ for i, v in enumerate(rows[0])}
    nonstr = {i: k for i, k in kinds.items() if k != 'str'}
    if not nonstr:
        return []                      # all strings here; the coercion is still correct
    body, param = _js_function(_REPORT_JS(), 'searchNorm')
    if body and 'String(' not in body:
        return [f'genomeData carries non-string columns {nonstr} and searchNorm does '
                f'not coerce -- every Genomes query will throw']
    return []


def main():
    failures = []
    failures += island_coercion(_REPORT_JS())
    failures += island_producer_has_a_number()
    for name, html, expected in CASES:
        got = undefined_handlers(html)
        if got != expected:
            failures.append(f"{name}: expected {expected}, got {got}")

    for name, js, expected in BRACE_CASES:
        got = unbalanced_braces(js)
        if got != expected:
            failures.append(f"braces/{name}: expected {expected}, got {got}")

    # check_report must name the offending script block, not just say "broken"
    msgs = check_report('<script>function a(){}\n}</script>')
    if not msgs or 'closes nothing' not in msgs[0]:
        failures.append(f"check_report did not report the stray brace: {msgs}")

    # a JSON island is data; its braces are not code
    island = '<script type="application/json" id="d">{"a": 1}</script>'
    if check_report(island):
        failures.append("check_report flagged a JSON data island")

    # and the shipped constant itself must parse
    bad = unbalanced_braces(_REPORT_JS())
    if bad is not None:
        failures.append(f"REPORT_JS has a brace that closes nothing at line {bad}")

    # check_report must also still describe an undefined handler
    msgs = check_report('<input onkeyup="nope()">')
    if not msgs or 'nope()' not in msgs[0]:
        failures.append(f"check_report did not describe the problem: {msgs}")

    for arg in sys.argv[1:]:
        path = Path(arg)
        if not path.exists():
            failures.append(f"{path}: no such report")
            continue
        html = path.read_text(encoding='utf-8')
        for problem in check_report(html):
            failures.append(f"{path.name}: {problem}")
        for problem in island_in_report_has_a_number(html):
            failures.append(f"{path.name}: {problem}")

    if failures:
        print("report JS checks failed:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"report JS checks OK ({len(CASES)} linter cases, "
          f"{len(ISLAND_CONSUMERS)} island consumers"
          f"{f', {len(sys.argv) - 1} report(s)' if len(sys.argv) > 1 else ''})")
    return 0


if __name__ == '__main__':
    sys.exit(main())
