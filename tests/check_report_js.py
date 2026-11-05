#!/usr/bin/env python3
"""Self-test for utils/report_lint, plus an optional check of a real report.

The linter runs at report-generation time, so its own correctness is what needs
covering here: a checker that silently stops detecting anything is worse than none.

Usage:
    python tests/check_report_js.py [report.html ...]

Exits 0 when the linter behaves and any supplied reports are clean.
"""
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


def main():
    failures = []
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
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'report_assets',
        str(Path(__file__).resolve().parent.parent / 'scripts' / 'viz' / 'report_assets.py'))
    ra = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ra)
    bad = unbalanced_braces(ra.REPORT_JS)
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
        for problem in check_report(path.read_text(encoding='utf-8')):
            failures.append(f"{path.name}: {problem}")

    if failures:
        print("report JS checks failed:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"report JS checks OK ({len(CASES)} linter cases"
          f"{f', {len(sys.argv) - 1} report(s)' if len(sys.argv) > 1 else ''})")
    return 0


if __name__ == '__main__':
    sys.exit(main())
