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
from utils.report_lint import undefined_handlers, check_report  # noqa: E402

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


def main():
    failures = []
    for name, html, expected in CASES:
        got = undefined_handlers(html)
        if got != expected:
            failures.append(f"{name}: expected {expected}, got {got}")

    # check_report must turn a finding into a readable message
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
