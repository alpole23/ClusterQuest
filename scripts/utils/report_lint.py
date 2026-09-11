"""Static checks on the assembled HTML report.

`tests/check_undefined.py` catches a Python call to a name nothing defines, but the
report's JavaScript is Python *string data* — REPORT_JS in viz/report_assets.py plus
inline fragments in viz/clustering.py and visualize_results.py — so `ast` sees it as
opaque text. That blind spot let a Genomes-tab search box ship wired to
`filterGenomes()`, a function defined nowhere: every keystroke threw a ReferenceError
and filtered nothing, silently, for the life of the feature.

These checks run against the finished HTML, which is the only place the JS fragments
are all present at once.
"""

import re

# `onclick="foo(...)"` — a bare identifier immediately followed by '('. Deliberately
# does not match `this.x()` or `window.x()`, whose resolution we cannot verify here.
_HANDLER = re.compile(
    r'\bon(?:click|keyup|keydown|keypress|change|input|submit|focus|blur)\s*=\s*'
    r'"\s*([A-Za-z_$][\w$]*)\s*\(')
_FUNC_DEF = re.compile(r'\bfunction\s+([A-Za-z_$][\w$]*)\s*\(')

# Callable without a `function NAME(` definition in the document.
_BROWSER_GLOBALS = {'alert', 'confirm', 'prompt', 'print', 'open', 'close', 'requestAnimationFrame'}


def undefined_handlers(html):
    """Names invoked from inline event handlers that the document never defines.

    Returns a sorted list of names; empty when the report is sound.
    """
    defined = set(_FUNC_DEF.findall(html)) | _BROWSER_GLOBALS
    used = set(_HANDLER.findall(html))
    return sorted(used - defined)


def check_report(html):
    """Return a list of human-readable problems with the assembled report."""
    problems = []
    missing = undefined_handlers(html)
    if missing:
        problems.append(
            "inline event handlers call JavaScript functions that are never defined: "
            + ", ".join(f"{n}()" for n in missing))
    return problems
