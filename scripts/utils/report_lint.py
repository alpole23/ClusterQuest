"""Static checks on the assembled HTML report.

`tests/check_undefined.py` catches a Python call to a name nothing defines, but the
report's JavaScript is Python *string data* — REPORT_JS in viz/report_assets.py plus
inline fragments in viz/clustering.py and visualize_results.py — so `ast` sees it as
opaque text. That blind spot let a Genomes-tab search box ship wired to
`filterGenomes()`, a function defined nowhere: every keystroke threw a ReferenceError
and filtered nothing, silently, for the life of the feature.

It did NOT catch the worse version of the same blind spot. On 2026-11-04 an orphaned
function body -- a `function` header deleted, its body and closing brace kept -- made
the whole <script> block a SyntaxError: `Uncaught SyntaxError: Unexpected token '}'`.
Every handler it defines was still textually present, so the undefined-handler check
passed while NO JavaScript ran at all: both search boxes, the gene-diagram hover
labels, the taxonomy tree and the JSON row renderer were dead, and the report rendered,
published and reported success. `unbalanced_braces` closes that gap.

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

# Executable script blocks only. A `type="application/json"` island is data -- the
# genome table ships as one -- and its braces are not code.
_SCRIPTS = re.compile(
    r'<script(?![^>]*type\s*=\s*"application/json")[^>]*>(.*?)</script>', re.S)

# Callable without a `function NAME(` definition in the document.
_BROWSER_GLOBALS = {'alert', 'confirm', 'prompt', 'print', 'open', 'close', 'requestAnimationFrame'}


def undefined_handlers(html):
    """Names invoked from inline event handlers that the document never defines.

    Returns a sorted list of names; empty when the report is sound.
    """
    defined = set(_FUNC_DEF.findall(html)) | _BROWSER_GLOBALS
    used = set(_HANDLER.findall(html))
    return sorted(used - defined)


def unbalanced_braces(js):
    """Line number of the first `}` that closes a block never opened, else None.

    A deliberately narrow check rather than a parser: it catches the one failure that
    actually shipped -- a statement block left without its `function` header, so the
    trailing brace closes nothing -- and needs no JavaScript engine, so it runs
    everywhere rather than skipping when a parser is absent. A skipped check would not
    have caught this.

    String literals, template literals, both comment forms and regex literals are
    skipped, because every one of them can legally contain a brace. Regex detection
    uses the standard heuristic: a `/` starts a literal only where a value cannot
    already have ended, which is after an operator, a comma, or an opening bracket.
    """
    depth = 0
    line = 1
    i = 0
    n = len(js)
    prev = ''                     # last significant character, for the regex test
    while i < n:
        c = js[i]
        if c == '\n':
            line += 1
            i += 1
            continue
        if c == '/' and i + 1 < n and js[i + 1] == '/':
            while i < n and js[i] != '\n':
                i += 1
            continue
        if c == '/' and i + 1 < n and js[i + 1] == '*':
            i += 2
            while i + 1 < n and not (js[i] == '*' and js[i + 1] == '/'):
                if js[i] == '\n':
                    line += 1
                i += 1
            i += 2
            continue
        if c in '"\'`':
            quote = c
            i += 1
            while i < n and js[i] != quote:
                if js[i] == '\\':
                    i += 1
                elif js[i] == '\n':
                    line += 1
                i += 1
            i += 1
            prev = quote
            continue
        if c == '/' and prev in '(,=:[!&|?{};+-*%~^<>' + '':
            i += 1                # regex literal
            while i < n and js[i] != '/':
                if js[i] == '\\':
                    i += 1
                elif js[i] == '[':
                    while i < n and js[i] != ']':
                        if js[i] == '\\':
                            i += 1
                        i += 1
                elif js[i] == '\n':
                    break         # unterminated: not a regex after all
                i += 1
            i += 1
            prev = '/'
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth < 0:
                return line
        if not c.isspace():
            prev = c
        i += 1
    return None


def check_report(html):
    """Return a list of human-readable problems with the assembled report."""
    problems = []
    missing = undefined_handlers(html)
    if missing:
        problems.append(
            "inline event handlers call JavaScript functions that are never defined: "
            + ", ".join(f"{n}()" for n in missing))

    # One SyntaxError anywhere disables EVERY script in the document, so this is
    # checked per <script> block to name the offending one.
    for idx, body in enumerate(_SCRIPTS.findall(html), start=1):
        bad = unbalanced_braces(body)
        if bad is not None:
            problems.append(
                f"script block {idx} has a closing brace that closes nothing, at its "
                f"line {bad} -- this is a SyntaxError and it disables ALL JavaScript "
                f"in the report, not just the surrounding function")
    return problems
