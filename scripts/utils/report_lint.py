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

import collections
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

# `onclick="toggleNode('node_abc')"` -- a handler whose first argument is a quoted
# string literal. Attributes here are double-quoted, so the literal inside is
# single-quoted.
_HANDLER_ARG = re.compile(
    r'\bon(?:click|keyup|keydown|keypress|change|input|submit|focus|blur)\s*=\s*'
    r'"\s*([A-Za-z_$][\w$]*)\s*\(\s*\'([^\']*)\'')
_ELEMENT_ID = re.compile(r'\bid\s*=\s*"([^"]*)"')

# Handlers whose first argument IS an element id. `dangling_handler_targets` infers
# this per document, so the set is not a list to keep up to date; it exists for the
# one shape inference cannot read -- a document in which EVERY one of a handler's
# targets is missing offers nothing to recognise the contract from, so the check
# would pass the most broken report of all. `toggleNode` is named because that is
# the failure this was written for.
_ID_ARG_HANDLERS = {'toggleNode'}


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


def dangling_handler_targets(html):
    """(handler, id) pairs where an inline handler names an element that is not there.

    `onclick="toggleNode(\'node_x\')"` with no `id="node_x"` anywhere is a click that
    throws a TypeError into the console and changes nothing on the page, while the
    report renders, publishes and reports success -- the same silent shape as a
    handler wired to a function nothing defines, which is why it belongs beside it.

    It shipped on every taxonomy node whose lineage stopped above species rank:
    `viz/taxonomy.py` put the handler on EVERY node header while the
    `<div class="node-children" id=...>` it targets was emitted only for nodes with
    children or species-rank genomes.

    A handler counts as id-taking when at least one of its literals is an id in this
    document, which keeps the rule general rather than a list of names to maintain;
    `_ID_ARG_HANDLERS` covers the all-missing case that inference cannot see. Ids are
    read from the whole document rather than from rendered elements only, so a
    handler pointing at something built at runtime is not flagged.
    """
    ids = set(_ELEMENT_ID.findall(html))
    literals = {}
    for name, arg in _HANDLER_ARG.findall(html):
        literals.setdefault(name, set()).add(arg)
    out = []
    for name, args in sorted(literals.items()):
        if name not in _ID_ARG_HANDLERS and not (args & ids):
            continue            # its argument is some other string, not an element id
        out += [(name, a) for a in sorted(args - ids)]
    return out


def duplicate_definitions(html):
    """(name, count) for every handler this document defines more than once.

    JavaScript lets a later `function f()` replace an earlier one without a word, so
    two definitions are not an error -- they are a coin toss decided by assembly
    order. `toggleNode` was defined twice for months: once in REPORT_JS and once in
    the taxonomy tree\'s own <script>. `visualize_results.py` emits the tree before
    REPORT_JS, so REPORT_JS won, and the guarded copy -- the one that checked its
    element existed -- was the one being discarded. Both copies looked fine on their
    own; only their number was the bug.

    Counted inside executable <script> blocks only, so a `function` that happens to
    appear in a JSON island or in prose is data, not a definition.
    """
    counts = collections.Counter()
    for body in _SCRIPTS.findall(html):
        counts.update(_FUNC_DEF.findall(body))
    return sorted((n, c) for n, c in counts.items() if c > 1)


def check_report(html):
    """Return a list of human-readable problems with the assembled report."""
    problems = []
    missing = undefined_handlers(html)
    if missing:
        problems.append(
            "inline event handlers call JavaScript functions that are never defined: "
            + ", ".join(f"{n}()" for n in missing))

    dangling = dangling_handler_targets(html)
    if dangling:
        shown = ", ".join(f"{n}('{t}')" for n, t in dangling[:5])
        more = f" (+{len(dangling) - 5} more)" if len(dangling) > 5 else ""
        problems.append(
            "inline event handlers name element ids that the document never defines: "
            f"{shown}{more} -- clicking one raises a TypeError into the console and "
            "leaves the page unchanged")

    for name, count in duplicate_definitions(html):
        problems.append(
            f"{name}() is defined {count} times in this report. A later definition "
            f"replaces an earlier one silently, so which copy runs is decided by the "
            f"order the page is assembled in -- delete all but one")

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
