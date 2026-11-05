#!/usr/bin/env python3
"""Shared report helpers must accept the shapes their producers actually emit.

Both failures below shipped in one commit, passed 29 tests, a pyflakes pass, an
undefined-name scan and a `-preview`, and then raised on the FIRST real run. Neither is
reachable without data: `check_undefined.py` sees that a name exists, not what type it
holds, and a helper shared by two producers is only exercised when both run.

1. `utils.gene_diagram._domain_label` read `gene['domains']` as a `;`-joined string.
   `utils.antismash_parser.extract_genes_from_region` sets it to antiSMASH's
   `sec_met_domain` qualifier, which is a **list**:
       AttributeError: 'list' object has no attribute 'strip'
   The consensus diagram sets no domain field at all, so the function must also return
   '' rather than raise. Three producers, three shapes, one consumer.

2. `analysis/biosynthetic_profile.py` sorted family ids with `key=lambda f: (len(f), f)`
   -- a key for string ids -- and joined them with `';'.join(...)`. `family.id` comes
   from SQLite as an **int**:
       TypeError: object of type 'int' has no len()

So this file pins the CONTRACT between producer and consumer: the shapes each producer
emits, and that the consumer survives all of them. When a producer changes shape, this
fails here instead of three stages into a resumed run.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))

FAILURES = []


def _code_only(src):
    """`src` with comments and docstrings removed, so a scan cannot match prose.

    The first version of this file matched its own explanation: the regex below rejects
    `key=lambda f: (len(f), f)`, and the fixed module describes that bug in a comment
    using those exact characters, so the check failed on correct code. A checker that
    reads source has to read CODE.

    Docstrings go too, since they carry the same prose. The test for one is "a string
    that opens a logical line" -- but ONLY outside brackets. Inside a bracket a line
    break is an NL token, so a dict key on its own line (`'same_core_gcfs':`, which is
    exactly what this file scans for) looks identical to a docstring. The first version
    dropped it and then reported the key as missing. Hence `depth`.

    Tokens are joined with single spaces, so every scan built on this must tolerate
    whitespace that the source does not contain.
    """
    import io
    import tokenize
    out = []
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except (tokenize.TokenError, IndentationError):
        return src                       # unparseable: scan the raw text rather than skip
    depth = 0
    opens_line = True                    # start of file opens a logical line
    for tok in toks:
        if tok.type == tokenize.COMMENT:
            continue
        if tok.type in (tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT):
            opens_line = True
            out.append(tok.string)
            continue
        if tok.type == tokenize.NL:
            out.append(tok.string)
            continue                     # a blank/continuation line changes nothing
        if tok.type == tokenize.STRING and depth == 0 and opens_line:
            opens_line = False
            continue                      # a docstring
        if tok.type == tokenize.OP:
            if tok.string in '([{':
                depth += 1
            elif tok.string in ')]}':
                depth = max(0, depth - 1)
        out.append(tok.string)
        opens_line = False
    return ' '.join(out)


def check(label, fn, want):
    """Compare fn() with want; a raise is a failure, not a traceback.

    `fn` is a thunk rather than a value because the bug this file exists for RAISES.
    Passing `_domain_label(...)` directly would abort the run at the first bad shape and
    report a traceback instead of naming every shape that is broken.
    """
    try:
        got = fn()
    except Exception as exc:                       # noqa: BLE001 -- that IS the failure
        FAILURES.append(f"{label}: raised {type(exc).__name__}: {exc}")
        return
    if got != want:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")


def check_domain_label():
    from utils.gene_diagram import _domain_label, _is_generic, mark_for

    # -- the three shapes the three producers actually emit -------------------
    # extract_gcf_representatives.py, via antismash_parser: clusterhmmer_hits carry
    # `hmm` (the domain NAME) and `identifier` (the accession); `domains` is a LIST.
    check('clusterhmmer preferred over domains',
          lambda: _domain_label({'clusterhmmer_hits': [{'hmm': 'PEP_mutase', 'score': 300},
                                               {'hmm': 'HMGL-like', 'score': 90}],
                         'domains': ['ATP-grasp (E-value: 1e-9)']}),
          'PEP_mutase')
    check('domains as a LIST (the AttributeError)',
          lambda: _domain_label({'domains': ['PEP_mutase']}),
          'PEP_mutase')
    # report_sections.py's consensus diagram sets neither field.
    check('no domain field at all returns empty',
          lambda: _domain_label({'product': 'x', 'locus_tag': 'y'}),
          '')

    # -- antiSMASH's sec_met_domain carries its evidence inline ---------------
    check('sec_met_domain trimmed to the name',
          lambda: _domain_label({'domains': [
              'ATP-grasp (E-value: 7.5e-08, bitscore: 23.1, seeds: 20, '
              'tool: rule-based-clusters)']}),
          'ATP-grasp')
    check('a name containing a space survives trimming',
          lambda: _domain_label({'domains': ['ATP grasp fold (E-value: 1e-9)']}),
          'ATP grasp fold')

    # -- degenerate inputs must not raise ------------------------------------
    check('empty list',
          lambda: _domain_label({'domains': []}),
          '')
    check('joined string form still accepted',
          lambda: _domain_label({'domains': 'PEP_mutase;HMGL-like'}),
          'PEP_mutase')
    check('blank entries skipped',
          lambda: _domain_label({'domains': ['', '   ', 'Fe-ADH']}),
          'Fe-ADH')
    check('em-dash placeholder is not a domain name',
          lambda: _domain_label({'domains': ['—']}),
          '')
    check('evidence with no name falls through',
          lambda: _domain_label({'domains': [' (E-value: 1e-9)', 'Fe-ADH']}),
          'Fe-ADH')
    check('unnamed hit falls through to domains',
          lambda: _domain_label({'clusterhmmer_hits': [{'hmm': '', 'description': '', 'score': 5}],
                         'domains': ['Aminotran_1_2']}),
          'Aminotran_1_2')
    check('unsortable score does not raise',
          lambda: _domain_label({'clusterhmmer_hits': [{'hmm': 'Ppd', 'score': None}]}),
          'Ppd')

    # -- the generic-product test the tooltip branches on ---------------------
    check('generic product recognised',
          lambda: _is_generic('hypothetical protein'),
          True)
    check('generic product is case- and space-insensitive',
          lambda: _is_generic('  HYPOTHETICAL PROTEIN '),
          True)
    check('a real product is not generic',
          lambda: _is_generic('phosphoenolpyruvate mutase'),
          False)
    check('None product is not generic',
          lambda: _is_generic(None),
          False)

    # -- mark_for: accessions arrive versioned, from a generator -------------
    check('versioned accession matches',
          lambda: mark_for(['PF13714.9']),
          'pepM')
    check('generator input accepted',
          lambda: mark_for(h['identifier'] for h in [{'identifier': 'PF00682'}]),
          'synthase')
    check('non-marker returns empty',
          lambda: mark_for(['PF99999']),
          '')
    check('empty input returns empty',
          lambda: mark_for([]),
          '')
    check('None entries skipped',
          lambda: mark_for([None, '', 'PF02775']),
          'Ppd')


def check_profile_id_types():
    """`same_core_gcfs` must be built from INT family ids and emitted as strings."""
    src = _code_only(
        (ROOT / 'scripts' / 'analysis' / 'biosynthetic_profile.py').read_text(
            encoding='utf-8'))

    # The exact key that raised. Reject it by shape, not by remembering not to retype it.
    if re.search(r'key\s*=\s*lambda\s+\w+\s*:\s*\(\s*len\s*\(', src):
        FAILURES.append(
            "biosynthetic_profile.py sorts family ids with a `(len(x), x)` key. "
            "family.id is an INTEGER from SQLite, so len() raises "
            "`TypeError: object of type 'int' has no len()`. Sort numerically.")

    # ';'.join over ints raises too, so the ids must be stringified at the join.
    m = re.search(r"'same_core_gcfs'\s*:\s*(?P<expr>[^,\n]+)", src)
    if not m:
        FAILURES.append("biosynthetic_profile.py no longer emits 'same_core_gcfs'")
    elif not re.search(r'\bstr\s*\(', m.group('expr')):
        FAILURES.append(
            f"biosynthetic_profile.py joins same_core_gcfs without str(): "
            f"{m.group('expr').strip()} -- ';'.join() over ints raises TypeError")

    # The report prefixes 'GCF-' onto each id, so the producer must emit BARE ids.
    sections = (ROOT / 'scripts' / 'viz' / 'report_sections.py').read_text(
        encoding='utf-8')
    if not re.search(r"""f['"]GCF-\{g\}['"]""", sections):
        FAILURES.append(
            "report_sections._core_cell no longer prefixes 'GCF-'; if it now expects "
            "pre-formatted ids, biosynthetic_profile.py must emit them that way")


def main():
    check_domain_label()
    check_profile_id_types()
    if FAILURES:
        for f in FAILURES:
            print(f"FAIL: {f}", file=sys.stderr)
        return 1
    print("ok: report helpers accept every shape their producers emit")
    return 0


if __name__ == '__main__':
    sys.exit(main())
