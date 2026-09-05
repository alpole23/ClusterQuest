#!/usr/bin/env python3
"""Verify each module's declared scripts-version dependencies still cover its imports.

Modules invoke their Python scripts as `python ${projectDir}/scripts/foo.py` —
interpolated paths, not declared `path` inputs — so Nextflow's task hash never sees
them and `-resume` would reuse output produced by code that has since changed. Each
process therefore declares `Utils.scriptsHash(projectDir, [...])`.

Nextflow hashes the *unevaluated* script source plus the input values, never the
rendered text, so the digest only invalidates anything when it arrives as a `val`
input — passed from the call site in subworkflows/. Older modules still embed it
as a comment in the script block, where it has no effect; both spellings are
accepted here so the two can coexist while modules are converted.

Those lists are written by hand, so they can drift when a script gains an import.
This checks that every module a process actually reaches is covered by its declared
list, and that the list does not name anything that does not exist.

Exit 0 if consistent, 1 otherwise.
"""
import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / 'scripts'
PKGS = {'utils', 'viz', 'analysis', 'clustering', 'genome', 'taxonomy', 'phylogeny'}


def transitive_imports(rel, seen=None):
    """Every scripts/-relative .py file reachable from `rel` by import."""
    seen = set() if seen is None else seen
    if rel in seen:
        return seen
    seen.add(rel)
    path = SCRIPTS / rel
    if not path.exists():
        return seen
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError:
        return seen
    for node in ast.walk(tree):
        mod = None
        if isinstance(node, ast.ImportFrom) and node.module:
            mod = node.module
        elif isinstance(node, ast.Import):
            mod = node.names[0].name
        if not mod or mod.split('.')[0] not in PKGS:
            continue
        as_file = SCRIPTS / (mod.replace('.', '/') + '.py')
        as_pkg = SCRIPTS / mod.replace('.', '/')
        if as_file.exists():
            transitive_imports(str(as_file.relative_to(SCRIPTS)), seen)
        elif as_pkg.is_dir():
            for f in sorted(as_pkg.rglob('*.py')):
                transitive_imports(str(f.relative_to(SCRIPTS)), seen)
    return seen


def covered(rel, declared):
    """True if `rel` falls under any declared file or package directory."""
    return any(rel == d or rel.startswith(d.rstrip('/') + '/') for d in declared)


def call_site_hashes():
    """process name -> declared list, for digests passed in as a `val` input."""
    out = {}
    for nf in sorted((ROOT / 'subworkflows').rglob('*.nf')) + \
              sorted((ROOT / 'workflows').rglob('*.nf')):
        text = nf.read_text()
        for m in re.finditer(r'\b([A-Z][A-Z0-9_]+)\s*\(', text):
            # Scan to the matching close paren rather than using a nested-paren
            # regex: call sites like GENBANK_TO_FASTA(x.collate(batchSize()), ...)
            # nest arbitrarily deep.
            depth, i = 0, m.end() - 1
            while i < len(text):
                if text[i] == '(':
                    depth += 1
                elif text[i] == ')':
                    depth -= 1
                    if depth == 0:
                        break
                i += 1
            call = text[m.end():i]
            h = re.search(r'scriptsHash\(projectDir,\s*\[([^\]]*)\]\)', call, re.S)
            if h:
                out.setdefault(m.group(1), []).extend(re.findall(r"'([^']+)'", h.group(1)))
    return out


def main():
    problems = []
    checked = 0
    call_sites = call_site_hashes()
    for nf in sorted((ROOT / 'modules').rglob('*.nf')):
        text = nf.read_text()
        blocks = [m.start() for m in re.finditer(r'^process ([A-Z_]+)', text, re.M)] + [len(text)]
        names = [re.match(r'process ([A-Z_]+)', text[i:]).group(1) for i in blocks[:-1]]
        for i, proc in enumerate(names):
            blk = text[blocks[i]:blocks[i + 1]]
            invoked = sorted(set(re.findall(r'scripts/([A-Za-z0-9_/]+\.py)', blk)))
            decl_m = re.search(r'scriptsHash\(projectDir,\s*\[([^\]]*)\]\)', blk)
            declared = re.findall(r"'([^']+)'", decl_m.group(1)) if decl_m else []
            if not declared and proc in call_sites:
                declared = call_sites[proc]
                decl_m = True

            if invoked and not decl_m:
                problems.append(f"{nf.name}:{proc} runs {invoked} but declares no scripts-version")
                continue
            if decl_m and not invoked:
                problems.append(f"{nf.name}:{proc} declares a scripts-version but runs no script")
                continue
            if not invoked:
                continue
            checked += 1

            for d in declared:
                if not (SCRIPTS / d).exists():
                    problems.append(f"{nf.name}:{proc} declares '{d}' which does not exist")

            needed = set()
            for entry in invoked:
                needed |= transitive_imports(entry)
            missing = sorted(r for r in needed if not covered(r, declared))
            if missing:
                problems.append(
                    f"{nf.name}:{proc} imports {missing} not covered by {declared}")

    if problems:
        print("script dependency declarations are out of date:")
        for p in problems:
            print(f"  - {p}")
        print("\nUpdate the scripts-version list in the affected process(es).")
        return 1
    print(f"script dependency declarations OK ({checked} processes checked)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
