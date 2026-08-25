"""Report calls to names that are never defined, imported, or bound in scope.

Catches the class of bug where a function is used but its import was dropped —
the failure otherwise only surfaces on a rarely-taken branch at runtime.

Usage: python tests/check_undefined.py [scripts_dir]
Prints nothing and exits 0 when clean.
"""
import ast, builtins, sys, pathlib

def check(path):
    tree = ast.parse(open(path).read(), path)
    defined = set(dir(builtins))
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            defined.add(node.name)
        elif isinstance(node, ast.Import):
            for a in node.names: defined.add((a.asname or a.name).split('.')[0])
        elif isinstance(node, ast.ImportFrom):
            for a in node.names: defined.add(a.asname or a.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            defined.add(node.id)
        elif isinstance(node, (ast.arg,)):
            defined.add(node.arg)
        elif isinstance(node, ast.ExceptHandler) and node.name:
            defined.add(node.name)
        elif isinstance(node, ast.comprehension):
            for n in ast.walk(node.target):
                if isinstance(n, ast.Name): defined.add(n.id)
    missing = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id not in defined:
                missing.setdefault(node.func.id, []).append(node.lineno)
    return missing

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else 'scripts')
for p in sorted(root.rglob('*.py')):
    if '__pycache__' in str(p): continue
    m = check(p)
    if m:
        print(f"{p}:")
        for name, lines in sorted(m.items()):
            print(f"    {name}  (lines {lines[:6]})")
