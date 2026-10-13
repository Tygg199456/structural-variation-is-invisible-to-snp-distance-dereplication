#!/usr/bin/env python3
"""
Guard against defects that snakemake surfaces late, or not at all.

Run as a pre-commit hook, or from the top of the Snakefile:

    import subprocess, sys
    subprocess.run([sys.executable, "scripts/_lint_scripts.py"], check=True)

CHECK 1 -- `from __future__` in a script: target
    Snakemake copies the script to .snakemake/scripts/tmpXXXX.py and inserts
    a preamble at the top of the file. A `from __future__` import is then no
    longer the first statement, and Python refuses to compile:

        SyntaxError: from __future__ imports must occur at the beginning
                     of the file

    This fires at runtime, often after a long upstream step has already run.
    Removal is safe on Python >= 3.9 for builtin generics (list[str],
    dict[str, str]); only PEP 604 unions (X | None) in annotations need
    >= 3.10.

CHECK 2 -- PEP 604 unions, which is what makes removal unsafe below 3.10.
    Reported as a warning with the running interpreter's version so the
    trade-off is visible rather than assumed.

CHECK 3 -- stubs still raising NotImplementedError, listed so the remaining
    work is explicit rather than discovered when a rule fails at 3am.
"""
import ast
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
PEP604 = re.compile(r":\s*[A-Za-z_][\w.\[\], ]*\s*\|\s*(None|[A-Z])")

errors, warnings, stubs = [], [], []

for path in sorted(HERE.glob("*.py")):
    if path.name.startswith("_"):
        continue
    text = path.read_text()

    if "from __future__" in text:
        for i, line in enumerate(text.splitlines(), 1):
            if line.strip().startswith("from __future__"):
                errors.append(
                    f"{path.name}:{i}: `from __future__` breaks the snakemake "
                    f"script wrapper. Remove it."
                )

    if PEP604.search(text) and sys.version_info < (3, 10):
        warnings.append(
            f"{path.name}: uses PEP 604 unions (X | None) and this "
            f"interpreter is {sys.version_info.major}.{sys.version_info.minor}"
        )

    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        errors.append(f"{path.name}: does not parse: {exc}")
        continue

    for node in ast.walk(tree):
        if (isinstance(node, ast.Raise) and isinstance(node.exc, ast.Call)
                and getattr(node.exc.func, "id", "") == "NotImplementedError"):
            stubs.append(path.name)
            break

if stubs:
    print(f"[lint] {len(stubs)} stub(s) still unimplemented:")
    for s in sorted(set(stubs)):
        print(f"       {s}")

for w in warnings:
    print(f"[lint] WARNING {w}", file=sys.stderr)

if errors:
    print(file=sys.stderr)
    for e in errors:
        print(f"[lint] ERROR {e}", file=sys.stderr)
    sys.exit(1)

print(f"[lint] OK -- no `from __future__` in any script: target")
