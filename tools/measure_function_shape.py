"""Reproducible function spans and argument counts, independent of CC analyzers."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess


def function_shapes(source: str, path: str) -> list[dict]:
    """Count physical definition spans and declared parameters, including self.

    Nested functions and methods are functions; class blocks are not additional
    observations. Decorators precede the measured span. *args/**kwargs each
    count once; annotations/defaults do not add parameters.
    """
    rows = []

    class Inventory(ast.NodeVisitor):
        def __init__(self):
            self.names = []

        def visit_ClassDef(self, node):
            self.names.append(node.name)
            self.generic_visit(node)
            self.names.pop()

        def visit_FunctionDef(self, node):
            args = node.args
            rows.append(
                dict(
                    path=path,
                    function=".".join([*self.names, node.name]),
                    line=node.lineno,
                    end_line=node.end_lineno,
                    lines=node.end_lineno - node.lineno + 1,
                    parameters=len(args.posonlyargs)
                    + len(args.args)
                    + len(args.kwonlyargs)
                    + bool(args.vararg)
                    + bool(args.kwarg),
                )
            )
            self.names.append(node.name)
            self.generic_visit(node)
            self.names.pop()

        visit_AsyncFunctionDef = visit_FunctionDef

    Inventory().visit(ast.parse(source, filename=path))
    return rows


def measure(root: Path, revision: str | None = None) -> dict:
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=root, text=True)

    commit = git("rev-parse", revision or "HEAD").strip()
    if revision:
        paths = git("ls-tree", "-r", "--name-only", commit, "pyar").splitlines()
        sources = [
            (path, git("show", f"{commit}:{path}"))
            for path in paths
            if path.endswith(".py")
        ]
    else:
        sources = [
            (str(path.relative_to(root)), path.read_text())
            for path in sorted((root / "pyar").rglob("*.py"))
        ]
    rows = [row for path, source in sources for row in function_shapes(source, path)]
    digest = hashlib.sha256()
    for path, source in sources:
        digest.update(path.encode() + b"\0" + source.encode() + b"\0")
    return dict(
        schema_version=1,
        source_commit=commit,
        working_tree=revision is None,
        source_sha256=digest.hexdigest(),
        definitions="physical def span; all declared arguments including self; nested functions included; no class blocks",
        functions=rows,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--revision")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(measure(args.root.resolve(), args.revision), indent=2) + "\n"
    if args.output:
        args.output.write_text(report)
    else:
        print(report, end="")


if __name__ == "__main__":
    main()
