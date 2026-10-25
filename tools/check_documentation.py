"""Validate documented install metadata and CLI flags against this checkout."""

from __future__ import annotations

import re
import shlex
import subprocess
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 fallback provided by PyAR dependencies
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SUBCOMMANDS = {
    "aggregate",
    "grow",
    "microsolvate",
    "solvate",
    "optimize",
    "energy",
    "engrad",
    "gradient",
    "engra",
    "frequency",
    "react",
    "scan-bond",
    "neb",
    "ts",
    "irc",
    "conformer",
    "clustering",
    "energies",
    "compare",
    "identify",
    "orient",
    "deduplicate",
    "select",
    "split",
    "trace",
    "run",
    "config",
    "profile",
    "inspect",
    "info",
    "backends",
    "doctor",
    "install",
}
NESTED_HELP = {
    ("pyar", "install", "status"),
    ("pyar", "config", "validate"),
    ("pyar", "profile", "create"),
    ("pyar", "profile", "list"),
}


def _shell_blocks(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        rst = re.match(r"\s*\.\. code-block:: (bash|console|shell)\s*$", lines[i])
        md = re.match(r"\s*```(?:bash|console|shell)\s*$", lines[i])
        if rst:
            i += 1
            block = []
            while i < len(lines) and (
                not lines[i].strip() or lines[i].startswith("   ")
            ):
                block.append(lines[i].strip())
                i += 1
        elif md:
            i += 1
            block = []
            while i < len(lines) and not lines[i].lstrip().startswith("```"):
                block.append(lines[i].strip())
                i += 1
            i += 1
        else:
            i += 1
            continue
        command = ""
        for line in block:
            if not line or line.startswith("#"):
                continue
            command += (" " if command else "") + line.rstrip("\\").strip()
            if line.endswith("\\"):
                continue
            yield command
            command = ""
        if command:
            yield command


def _python_blocks(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        match = re.match(r"\s*\.\. code-block:: python\s*$", lines[i])
        if not match:
            i += 1
            continue
        i += 1
        block = []
        while i < len(lines) and (not lines[i].strip() or lines[i].startswith("   ")):
            block.append(lines[i].strip())
            i += 1
        yield "\n".join(block).strip()


def _route(tokens):
    executable = tokens[0]
    if executable == "pyar":
        if len(tokens) > 2 and tuple(tokens[:3]) in NESTED_HELP:
            return tokens[:3] + ["--help"]
        if len(tokens) > 1 and tokens[1] in SUBCOMMANDS:
            return tokens[:2] + ["--help"]
        return ["pyar", "--help"]
    if executable == "pyar-cli":
        if len(tokens) > 1 and tokens[1] in SUBCOMMANDS:
            return tokens[:2] + ["--help"]
        return ["pyar-cli", "--help"]
    if executable.startswith("pyar-"):
        return [executable, "--help"]
    return None


def _options(help_text):
    return set(re.findall(r"(?<![\w-])--?[A-Za-z][A-Za-z0-9-]*", help_text))


def _documented_commands(paths=None):
    paths = paths or [ROOT / "README.md", *sorted(DOCS.glob("*.rst"))]
    for path in paths:
        for command in _shell_blocks(path):
            try:
                tokens = shlex.split(command)
            except ValueError:
                continue
            if tokens:
                yield path, command, tokens


def check(root=ROOT, run=subprocess.run):
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))[
        "project"
    ]
    install_doc = (root / "docs/installation.rst").read_text(encoding="utf-8")
    minimum = re.fullmatch(r">=([0-9]+\.[0-9]+)", project["requires-python"])
    if not minimum or f"{minimum.group(1)} or newer" not in install_doc:
        raise ValueError(
            "Installation guide Python support does not match project metadata"
        )
    ci = (root / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    matrix = re.search(r"python-version:\s*\[([^]]+)\]", ci)
    tested_versions = (
        re.findall(r"[\"'](3\.\d+)[\"']", matrix.group(1)) if matrix else []
    )
    classified = sorted(
        re.findall(
            r"Programming Language :: Python :: (3\.\d+)",
            "\n".join(project["classifiers"]),
        )
    )
    if not tested_versions or sorted(tested_versions) != classified:
        raise ValueError("Python classifiers and CI support matrix must agree")
    support_range = f"{tested_versions[0]}–{tested_versions[-1]}"
    if any(
        support_range not in text
        for text in (install_doc, (root / "README.md").read_text(encoding="utf-8"))
    ):
        raise ValueError(
            "README and installation guide must state the classified, CI-tested Python range"
        )
    extras = set(project["optional-dependencies"])
    for extra in extras:
        if f"``{extra}``" not in install_doc:
            raise ValueError(
                f"Installation guide does not document declared extra [{extra}]"
            )
    groups = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8")).get(
        "dependency-groups", {}
    )
    for group in groups:
        if f"``{group}``" not in install_doc:
            raise ValueError(
                f"Installation guide does not document dependency group [{group}]"
            )
    scripts = set(project["scripts"])
    cli_reference = (root / "docs/cli_reference.rst").read_text(encoding="utf-8")
    documented_scripts = set(re.findall(r"``(pyar(?:-[a-z0-9-]+)?)``", cli_reference))
    if scripts != documented_scripts:
        missing = sorted(scripts - documented_scripts)
        stale = sorted(documented_scripts - scripts)
        raise ValueError(
            f"CLI reference differs from project.scripts; missing={missing}, stale={stale}"
        )
    errors = []
    api_import_block = next(
        (
            block
            for block in _python_blocks(root / "docs/architecture.rst")
            if "from pyar import Molecule" in block
        ),
        None,
    )
    if api_import_block is None:
        errors.append(
            "Architecture guide is missing its declared public-import example"
        )
    else:
        imports = run(
            [sys.executable, "-c", api_import_block],
            cwd=root,
            text=True,
            capture_output=True,
        )
        if imports.returncode:
            errors.append(
                f"Architecture public imports no longer work: {imports.stderr.strip()}"
            )
    help_cache = {}
    examples = list(
        _documented_commands(
            [root / "README.md", *sorted((root / "docs").rglob("*.rst"))]
        )
    )
    for path, command, tokens in examples:
        for group in re.findall(r"(?:pyar-chem|\.)\[([^]]+)\]", command):
            for extra in group.split(","):
                if extra not in extras:
                    errors.append(
                        f"{path.relative_to(root)}: undeclared extra [{extra}]"
                    )
        route = _route(tokens)
        if route is None:
            continue
        route_key = tuple(route)
        if route_key not in help_cache:
            result = run(route, cwd=root, text=True, capture_output=True)
            help_cache[route_key] = result
        result = help_cache[route_key]
        if result.returncode:
            errors.append(
                f"{path.relative_to(root)}: cannot inspect {' '.join(route)}: {result.stderr.strip()}"
            )
            continue
        documented = {t.split("=", 1)[0] for t in tokens if re.match(r"^-+[A-Za-z]", t)}
        unknown = sorted(
            flag
            for flag in documented
            if flag not in _options(result.stdout + result.stderr)
        )
        if unknown:
            errors.append(
                f"{path.relative_to(root)}: {' '.join(route[:-1])} has no documented option(s) {', '.join(unknown)}"
            )
    preflight = run(
        ["pyar", "conformer", "CCO", "--check"],
        cwd=root,
        text=True,
        capture_output=True,
    )
    if preflight.returncode or "no calculations were performed" not in preflight.stdout:
        errors.append(
            "Documented first-run fixture preflight did not complete read-only"
        )
    if errors:
        raise ValueError("Documentation contract failures:\n- " + "\n- ".join(errors))
    cli_examples = sum(1 for _, _, tokens in examples if _route(tokens))
    print(
        f"Documentation contract passed: {cli_examples} CLI examples, package metadata, extras, scripts and first-run preflight."
    )


if __name__ == "__main__":
    check()
