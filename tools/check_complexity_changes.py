"""Compare changed first-party PyAR functions with their Git baseline."""

from __future__ import annotations

import argparse
import ast
import hashlib
import io
import json
import re
import subprocess
import sys
import tarfile
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.measure_complexity import (
    MeasurementError,
    _named_functions,
    measure_complexity,
)


TOOL_VERSION = "1"
_EXCEPTION_REASON = re.compile(r"^Complexity-Increase-Reason:\s*(\S.*)$", re.I | re.M)
_EXCEPTION_OWNER = re.compile(r"^Complexity-Increase-Owner:\s*(\S.*)$", re.I | re.M)


class ComplexityChangeError(ValueError):
    """Raised when a comparable source snapshot cannot be measured."""


def _git(repo: Path, *args: str, text: bool = True):
    return subprocess.run(
        ["git", *args], cwd=repo, check=True, capture_output=True, text=text
    ).stdout


def _verify_commit(repo: Path, revision: str) -> str:
    try:
        resolved = _git(repo, "rev-parse", "--verify", f"{revision}^{{commit}}").strip()
    except subprocess.CalledProcessError as exc:
        raise ComplexityChangeError(f"Cannot resolve Git commit {revision!r}") from exc
    if resolved != revision:
        raise ComplexityChangeError(
            f"Revision is not an exact commit SHA: {revision!r}"
        )
    return resolved


def _extract_python_archive(repo: Path, revision: str, destination: Path) -> None:
    archive_bytes = _git(repo, "archive", "--format=tar", revision, "pyar", text=False)
    destination.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:") as archive:
            for member in archive:
                relative = Path(member.name)
                if relative.is_absolute() or ".." in relative.parts:
                    raise ComplexityChangeError("Git archive contains an unsafe path")
                if not member.isfile() or relative.suffix != ".py":
                    continue
                target = destination / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                content = archive.extractfile(member)
                if content is None:
                    raise ComplexityChangeError(
                        f"Cannot read archived source {member.name}"
                    )
                target.write_bytes(content.read())
    except tarfile.TarError as exc:
        raise ComplexityChangeError(
            f"Cannot read source archive for {revision}"
        ) from exc
    if not (destination / "pyar").is_dir():
        raise ComplexityChangeError(f"Commit {revision} has no first-party pyar source")


def _changed_python_files(repo: Path, base_sha: str, head_sha: str) -> list[str]:
    try:
        output = _git(
            repo,
            "diff",
            "--name-only",
            "--diff-filter=ACMRTD",
            base_sha,
            head_sha,
            "--",
            "pyar",
        )
    except subprocess.CalledProcessError as exc:
        raise ComplexityChangeError(
            "Cannot determine changed PyAR source files"
        ) from exc
    return sorted(
        path
        for path in output.splitlines()
        if path.endswith(".py") and "AIMNet2" not in Path(path).parts
    )


def _function_inventory(source_root: Path, rows: list[dict]) -> list[dict]:
    complexity = {
        (row["file"], row["qualified_name"], row["line"]): row["cc_radon"]
        for row in rows
    }
    result = []
    package = source_root / "pyar"
    for path in sorted(package.rglob("*.py")):
        if (
            "AIMNet2" in path.relative_to(source_root).parts
            or "__pycache__" in path.parts
        ):
            continue
        relative = path.relative_to(source_root).as_posix()
        source = path.read_text(encoding="utf-8")
        module = ast.parse(source, filename=relative)
        lines = source.splitlines(keepends=True)
        for node, qualified_name in _named_functions(module):
            start = min(
                [decorator.lineno for decorator in node.decorator_list] or [node.lineno]
            )
            body = "".join(lines[start - 1 : node.end_lineno])
            key = (relative, qualified_name, node.lineno)
            if key not in complexity:
                raise ComplexityChangeError(
                    f"No Radon metric for {relative}:{node.lineno}:{qualified_name}"
                )
            result.append(
                {
                    "file": relative,
                    "qualified_name": qualified_name,
                    "line": node.lineno,
                    "end_line": node.end_lineno,
                    "source_sha256": hashlib.sha256(body.encode("utf-8")).hexdigest(),
                    "cc_radon": complexity[key],
                }
            )
    return result


def _changed_functions(
    before: list[dict], after: list[dict], changed_files: set[str]
) -> dict:
    by_snapshot = []
    for rows in (before, after):
        grouped = defaultdict(list)
        for row in rows:
            if row["file"] in changed_files:
                grouped[(row["file"], row["qualified_name"])].append(row)
        for group in grouped.values():
            group.sort(key=lambda row: (row["line"], row["end_line"]))
        by_snapshot.append(grouped)

    before_by_key, after_by_key = by_snapshot
    modified, added, removed = [], [], []
    for key in sorted(set(before_by_key) | set(after_by_key)):
        old = before_by_key.get(key, [])
        new = after_by_key.get(key, [])
        old_hashes = Counter(row["source_sha256"] for row in old)
        new_hashes = Counter(row["source_sha256"] for row in new)
        unchanged_hashes = old_hashes & new_hashes
        old_unmatched = old_hashes - unchanged_hashes
        new_unmatched = new_hashes - unchanged_hashes

        def unmatched_rows(rows, unmatched_hashes):
            remaining = []
            for row in rows:
                digest = row["source_sha256"]
                if unmatched_hashes[digest]:
                    remaining.append(row)
                    unmatched_hashes[digest] -= 1
            return remaining

        old_remaining = unmatched_rows(old, old_unmatched)
        new_remaining = unmatched_rows(new, new_unmatched)
        old_remaining.sort(key=lambda row: row["line"])
        new_remaining.sort(key=lambda row: row["line"])
        paired = min(len(old_remaining), len(new_remaining))
        for index in range(paired):
            modified.append(
                {"before": old_remaining[index], "after": new_remaining[index]}
            )
        removed.extend(old_remaining[paired:])
        added.extend(new_remaining[paired:])
    return {"modified": modified, "added": added, "removed": removed}


def _decision_totals(rows: list[dict], changed_files: set[str]) -> dict:
    totals = {path: 0 for path in sorted(changed_files)}
    for row in rows:
        if row["file"] in totals:
            totals[row["file"]] += row["cc_radon"] - 1
    return totals


def _exception_from_event(path: Path | None) -> dict:
    if path is None:
        return {"reason": None, "owner": None}
    try:
        event = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ComplexityChangeError(f"Cannot read GitHub event JSON: {path}") from exc
    body = event.get("pull_request", {}).get("body") or ""
    reason = _EXCEPTION_REASON.search(body)
    owner = _EXCEPTION_OWNER.search(body)
    return {
        "reason": reason.group(1).strip() if reason else None,
        "owner": owner.group(1).strip() if owner else None,
    }


def analyze_complexity_change(
    repo: str | Path,
    base_sha: str,
    head_sha: str,
    event_json: str | Path | None = None,
    enforce: bool = False,
) -> dict:
    """Compare exact Git snapshots and report decision growth in changed files."""
    repository = Path(repo).resolve()
    base_sha = _verify_commit(repository, base_sha)
    head_sha = _verify_commit(repository, head_sha)
    changed_files = set(_changed_python_files(repository, base_sha, head_sha))
    exception = _exception_from_event(Path(event_json) if event_json else None)

    with (
        tempfile.TemporaryDirectory(prefix="pyar-complexity-base-") as base_dir,
        tempfile.TemporaryDirectory(prefix="pyar-complexity-head-") as head_dir,
    ):
        base_root, head_root = Path(base_dir), Path(head_dir)
        _extract_python_archive(repository, base_sha, base_root)
        _extract_python_archive(repository, head_sha, head_root)
        base_rows, base_manifest = measure_complexity(base_root)
        head_rows, head_manifest = measure_complexity(head_root)
        base_inventory = _function_inventory(base_root, base_rows)
        head_inventory = _function_inventory(head_root, head_rows)

    base_totals = _decision_totals(base_inventory, changed_files)
    head_totals = _decision_totals(head_inventory, changed_files)
    decision_delta = sum(head_totals.values()) - sum(base_totals.values())
    function_changes = _changed_functions(base_inventory, head_inventory, changed_files)
    exception_supplied = bool(exception["reason"] and exception["owner"])
    if decision_delta <= 0:
        status = "pass"
    elif enforce and exception_supplied:
        status = "pass_with_review_note"
    elif enforce:
        status = "fail"
    else:
        status = "report_only_increase"

    return {
        "schema_version": 1,
        "tool": "pyar-changed-code-complexity-ratchet",
        "tool_version": TOOL_VERSION,
        "gate_tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "base_sha": base_sha,
        "head_sha": head_sha,
        "metric": "sum of Radon function-level (CC - 1) across all functions in changed first-party PyAR Python files",
        "radon_version": head_manifest["radon_version"],
        "python_version": head_manifest["python_version"],
        "measurement_tool_sha256": head_manifest["measurement_tool_sha256"],
        "changed_files": sorted(changed_files),
        "base_source_tree_sha256": base_manifest["source_tree_sha256"],
        "head_source_tree_sha256": head_manifest["source_tree_sha256"],
        "base_function_count": len(base_inventory),
        "head_function_count": len(head_inventory),
        "base_decisions_by_file": base_totals,
        "head_decisions_by_file": head_totals,
        "base_decision_total": sum(base_totals.values()),
        "head_decision_total": sum(head_totals.values()),
        "decision_delta": decision_delta,
        "changed_functions": function_changes,
        "coverage_status": "unavailable_for_this_function_level_gate; no coverage or CRAP score is inferred",
        "exception": exception,
        "enforced": enforce,
        "status": status,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--base", required=True)
    parser.add_argument("--head", required=True)
    parser.add_argument("--event-json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--enforce", action="store_true")
    arguments = parser.parse_args(argv)
    try:
        report = analyze_complexity_change(
            arguments.repo,
            arguments.base,
            arguments.head,
            arguments.event_json,
            enforce=arguments.enforce,
        )
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    except (
        ComplexityChangeError,
        MeasurementError,
        OSError,
        subprocess.CalledProcessError,
        SyntaxError,
        ValueError,
    ) as exc:
        parser.error(str(exc))
    print(
        f"Changed-file decision delta: {report['decision_delta']}; "
        f"status={report['status']}; report={arguments.output}"
    )
    return 1 if report["status"] == "fail" else 0


if __name__ == "__main__":
    sys.exit(main())
