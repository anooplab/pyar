"""Create a source-bound Radon cyclomatic-complexity inventory."""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import textwrap
from collections import defaultdict
from importlib.metadata import version
from pathlib import Path

from radon.complexity import cc_visit


TOOL_VERSION = "1"
CSV_FIELDS = (
    "file",
    "qualified_name",
    "line",
    "end_line",
    "cc_radon",
    "cc_measurement_mode",
    "file_sha256",
    "source_sha",
    "source_tree_sha256",
    "radon_version",
    "tool_version",
)


class MeasurementError(ValueError):
    """Raised when a complete source-bound inventory cannot be produced."""


def _production_files(source_root: Path) -> list[Path]:
    package_root = source_root / "pyar"
    if not package_root.is_dir():
        raise MeasurementError(f"Source root has no pyar package: {source_root}")
    return sorted(
        path
        for path in package_root.rglob("*.py")
        if "AIMNet2" not in path.relative_to(source_root).parts
        and "__pycache__" not in path.parts
    )


def _named_functions(node: ast.AST, scope: tuple[str, ...] = ()):
    for child in ast.iter_child_nodes(node):
        if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
            qualified_name = ".".join((*scope, child.name))
            yield child, qualified_name
            yield from _named_functions(child, (*scope, child.name))
        elif isinstance(child, ast.ClassDef):
            yield from _named_functions(child, (*scope, child.name))
        else:
            yield from _named_functions(child, scope)


def _radon_function_blocks(blocks):
    """Yield function blocks, descending through Radon's class/closure blocks."""
    for block in blocks:
        if hasattr(block, "methods"):
            yield from _radon_function_blocks(block.methods)
            yield from _radon_function_blocks(block.inner_classes)
            continue
        yield block
        yield from _radon_function_blocks(block.closures)


def _git_head(source_root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _percentile(values: list[int], fraction: float) -> float | None:
    ordered = sorted(values)
    if not ordered:
        return None
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def _complexity_summary(values: list[int]) -> dict:
    if not values:
        return {"n": 0}
    return {
        "n": len(values),
        "min": min(values),
        "median": statistics.median(values),
        "mean": statistics.mean(values),
        "p75_linear": _percentile(values, 0.75),
        "p90_linear": _percentile(values, 0.90),
        "p95_linear": _percentile(values, 0.95),
        "p99_linear": _percentile(values, 0.99),
        "max": max(values),
        "counts_gt": {
            str(n): sum(value > n for value in values) for n in (10, 15, 20, 30, 50)
        },
        "counts_gte": {
            str(n): sum(value >= n for value in values) for n in (10, 15, 20, 30, 50)
        },
    }


def reconcile_function_inventories(
    previous_rows: list[dict], current_rows: list[dict]
) -> dict:
    """Compare function identities and Radon CC without treating line moves as new code."""
    previous = defaultdict(list)
    current = defaultdict(list)
    for row in previous_rows:
        if row.get("primary", row["file"].startswith("pyar/")):
            previous[(row["file"], row["qualified_name"])].append(row)
    for row in current_rows:
        if row.get("primary", row["file"].startswith("pyar/")):
            current[(row["file"], row["qualified_name"])].append(row)
    for group in (*previous.values(), *current.values()):
        group.sort(key=lambda row: (row.get("line", 0), row.get("end_line", 0)))

    changed = []
    added = []
    removed = []
    ambiguous = []
    unchanged_count = 0
    matched_count = 0
    for key in sorted(set(previous) | set(current)):
        before_rows = previous.get(key, [])
        after_rows = current.get(key, [])
        if not before_rows:
            added.extend(after_rows)
            continue
        if not after_rows:
            removed.extend(before_rows)
            continue

        if len(before_rows) == len(after_rows) == 1:
            pairs = [(before_rows[0], after_rows[0])]
        else:
            by_location = defaultdict(list)
            for row in after_rows:
                by_location[(row.get("line"), row.get("end_line"))].append(row)
            pairs = []
            unmatched_before = []
            for before in before_rows:
                location = (before.get("line"), before.get("end_line"))
                candidates = by_location.get(location, [])
                if len(candidates) == 1:
                    pairs.append((before, candidates.pop()))
                else:
                    unmatched_before.append(before)
            unmatched_after = [row for group in by_location.values() for row in group]
            if unmatched_before or unmatched_after:
                ambiguous.append(
                    {
                        "file": key[0],
                        "qualified_name": key[1],
                        "before": unmatched_before,
                        "after": unmatched_after,
                    }
                )

        for before, after in pairs:
            matched_count += 1
            if before.get("cc_radon") != after.get("cc_radon") or before.get(
                "line"
            ) != after.get("line"):
                changed.append(
                    {
                        "file": key[0],
                        "qualified_name": key[1],
                        "before_line": before.get("line"),
                        "after_line": after.get("line"),
                        "cc_before": before.get("cc_radon"),
                        "cc_after": after.get("cc_radon"),
                        "cc_delta": (
                            None
                            if before.get("cc_radon") is None
                            or after.get("cc_radon") is None
                            else after["cc_radon"] - before["cc_radon"]
                        ),
                    }
                )
            else:
                unchanged_count += 1
    return {
        "baseline_source_sha": sorted(
            {row.get("source_sha") for row in previous_rows if row.get("source_sha")}
        ),
        "current_source_sha": sorted(
            {row.get("source_sha") for row in current_rows if row.get("source_sha")}
        ),
        "baseline_function_count": sum(map(len, previous.values())),
        "current_function_count": sum(map(len, current.values())),
        "matched_function_count": matched_count,
        "unchanged_function_count": unchanged_count,
        "changed_functions": changed,
        "added_functions": [
            {
                field: row.get(field)
                for field in ("file", "qualified_name", "line", "cc_radon")
            }
            for row in added
        ],
        "removed_functions": [
            {
                field: row.get(field)
                for field in ("file", "qualified_name", "line", "cc_radon")
            }
            for row in removed
        ],
        "ambiguous_identity_groups": ambiguous,
    }


def measure_complexity(source_root: str | Path) -> tuple[list[dict], dict]:
    """Measure each named function/method in first-party PyAR production code."""
    root = Path(source_root).resolve()
    paths = _production_files(root)
    radon_version = version("radon")
    file_hashes: dict[str, str] = {}
    file_contents: dict[str, str] = {}
    for path in paths:
        relative = path.relative_to(root).as_posix()
        raw = path.read_bytes()
        file_hashes[relative] = hashlib.sha256(raw).hexdigest()
        file_contents[relative] = raw.decode("utf-8")

    tree_digest = hashlib.sha256()
    for relative, digest in sorted(file_hashes.items()):
        tree_digest.update(relative.encode("utf-8"))
        tree_digest.update(b"\0")
        tree_digest.update(digest.encode("ascii"))
        tree_digest.update(b"\n")
    source_tree_sha256 = tree_digest.hexdigest()
    source_sha = _git_head(root)
    measurement_tool_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    rows = []
    fallback_functions = []
    for relative, source in file_contents.items():
        module = ast.parse(source, filename=relative)
        radon_by_location = {
            (block.name, block.lineno): block.complexity
            for block in _radon_function_blocks(cc_visit(source))
        }
        for node, qualified_name in _named_functions(module):
            key = (node.name, node.lineno)
            measurement_mode = "module"
            complexity = radon_by_location.get(key)
            if complexity is None:
                segment = ast.get_source_segment(source, node)
                isolated_source = textwrap.dedent(segment or "")
                isolated_blocks = [
                    block
                    for block in _radon_function_blocks(cc_visit(isolated_source))
                    if block.name == node.name and block.lineno == 1
                ]
                if not isolated_blocks:
                    raise MeasurementError(
                        f"Radon did not return a function block for {relative}:{node.lineno} {qualified_name}"
                    )
                complexity = isolated_blocks[0].complexity
                measurement_mode = "isolated_ast_segment"
                fallback_functions.append(f"{relative}:{node.lineno}:{qualified_name}")
            rows.append(
                {
                    "file": relative,
                    "qualified_name": qualified_name,
                    "line": node.lineno,
                    "end_line": node.end_lineno,
                    "cc_radon": complexity,
                    "cc_measurement_mode": measurement_mode,
                    "file_sha256": file_hashes[relative],
                    "source_sha": source_sha or "",
                    "source_tree_sha256": source_tree_sha256,
                    "radon_version": radon_version,
                    "tool_version": TOOL_VERSION,
                }
            )
    rows.sort(key=lambda row: (row["file"], row["line"], row["qualified_name"]))
    manifest = {
        "schema_version": 1,
        "measurement_tool": "pyar-complexity-inventory",
        "tool_version": TOOL_VERSION,
        "measurement_tool_sha256": measurement_tool_sha256,
        "radon_version": radon_version,
        "python_version": platform.python_version(),
        "source_sha": source_sha,
        "source_tree_sha256": source_tree_sha256,
        "source_tree_hash_method": "sha256 of sorted relative_path + NUL + file_sha256 + newline",
        "cc_definition": "Radon cc_visit complexity for each named function or method",
        "function_identity": "file, qualified_name, line, end_line",
        "scope": "first-party pyar/**/*.py production functions and methods",
        "exclusions": ["pyar/AIMNet2/**", "__pycache__/**"],
        "function_count": len(rows),
        "isolated_segment_fallback_count": len(fallback_functions),
        "isolated_segment_fallback_functions": fallback_functions,
        "cc_radon": _complexity_summary([row["cc_radon"] for row in rows]),
        "files": file_hashes,
    }
    return rows, manifest


def write_inventory(
    source_root: str | Path,
    output_dir: str | Path,
    compare_baseline_json: str | Path | None = None,
) -> dict:
    """Write deterministic CSV and JSON outputs and return the manifest."""
    rows, manifest = measure_complexity(source_root)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / "function_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if compare_baseline_json is not None:
        baseline_path = Path(compare_baseline_json)
        baseline_rows = json.loads(baseline_path.read_text(encoding="utf-8"))
        if not isinstance(baseline_rows, list):
            raise MeasurementError("Baseline JSON must contain a function-row list")
        reconciliation = reconcile_function_inventories(baseline_rows, rows)
        reconciliation["baseline_artifact"] = baseline_path.name
        reconciliation["baseline_artifact_sha256"] = hashlib.sha256(
            baseline_path.read_bytes()
        ).hexdigest()
        (destination / "baseline_reconciliation.json").write_text(
            json.dumps(reconciliation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return manifest


def main(argv: list[str] | None = None) -> int:
    repository_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=repository_root)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--compare-baseline-json", type=Path)
    arguments = parser.parse_args(argv)
    try:
        manifest = write_inventory(
            arguments.source_root,
            arguments.output_dir,
            compare_baseline_json=arguments.compare_baseline_json,
        )
    except (MeasurementError, OSError, SyntaxError, ValueError) as exc:
        parser.error(str(exc))
    print(
        f"Measured {manifest['function_count']} functions; "
        f"source_sha={manifest['source_sha'] or 'unavailable'}; "
        f"source_tree_sha256={manifest['source_tree_sha256']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
