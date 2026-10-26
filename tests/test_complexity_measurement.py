import ast
import csv
import json
import subprocess
import textwrap
from pathlib import Path

import pytest
from radon.complexity import cc_visit

from tools.measure_complexity import (
    MeasurementError,
    measure_complexity,
    reconcile_function_inventories,
    write_inventory,
)


def _fixture_source(root):
    package = root / "pyar"
    package.mkdir(parents=True)
    (package / "sample.py").write_text(
        """def linear(value):
    return value

def decision(value):
    if value and value > 0:
        return True
    return False

async def asynchronous(value):
    return value

class Worker:
    def method(self, value):
        if value:
            return value
        return None

def outer():
    def inner(value):
        if value:
            return value
        return None
    class Local:
        def method(self, value):
            if value:
                return value
            return None
    return inner(True)
""",
        encoding="utf-8",
    )
    vendored = package / "AIMNet2"
    vendored.mkdir()
    (vendored / "vendor.py").write_text(
        "def excluded():\n    if True:\n        return 1\n", encoding="utf-8"
    )
    cache = package / "__pycache__"
    cache.mkdir()
    (cache / "generated.py").write_text(
        "def generated():\n    return 1\n", encoding="utf-8"
    )


def test_inventory_matches_independent_function_and_decision_counts(tmp_path):
    _fixture_source(tmp_path)

    rows, manifest = measure_complexity(tmp_path)

    by_name = {row["qualified_name"]: row["cc_radon"] for row in rows}
    assert by_name == {
        "linear": 1,
        "decision": 3,
        "asynchronous": 1,
        "Worker.method": 2,
        "outer": 1,
        "outer.inner": 2,
        "outer.Local.method": 2,
    }
    independent_count = sum(
        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        for node in ast.walk(ast.parse((tmp_path / "pyar/sample.py").read_text()))
    )
    assert len(rows) == independent_count == manifest["function_count"]
    assert all("AIMNet2" not in row["file"] for row in rows)
    assert all(
        row["source_tree_sha256"] == manifest["source_tree_sha256"] for row in rows
    )
    assert len(manifest["source_tree_sha256"]) == 64
    assert len(manifest["measurement_tool_sha256"]) == 64
    assert manifest["isolated_segment_fallback_count"] == 1


def test_inventory_outputs_are_deterministic_and_source_bound(tmp_path):
    _fixture_source(tmp_path)
    first = tmp_path / "first"
    second = tmp_path / "second"
    baseline = tmp_path / "baseline.json"
    baseline.write_text(
        json.dumps(
            [
                {
                    "file": "pyar/sample.py",
                    "qualified_name": "linear",
                    "line": 1,
                    "cc_radon": 1,
                    "primary": True,
                    "source_sha": "fixture-baseline",
                }
            ]
        ),
        encoding="utf-8",
    )

    first_manifest = write_inventory(tmp_path, first, compare_baseline_json=baseline)
    second_manifest = write_inventory(tmp_path, second, compare_baseline_json=baseline)

    assert first_manifest == second_manifest
    assert (first / "function_metrics.csv").read_bytes() == (
        second / "function_metrics.csv"
    ).read_bytes()
    assert (first / "manifest.json").read_bytes() == (
        second / "manifest.json"
    ).read_bytes()
    with (first / "function_metrics.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == first_manifest["function_count"]
    assert json.loads((first / "manifest.json").read_text()) == first_manifest
    first_reconciliation = (first / "baseline_reconciliation.json").read_bytes()
    assert (
        first_reconciliation == (second / "baseline_reconciliation.json").read_bytes()
    )
    assert json.loads(first_reconciliation)["baseline_source_sha"] == [
        "fixture-baseline"
    ]

    before_digest = first_manifest["source_tree_sha256"]
    (tmp_path / "pyar/sample.py").write_text(
        "def changed():\n    return 1\n", encoding="utf-8"
    )
    _, changed_manifest = measure_complexity(tmp_path)
    assert changed_manifest["source_tree_sha256"] != before_digest


def test_repository_inventory_binds_to_git_head_and_matches_ast_count():
    repository = Path(__file__).resolve().parents[1]
    rows, manifest = measure_complexity(repository)
    expected_sha = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    expected_count = 0
    for path in (repository / "pyar").rglob("*.py"):
        relative = path.relative_to(repository)
        if "AIMNet2" in relative.parts or "__pycache__" in relative.parts:
            continue
        expected_count += sum(
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
        )

    assert manifest["source_sha"] == expected_sha
    assert manifest["function_count"] == len(rows) == expected_count
    assert {row["source_sha"] for row in rows} == {expected_sha}
    assert {row["source_tree_sha256"] for row in rows} == {
        manifest["source_tree_sha256"]
    }

    source_cache = {}
    node_cache = {}
    for row in rows:
        relative = row["file"]
        if relative not in source_cache:
            source_cache[relative] = (repository / relative).read_text(encoding="utf-8")
        source = source_cache[relative]
        if relative not in node_cache:
            tree = ast.parse(source, filename=relative)
            node_cache[relative] = {
                (node.name, node.lineno): node
                for node in ast.walk(tree)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            }
        node = node_cache[relative][(row["qualified_name"].split(".")[-1], row["line"])]
        isolated_source = textwrap.dedent(ast.get_source_segment(source, node))
        independent_blocks = [
            block
            for block in cc_visit(isolated_source)
            if block.name == node.name and block.lineno == 1
        ]
        assert len(independent_blocks) == 1, (relative, row["qualified_name"])
        assert row["cc_radon"] == independent_blocks[0].complexity, (
            relative,
            row["qualified_name"],
        )


def test_inventory_requires_a_pyar_source_root(tmp_path):
    with pytest.raises(MeasurementError, match="no pyar package"):
        measure_complexity(tmp_path)


def test_reconciliation_reports_changed_added_and_removed_definitions():
    previous = [
        {
            "file": "pyar/a.py",
            "qualified_name": "kept",
            "line": 1,
            "cc_radon": 4,
            "primary": True,
            "source_sha": "old",
        },
        {
            "file": "pyar/a.py",
            "qualified_name": "removed",
            "line": 5,
            "cc_radon": 2,
            "primary": True,
            "source_sha": "old",
        },
        {
            "file": "tests/test_a.py",
            "qualified_name": "test_x",
            "line": 1,
            "cc_radon": 1,
            "primary": False,
            "source_sha": "old",
        },
    ]
    current = [
        {
            "file": "pyar/a.py",
            "qualified_name": "kept",
            "line": 3,
            "cc_radon": 5,
            "source_sha": "new",
        },
        {
            "file": "pyar/a.py",
            "qualified_name": "added",
            "line": 8,
            "cc_radon": 3,
            "source_sha": "new",
        },
    ]

    result = reconcile_function_inventories(previous, current)

    assert result["baseline_source_sha"] == ["old"]
    assert result["current_source_sha"] == ["new"]
    assert result["matched_function_count"] == 1
    assert result["changed_functions"] == [
        {
            "file": "pyar/a.py",
            "qualified_name": "kept",
            "before_line": 1,
            "after_line": 3,
            "cc_before": 4,
            "cc_after": 5,
            "cc_delta": 1,
        }
    ]
    assert [
        (row["qualified_name"], row["cc_radon"]) for row in result["added_functions"]
    ] == [("added", 3)]
    assert [row["qualified_name"] for row in result["removed_functions"]] == ["removed"]


def test_reconciliation_does_not_collapse_property_accessors():
    accessors_before = [
        {
            "file": "pyar/model.py",
            "qualified_name": "Model.value",
            "line": line,
            "end_line": line + 1,
            "cc_radon": 1,
            "primary": True,
        }
        for line in (10, 20)
    ]
    accessors_after = [dict(row, source_sha="new") for row in accessors_before]

    result = reconcile_function_inventories(accessors_before, accessors_after)

    assert result["baseline_function_count"] == 2
    assert result["current_function_count"] == 2
    assert result["matched_function_count"] == 2
    assert result["unchanged_function_count"] == 2
    assert result["ambiguous_identity_groups"] == []
