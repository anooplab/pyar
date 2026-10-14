"""Capture audit fixtures from existing tests without changing workflow code.

Run from the repository root with the repository's test environment:
python docs/run_layout_v2/capture_baseline.py --output /tmp/pyar-layout-fixtures
The destination must not exist. Raw snapshots retain their original bytes;
the portable inventory contains only relative paths, sizes and SHA-256 hashes.
Snapshots are audit evidence, not automatically portable restart directories.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[2]
BASELINE = "8ab5bcca43d43fbcf69a76644fad79e1c233b37e"
sys.path.insert(0, str(ROOT))

TESTS = [
    "tests/test_grow.py::test_interrupted_restart_reuses_pool",
    "tests/test_grow.py::test_geometry_end_to_end",
    "tests/test_microsolvation.py::test_restart_resumes_at_next_unfinished_solvent_step",
    "tests/test_aggregate_state.py::AggregateRunStateTests::test_create_and_load_roundtrips_remaining_pathways",
    "tests/test_aggregator.py::AggregatorTests::test_aggregate_resumes_only_uncompleted_persisted_pathways",
    "tests/test_solvation.py::SolvationWorkflowTests::test_solvation_resumes_from_saved_state",
    "tests/test_conformer_workflow.py::ConformerWorkflowTests::test_conformer_search_writes_ranked_rdkit_outputs",
    "tests/test_reaction_state.py::ReactionRunStateTests::test_create_writes_json_and_roundtrips_pending_geometry_metadata",
    "tests/test_reaction_characterization.py::test_characterization_uses_latest_reactant_trace_and_records_route",
    "tests/test_neb.py::test_scan_workflow_continues_using_standardized_neb_handoff",
    "tests/test_neb.py::test_split_endpoint_stages_validate_without_reoptimizing_or_mutating_upstream",
]


class Capture:
    def __init__(self, output):
        self.output = output
        self.nodeid = None
        self.records = []
        self.results = []

    def snapshot(self, source, event):
        if self.nodeid is None or not source.is_dir():
            return
        files = sorted(p for p in source.rglob("*") if p.is_file() and not p.is_symlink())
        if not files:
            return
        number = len(self.records)
        target = self.output / "snapshots" / f"{number:03d}"
        entries = []
        for path in files:
            relative = path.relative_to(source)
            destination = target / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, destination)
            raw = destination.read_bytes()
            entries.append({"path": relative.as_posix(), "size": len(raw),
                            "sha256": hashlib.sha256(raw).hexdigest()})
        self.records.append({"snapshot": f"snapshots/{number:03d}", "test": self.nodeid,
                             "event": event, "files": entries})

    def pytest_runtest_setup(self, item):
        self.nodeid = item.nodeid

    def pytest_runtest_logreport(self, report):
        if report.when == "call" or report.failed:
            self.results.append({"test": report.nodeid, "phase": report.when,
                                 "outcome": report.outcome})

    def pytest_runtest_teardown(self, item):
        temporary = item.funcargs.get("tmp_path")
        if temporary:
            self.snapshot(temporary, "before-test-teardown")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    baseline_check = subprocess.run(
        ["git", "diff", "--quiet", BASELINE, "--", "pyar", "tests", "pyproject.toml"],
        cwd=ROOT, check=False,
    )
    if baseline_check.returncode != 0:
        parser.error("Baseline source/tests/metadata differ or cannot be verified; "
                     "capture from the audited baseline checkout")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    capture = Capture(output)
    original_cleanup = tempfile.TemporaryDirectory.cleanup

    def cleanup(directory):
        capture.snapshot(Path(directory.name), "before-temporary-directory-cleanup")
        original_cleanup(directory)

    patches = [patch.object(tempfile.TemporaryDirectory, "cleanup", cleanup)]
    for module_name, class_name, directory_attribute in [
        ("grow", "GrowRunState", "directory"),
        ("microsolvation", "MicrosolvationRunState", "directory"),
        ("aggregate", "AggregateRunState", "aggregate_directory"),
        ("reaction", "ReactionRunState", "reaction_directory"),
        ("solvation", "SolvationRunState", "solvation_directory"),
    ]:
        owner = getattr(importlib.import_module(f"pyar.state.{module_name}"), class_name)
        original = owner.save

        def save(instance, original=original, attribute=directory_attribute, name=module_name):
            result = original(instance)
            capture.snapshot(Path(getattr(instance, attribute)), f"{name}-state-save")
            return result

        patches.append(patch.object(owner, "save", save))
    for instrument in patches:
        instrument.start()
    try:
        code = pytest.main(["-q", "--basetemp", str(output / "test-work"), *TESTS], plugins=[capture])
    finally:
        for instrument in reversed(patches):
            instrument.stop()
        inventory = {"baseline_commit": BASELINE,
                     "scope": "existing tests; mocked calculators and deterministic geometry-only runs",
                     "results": capture.results, "snapshots": capture.records}
        (output / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
