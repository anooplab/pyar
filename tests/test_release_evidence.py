import copy
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "release_evidence", Path(__file__).parents[1] / "tools/release_evidence.py"
)
evidence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evidence)


@pytest.mark.parametrize(
    "mutation",
    ["none", "sha", "version", "missing", "failed", "stale", "hash", "extra"],
)
def test_candidate_evidence_rejects_invalid_results(tmp_path, mutation):
    for name in ("pyar.whl", "pyar.tar.gz"):
        (tmp_path / name).write_bytes(b"artifact")
    record = {
        "schema_version": 1,
        "sha": "candidate",
        "version": "1.5.1",
        "checks": {"tests": {"sha": "candidate", "conclusion": "success"}},
        "artifacts": {p.name: evidence.digest(p) for p in tmp_path.iterdir()},
    }
    record = copy.deepcopy(record)
    if mutation in ("sha", "version"):
        record[mutation] = "wrong"
    elif mutation == "missing":
        record["checks"] = {}
    elif mutation == "failed":
        record["checks"]["tests"]["conclusion"] = "failure"
    elif mutation == "stale":
        record["checks"]["tests"]["sha"] = "old"
    elif mutation == "hash":
        (tmp_path / "pyar.whl").write_bytes(b"tampered")
    elif mutation == "extra":
        (tmp_path / "old.whl").write_bytes(b"stale")
    if mutation == "none":
        evidence.validate(record, "candidate", "1.5.1", tmp_path, ["tests"])
    else:
        with pytest.raises(ValueError):
            evidence.validate(record, "candidate", "1.5.1", tmp_path, ["tests"])


@pytest.mark.parametrize("outcome", ["xfail", "skip", "empty"])
def test_omission_gate_rejects_unreviewed_outcomes(tmp_path, outcome):
    spec = importlib.util.spec_from_file_location(
        "outcomes", Path(__file__).parents[1] / "tools/check_test_outcomes.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    inventory = tmp_path / "inventory.json"
    inventory.write_text("{}")
    report = tmp_path / "report.xml"
    case = (
        ""
        if outcome == "empty"
        else f'<testcase classname="tests.test_unreviewed.Example" name="case"><skipped type="pytest.{outcome}" message="new omission"/></testcase>'
    )
    report.write_text(f"<testsuite>{case}</testsuite>")
    with pytest.raises(ValueError):
        module.check(report, inventory)
