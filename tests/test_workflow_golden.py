"""Golden artifacts recorded from the exact pre-decomposition source commit."""

import json
from pathlib import Path

import pytest

from workflow_golden_support import run_scenario


@pytest.mark.parametrize(
    "kind",
    [
        "conformer",
        "aggregate",
        "grow",
        "microsolvation",
        "microsolvation_resampling",
        "reaction",
        "scan_inter",
        "scan_intra",
    ],
)
def test_workflow_matches_original_complete_evidence(
    kind, tmp_path, monkeypatch, deterministic_workflow_backend
):
    actual = run_scenario(
        kind, tmp_path / "run", monkeypatch, deterministic_workflow_backend
    )
    expected = json.loads(
        (Path(__file__).parent / "fixtures/golden" / kind / "evidence.json").read_text()
    )
    assert actual == expected
