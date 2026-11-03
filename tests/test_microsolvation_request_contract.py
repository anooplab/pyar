"""Preserve constructor ordering, normalization, ownership and restart fingerprints."""

from collections.abc import Mapping
from dataclasses import FrozenInstanceError
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import pyar.backend_capabilities as capabilities
import pyar.molecule_merge as merge
from pyar.microsolvation.request import (
    MicrosolvationRequest,
    MicrosolvationRequestError,
)
from microsolvation_request_support import CASES, fingerprint, request_values, water

ORIGINAL = json.loads(
    (
        Path(__file__).parent / "fixtures/microsolvation_request_contract.json"
    ).read_text()
)


@pytest.mark.parametrize("case", CASES)
def test_complete_serialization_fingerprints_match_original(case):
    request = MicrosolvationRequest(**request_values(case))
    assert fingerprint(request) == ORIGINAL["cases"][case]
    if case == "sites":
        assert request.site == (1, 2)
        assert request.site_groups == (("hydrogens", (1, 2)),)


@pytest.mark.parametrize(
    "options,message",
    [
        ({"count": 0, "objective": "bad"}, "count must be a positive integer"),
        (
            {"number_of_orientations": True, "maximum_number_of_seeds": 0},
            "orientations must be a positive integer",
        ),
        (
            {"maximum_number_of_seeds": False, "surface_points_per_atom": 0},
            "maximum_number_of_seeds must be a positive integer",
        ),
        (
            {"probe_radius": -1, "shell_tolerance": 0},
            "probe radius must be finite and non-negative",
        ),
        ({"objective": "bad", "wall": "bad"}, "objective must be 'shell' or 'cluster'"),
        (
            {"wall": "bad", "placement_score": "bad"},
            "wall must be 'none', 'outer-sphere', or 'auto'",
        ),
        (
            {"placement_score": "bad", "resampling_moves": -1},
            "placement score must be 'geometric' or 'energy'",
        ),
        (
            {"placement_score": "energy", "site": []},
            "energy placement scoring requires --backend",
        ),
        (
            {"objective": "cluster", "resampling_moves": 1},
            "cluster objective requires --placement-score energy",
        ),
        (
            {"resampling_moves": -1, "site": []},
            "resampling moves must be a non-negative integer",
        ),
        (
            {"resampling_moves": 1, "site": []},
            "whole-cluster resampling requires --objective cluster",
        ),
        (
            {
                "objective": "cluster",
                "placement_score": "energy",
                "resampling_moves": 1,
                "backend_parameters": {"software": "xtb"},
            },
            "energy-gradient backend and geometric optimizer",
        ),
        (
            {"resampling_evaluation_budget": 1, "site": []},
            "resampling budgets require --resampling-moves",
        ),
        (
            {"site": [1, 1], "site_groups": [("bad", [])]},
            "--site requires one or more distinct",
        ),
        (
            {"site": [3], "site_groups": [("bad", [])]},
            "Solute site index 3 is out of range",
        ),
        (
            {"site_groups": [("same", [0]), ("same", [1])], "contact_distance": -1},
            "site groups require unique names",
        ),
    ],
)
def test_constructor_preserves_first_error(options, message):
    values = request_values("default")
    values.update(options)
    with pytest.raises(MicrosolvationRequestError, match=message):
        MicrosolvationRequest(**values)


def test_each_molecule_finishes_validation_before_the_next_molecule():
    values = request_values("default")
    values["solute"].multiplicity = 2
    values["solvent"].coordinates[0, 0] = np.nan
    with pytest.raises(
        MicrosolvationRequestError, match="Invalid charge/multiplicity for solute"
    ):
        MicrosolvationRequest(**values)


def test_geometry_validation_precedes_unresolved_state():
    values = request_values("default")
    values["solute"].coordinates[0, 0] = np.nan
    values["solute"].charge = None
    with pytest.raises(
        MicrosolvationRequestError,
        match="solute must contain atoms with finite coordinates",
    ):
        MicrosolvationRequest(**values)


class UnreadableBackend(Mapping):
    def __init__(self, failure):
        self.failure = failure
        self.read = False

    def __len__(self):
        return 1

    def __getitem__(self, key):
        raise AssertionError("unexpected key access")

    def __iter__(self):
        self.read = True
        raise self.failure


@pytest.mark.parametrize("probe", [-1, 1.4])
def test_backend_copy_retains_position_between_surface_and_objective_validation(probe):
    failure = RuntimeError("cannot copy backend mapping")
    backend = UnreadableBackend(failure)
    values = request_values("default")
    values.update(probe_radius=probe, objective="bad", backend_parameters=backend)
    if probe < 0:
        with pytest.raises(MicrosolvationRequestError, match="probe radius"):
            MicrosolvationRequest(**values)
        assert not backend.read
    else:
        with pytest.raises(RuntimeError) as caught:
            MicrosolvationRequest(**values)
        assert caught.value is failure and backend.read


def test_frozen_request_retains_outer_copy_and_nested_backend_ownership():
    nested = {"source": "caller"}
    backend = {"metadata": nested}
    request = MicrosolvationRequest(
        water("solute"), water("solvent"), 1, backend_parameters=backend
    )
    assert request.backend_parameters is not backend
    assert request.backend_parameters["metadata"] is nested
    backend["new"] = "caller-only"
    assert "new" not in request.backend_parameters
    nested["source"] = "shared"
    assert request.backend_parameters["metadata"]["source"] == "shared"
    with pytest.raises(FrozenInstanceError):
        request.count = 3


def test_capability_rejection_identifies_first_unsupported_grown_state(monkeypatch):
    monkeypatch.setattr(
        capabilities,
        "get_backend_capabilities",
        lambda name: SimpleNamespace(supports_charge=False, supports_multiplicity=True),
    )
    values = request_values("default")
    values["backend_parameters"] = {"software": "xtb"}
    values["solvent"].charge = 1
    values["solvent"].multiplicity = 2
    with pytest.raises(
        MicrosolvationRequestError, match="combined electronic state after 1 solvent"
    ):
        MicrosolvationRequest(**values)


def test_spin_combination_failure_precedes_later_capability_decisions(monkeypatch):
    combinations = []

    def wrong_spin(first, second):
        combinations.append((first, second))
        return 2

    monkeypatch.setattr(merge, "combine_multiplicity", wrong_spin)
    with pytest.raises(
        MicrosolvationRequestError, match="invalid cluster state after 1 solvent"
    ):
        MicrosolvationRequest(**request_values("default"))
    assert combinations == [(1, 1)]


def test_spin_combination_is_not_advanced_past_requested_count(monkeypatch):
    calls = []
    original = merge.combine_multiplicity

    def combine(first, second):
        calls.append((first, second))
        return original(first, second)

    monkeypatch.setattr(merge, "combine_multiplicity", combine)
    MicrosolvationRequest(**request_values("default"))
    assert calls == [(1, 1), (1, 1)]


def test_request_validation_and_serialization_do_not_launch_calculations(monkeypatch):
    monkeypatch.setattr(
        "subprocess.Popen",
        lambda *a, **k: pytest.fail("backend during request validation"),
    )
    assert (
        MicrosolvationRequest(**request_values("automatic_shell")).to_state_dict()[
            "confinement"
        ]["mode"]
        == "case-dependent-outer-sphere-v1"
    )
