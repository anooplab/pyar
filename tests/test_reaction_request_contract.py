"""Behavioral evidence for resolver extraction, recorded before implementation."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

import pyar.cli_validation as validation
import pyar.reaction_request as resolver
from reaction_request_support import CASES, snapshot

FIXTURE = json.loads(
    (Path(__file__).parent / "fixtures/reaction_request_contract.json").read_text()
)


@pytest.fixture
def inputs(tmp_path):
    first, second = tmp_path / "first.xyz", tmp_path / "second.xyz"
    first.write_text("2\nfixture\nH 0 0 0\nH .74 0 0\n")
    second.write_text("1\nfixture\nH 0 0 0\n")
    return [str(first), str(second)]


@pytest.mark.parametrize("case", CASES)
def test_complete_resolved_settings_and_restart_payload_match_original(
    inputs, tmp_path, case
):
    options = {"backend": "xtb", "bias_max": 100, **deepcopy(CASES[case])}
    original = deepcopy(options)
    files = list(inputs)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    request = resolver.resolve_reaction_request(inputs, options)
    assert snapshot(request, tmp_path) == FIXTURE["cases"][case]
    assert options == original and inputs == files
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


@pytest.mark.parametrize(
    "options, message",
    [
        ({}, "A backend is required."),
        ({"backend": "xtb"}, "requires --bias-max"),
        (
            {
                "backend": "xtb",
                "bias_max": 100,
                "geometry_optimizer": "native",
                "bias_controller": "invalid",
            },
            "requires --geometry-optimizer geometric",
        ),
        (
            {"backend": "xtb", "bias_max": float("nan"), "bias_controller": "invalid"},
            "Unsupported bias-controller policy",
        ),
        (
            {"backend": "xtb", "bias_max": -1, "orientations": 0},
            "--bias-max must be finite",
        ),
        (
            {"backend": "xtb", "bias_max": 100, "bias_min": -1, "orientations": 0},
            "--bias-min must be finite",
        ),
        (
            {
                "backend": "xtb",
                "bias_max": 100,
                "bias_controller": "fixed",
                "orientations": 0,
            },
            "requires --bias-min",
        ),
        (
            {
                "backend": "xtb",
                "bias_max": 100,
                "orientations": 0,
                "proximity_factor": 0,
            },
            "--orientations must be a positive integer",
        ),
        (
            {
                "backend": "xtb",
                "bias_max": 100,
                "proximity_factor": 0,
                "release_retry_limit": -1,
            },
            "--proximity-factor must be finite",
        ),
        (
            {
                "backend": "xtb",
                "bias_max": 100,
                "release_retry_limit": -1,
                "release_distance_fraction": 2,
            },
            "--release-retry-limit must be nonnegative",
        ),
        (
            {"backend": "xtb", "bias_max": 100, "release_distance_fraction": 2},
            "release tolerances must be positive",
        ),
    ],
)
def test_control_rejection_precedes_input_loading(monkeypatch, options, message):
    monkeypatch.setattr(
        resolver,
        "load_structures",
        lambda *a: pytest.fail("loaded molecules before control validation"),
    )
    with pytest.raises(ValueError, match=message):
        resolver.resolve_reaction_request(["missing-a.xyz", "missing-b.xyz"], options)


@pytest.mark.parametrize(
    "options, message",
    [
        (
            {"site": [99, 0], "softmin_beta": -1},
            "Site index 99 is out of range for reactant A",
        ),
        (
            {"softmin_beta": -1, "bias_potential": "invalid"},
            "Soft-min beta must be a finite positive number",
        ),
        (
            {"bias_potential": "invalid", "nprocs": 0},
            "--bias-potential must be afir or softmin",
        ),
        ({"nprocs": 0, "xtb_model": "invalid"}, "--nprocs must be a positive integer"),
        (
            {"backend": "orca", "xtb_model": "gfn2", "method": ""},
            "--xtb-model requires --backend xtb",
        ),
        ({"method": "BP86", "basis": " "}, "does not support --method or --basis"),
        ({"backend": "orca", "method": "", "basis": ""}, "--method must not be empty"),
        (
            {"backend": "orca", "method": "g-xTB", "opt_threshold": "invalid"},
            "no supported Cartesian provider",
        ),
        (
            {"backend": "orca", "method": "r2scan-3c", "opt_threshold": "invalid"},
            "composite methods are not supported",
        ),
        (
            {
                "backend": "orca",
                "method": "XTB2",
                "basis": "def2-SVP",
                "opt_threshold": "invalid",
            },
            "xTB methods do not accept --basis",
        ),
        (
            {
                "backend": "gaussian",
                "method": "B3LYP",
                "basis": "6-31G",
                "opt_threshold": "invalid",
            },
            "native Gaussian optimization route is not yet qualified",
        ),
    ],
)
def test_loaded_request_validation_preserves_first_error(inputs, options, message):
    with pytest.raises(ValueError, match=message):
        resolver.resolve_reaction_request(
            inputs, {"backend": "xtb", "bias_max": 100, **options}
        )


def test_singlet_scf_rejection_precedes_xtb_model_validation(inputs):
    with pytest.raises(ValueError, match="singlet reactions require --scftype rhf"):
        resolver.resolve_reaction_request(
            inputs,
            dict(
                backend="xtb",
                bias_max=100,
                charge=[0, -1],
                scftype=["uhf"],
                xtb_model="invalid",
            ),
        )


def test_input_shape_precedes_all_options():
    with pytest.raises(ValueError, match="exactly two XYZ reactants"):
        resolver.resolve_reaction_request([], {})


def test_spin_inference_exit_is_not_reclassified(inputs, monkeypatch):
    failure = SystemExit("inference exit")

    def reject(*args):
        raise failure

    monkeypatch.setattr(validation, "_infer_default_multiplicities", reject)
    with pytest.raises(SystemExit) as caught:
        resolver.resolve_reaction_request(inputs, dict(backend="xtb", bias_max=100))
    assert caught.value is failure


def test_spin_validation_exit_retains_translation_and_cause(inputs, monkeypatch):
    failure = SystemExit("validation exit")

    def reject(*args):
        raise failure

    monkeypatch.setattr(validation, "_validate_backend_spin_inputs", reject)
    with pytest.raises(
        ValueError, match="Invalid charge/multiplicity: validation exit"
    ) as caught:
        resolver.resolve_reaction_request(inputs, dict(backend="xtb", bias_max=100))
    assert caught.value.__cause__ is failure


def test_loading_failure_identity_is_preserved(monkeypatch):
    failure = OSError("unreadable geometry")

    def reject(*args):
        raise failure

    monkeypatch.setattr(resolver, "load_structures", reject)
    with pytest.raises(OSError) as caught:
        resolver.resolve_reaction_request(
            ["first.xyz", "second.xyz"], dict(backend="xtb", bias_max=100)
        )
    assert caught.value is failure


def test_resolution_never_preflights_or_launches_a_backend(inputs, monkeypatch):
    monkeypatch.setattr(
        resolver,
        "preflight",
        lambda *a, **k: pytest.fail("preflight during resolution"),
    )
    monkeypatch.setattr(
        "subprocess.Popen", lambda *a, **k: pytest.fail("backend during resolution")
    )
    request = resolver.resolve_reaction_request(
        inputs, dict(backend="xtb", bias_max=100)
    )
    assert request.restart_request == FIXTURE["cases"]["adaptive"]["restart_request"]
