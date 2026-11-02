"""Behavioral evidence for responsibility extraction, checked against its baseline.

These assertions protect software evidence and reporting, not real-backend science.
"""

import hashlib
import json
import subprocess
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pyar import optimization_service
from pyar.recovery import controller
from pyar.recovery.models import (
    ExecutionObservation,
    FailureCode,
    FailureDiagnostic,
    GeometryFrame,
    Outcome,
    ProcessRecord,
)
from pyar.recovery.parsers import NativeObservation
from pyar.recovery.policy import RecoveryPolicy


@pytest.mark.parametrize(
    "changes, code, trusted, diagnostic, exception, expected",
    [
        ({}, 0, True, None, None, Outcome.CONVERGED),
        ({}, 0, False, None, None, Outcome.INVALID_RESULT),
        ({}, None, True, None, None, Outcome.INVALID_RESULT),
        ({}, 7, True, None, None, Outcome.FAILED_FINAL),
        ({"stationary": False}, 0, True, None, None, Outcome.INVALID_RESULT),
        ({"terminated": False}, 0, True, None, None, Outcome.INVALID_RESULT),
        ({"accepted": False}, 0, True, None, None, Outcome.INVALID_RESULT),
        ({"frame": None}, 0, True, None, None, Outcome.INVALID_RESULT),
        (
            {"converged": False, "limited": True},
            1,
            True,
            None,
            None,
            Outcome.INCOMPLETE_USABLE,
        ),
        (
            {"converged": False, "limited": True, "accepted": False},
            1,
            True,
            None,
            None,
            Outcome.FAILED_RECOVERABLE,
        ),
        (
            {"converged": False, "limited": True},
            1,
            True,
            None,
            RuntimeError("adapter"),
            Outcome.FAILED_RECOVERABLE,
        ),
        ({}, 0, True, FailureCode.CALC_SCF_MAX_STEPS, None, Outcome.FAILED_FINAL),
        (
            {"limited": True},
            0,
            False,
            FailureCode.CALC_SCF_MAX_STEPS,
            None,
            Outcome.INVALID_RESULT,
        ),
        ({}, 0, True, FailureCode.USER_INPUT_INVALID, None, Outcome.FAILED_FINAL),
    ],
)
def test_observation_policy_preserves_evidence_precedence(
    changes, code, trusted, diagnostic, exception, expected
):
    frame = GeometryFrame(("H",), ((0.0, 0.0, 0.0),), -0.5, -0.5)
    observed = NativeObservation("6.7.1", True, True, False, True, frame, True)
    observed = replace(observed, **changes)
    if diagnostic is not None:
        observed = replace(
            observed, diagnostics=(FailureDiagnostic(diagnostic, "backend", "fixture"),)
        )
    request = SimpleNamespace(optimizer="native", backend="xtb")
    processes = () if code is None else (ProcessRecord("xtb", code),)
    execution = ExecutionObservation("0" * 64, (), processes, True, True)
    result = optimization_service._classify_observation(
        request, execution, exception, (), trusted, observed
    )
    assert result.status is expected
    assert result.exception is exception
    assert result.accepted_frame == (
        observed.frame if observed.accepted and trusted else None
    )
    assert result.observed_frame == observed.frame
    assert [gate.axis for gate in result.gates] == [
        "backend",
        "scf",
        "optimizer",
        "stationarity",
        "objective",
        "constraints",
        "workflow",
    ]
    assert result.gates[-1].state.value == "unknown"


@pytest.mark.parametrize("failure", ["nonzero", "no_receipt", "timeout"])
def test_unreceipted_attempt_charges_reserved_budget_and_never_replays(
    tmp_path, monkeypatch, failure
):
    molecule = SimpleNamespace(
        name="h2",
        atoms_list=["H", "H"],
        coordinates=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        charge=0,
        multiplicity=1,
        scftype="rhf",
        fragments=[],
    )
    settings = dict(
        software="xtb",
        geometry_optimizer="geometric",
        xtb_model="gfn2",
        opt_cycles=3,
        opt_threshold="normal",
        gamma=None,
    )
    policy = RecoveryPolicy("safe", max_evaluations=10)
    calls = []

    def execute(directory, timeout):
        calls.append(directory)
        (directory / "worker.out").write_text("interrupted worker\n")
        if failure == "timeout":
            raise subprocess.TimeoutExpired("worker", timeout)
        return 9 if failure == "nonzero" else 0

    monkeypatch.setattr(controller, "_execute_attempt", execute)
    first = controller.run_recovery(molecule, settings, policy, directory=tmp_path)
    assert first.stop_reason == (
        "wall_budget_exhausted"
        if failure == "timeout"
        else "worker_failed_without_receipt"
    )
    assert first.reserved_cycles == 3
    assert first.backend_evaluations == 10
    assert first.result.attempts[0].worker_returncode == (
        None if failure == "timeout" else 9 if failure == "nonzero" else 0
    )
    replay = controller.run_recovery(
        molecule, settings, policy, directory=tmp_path, resume=True
    )
    assert replay.stop_reason == "interrupted_attempt_without_receipt"
    assert replay.reserved_cycles == 3 and replay.backend_evaluations == 10
    assert len(calls) == 1
    assert not (tmp_path / "job_h2/result_h2.xyz").exists()


def test_conformer_meta_analysis_preserves_noise_pairs_and_macro_weights(
    tmp_path, monkeypatch
):
    from pyar.scripts import meta_analyze_clustering_benchmarks as analysis

    data = tmp_path / "mpconf196gen"
    data.mkdir()
    source = data / "fixture_crest_conformers.xyz"
    source.write_bytes(b"controlled ensemble boundary")
    run = tmp_path / "fixture"
    (run / "results").mkdir(parents=True)
    (run / "manifest.json").write_text(
        json.dumps({"ensemble_sha256": hashlib.sha256(source.read_bytes()).hexdigest()})
    )
    (run / "reference_basins.json").write_text(
        json.dumps(
            {
                "structures": [
                    {"name": f"frame_{i:04d}", "basin_id": basin}
                    for i, basin in enumerate([0, 0, 1])
                ]
            }
        )
    )
    conditions = [
        {
            "algorithm_requested": "hybrid",
            "feature_requested": "fixture",
            "basins_retained": 2,
            "diagnostics": {
                "labels": labels,
                "algorithm_used": "auto",
                "feature_used": "fixture",
            },
        }
        for labels in ([0, 0, -1], [0, 1, 1])
    ]
    (run / "results/comparison.json").write_text(
        json.dumps(
            {"seed_budget": 2, "reference_basin_count": 2, "conditions": conditions}
        )
    )
    monkeypatch.setattr(analysis, "CONFORMERS", ("fixture",))
    monkeypatch.setattr(analysis, "FEATURES", ("fixture",))
    monkeypatch.setattr(analysis, "ALGORITHMS", ("auto",))
    monkeypatch.setattr(
        analysis,
        "load_ensemble",
        lambda _: [{"atoms": object(), "name": f"frame_{i:04d}"} for i in range(3)],
    )
    monkeypatch.setattr(
        analysis, "_as_molecule", lambda _atoms, name: SimpleNamespace(name=name)
    )
    monkeypatch.setattr(
        analysis, "_feature_roc_auc", lambda *_: (0.75, "fixture", [], None)
    )
    confusion, macro, auc, auc_macro = analysis.analyze_conformers(tmp_path, tmp_path)
    assert [(r["tp"], r["fp"], r["fn"], r["tn"]) for r in confusion] == [
        (1, 0, 0, 2),
        (0, 1, 1, 1),
    ]
    assert [r["algorithm_requested"] for r in confusion] == ["auto", "auto"]
    assert macro[0]["macro_sensitivity"] == 0.5
    assert macro[0]["macro_specificity"] == 0.75
    assert macro[0]["macro_seed_budget_basin_efficiency"] == 1.0
    assert (
        auc[0]["positive_same_basin_pairs"] == 1
        and auc[0]["negative_different_basin_pairs"] == 2
    )
    assert auc_macro[0]["macro_pairwise_same_basin_auc"] == 0.75


@pytest.mark.parametrize("length", [0, 1, 3, 5])
def test_wall_force_matches_finite_differences_and_zero_total_force(length):
    from pyar.microsolvation.wall import outer_sphere_energy_gradient

    # Analytical, dimensionless fixture; finite-difference step and tolerances
    # only verify this expression, not a calibrated chemical wall model.
    xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [float(length), 1.0, 0.0]])
    masses = np.array([1.0, 2.0, 3.0])
    options = dict(radius_angstrom=2.0, force_constant_hartree_per_angstrom2=0.01)
    _, gradient = outer_sphere_energy_gradient(xyz, masses, [0, 1], **options)
    np.testing.assert_allclose(gradient.sum(axis=0), 0.0, atol=1e-14)
    for atom in range(3):
        for axis in range(3):
            plus, minus = xyz.copy(), xyz.copy()
            plus[atom, axis] += 1e-5
            minus[atom, axis] -= 1e-5
            upper = outer_sphere_energy_gradient(plus, masses, [0, 1], **options)[0]
            lower = outer_sphere_energy_gradient(minus, masses, [0, 1], **options)[0]
            assert gradient[atom, axis] == pytest.approx(
                (upper - lower) / 2e-5, abs=1e-11
            )
