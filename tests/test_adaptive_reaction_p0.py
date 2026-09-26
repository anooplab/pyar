"""P0 contracts for the adaptive reaction-bias workflow."""

import json
from unittest.mock import Mock, patch

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.workflows import reaction


def _reactants():
    return (
        Molecule(["H"], np.zeros((1, 3)), name="a"),
        Molecule(["H"], np.ones((1, 3)), name="b"),
    )


def test_adaptive_initialization_uses_bias_max_as_single_schedule(tmp_path, monkeypatch):
    reactant_a, reactant_b = _reactants()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        reaction.trial_generation,
        "create_trial_geometries",
        lambda *args: [],
    )
    monkeypatch.setattr(
        reaction,
        "build_gamma_schedule",
        Mock(side_effect=AssertionError("adaptive mode must not build a gamma ladder")),
    )

    _, _, state, gamma_list, orientations, _ = reaction.initialize_reaction_run(
        reactant_a,
        reactant_b,
        None,
        300.0,
        2,
        {"software": "xtb", "bias_controller": "adaptive", "geometry_optimizer": "geometric"},
        None,
        2.3,
    )

    assert list(gamma_list) == [300.0]
    assert orientations == []
    assert state.data["request"]["bias_mode"] == "adaptive"
    assert state.data["gamma_schedule"] == [300.0]


def test_main_cli_dispatches_adaptive_without_bias_min(monkeypatch):
    import pyar.cli as cli

    workflow = Mock()
    cli_result = Mock(status="completed")
    workflow.react.return_value = cli_result
    monkeypatch.setattr(cli, "_log_workflow_result", Mock())
    cli._run_reaction_workflow(
        workflow,
        {"bias_min": None, "bias_max": 300.0, "bias_controller": "adaptive"},
        [object(), object()],
        1,
        {"software": "xtb"},
        None,
    )

    assert workflow.react.call_args.args[2:4] == (None, 300.0)


def test_adaptive_reaction_runs_one_named_cycle(tmp_path, monkeypatch):
    state = Mock()
    state.data = {"products": []}
    monkeypatch.setattr(
        reaction,
        "initialize_reaction_run",
        lambda *args: (str(tmp_path), str(tmp_path), state, [250.0], [], str(tmp_path)),
    )
    optimize_all = Mock(return_value=[])
    monkeypatch.setattr(reaction, "optimize_all", optimize_all)
    monkeypatch.setattr(reaction.os, "chdir", lambda *args: None)

    result = reaction.react(
        object(), object(), None, 250.0, 1,
        {"software": "xtb", "bias_controller": "adaptive"}, None, 2.3,
    )

    assert result.status == "completed_no_candidates"
    assert state.complete_cycle.call_args.args == (250.0, [])
    optimize_all.assert_called_once()
    assert optimize_all.call_args.args[0] == "adaptive"


def test_unbiased_relaxation_forces_native_and_removes_controller_state(monkeypatch, tmp_path):
    molecule = Molecule(["H"], np.zeros((1, 3)), name="candidate")
    captured = {}

    def fake_optimise(candidate, params):
        captured.update(params)
        return True

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(reaction, "optimise", fake_optimise)
    reaction.relax_without_afir_bias(
        molecule,
        {
            "software": "xtb",
            "geometry_optimizer": "geometric",
            "gamma": 250.0,
            "bias_controller": "adaptive",
            "bias_controller_restart": True,
            "bias_alpha_min": 0.01,
            "bias_potential": "softmin",
            "softmin_beta": 2.0,
            "trace_enabled": True,
        },
    )

    assert captured["geometry_optimizer"] == "native"
    assert captured["gamma"] == 0.0
    assert captured["trace_enabled"] is False
    assert captured["reaction_trace"] is False
    assert captured["opt_threshold"] == reaction.defualt_parameters.values["opt_threshold"]
    assert captured["opt_cycles"] == reaction.defualt_parameters.values["opt_cycles"]
    assert not any(key.startswith("bias_") for key in captured)
    assert "softmin_beta" not in captured


def test_release_retry_uses_a_fresh_unbiased_relaxation_job(monkeypatch, tmp_path):
    molecule = Molecule(["H"], np.zeros((1, 3)), name="candidate")
    job_names = []
    monkeypatch.chdir(tmp_path)

    def fake_optimise(candidate, params):
        job_names.append(candidate.name)
        return True

    monkeypatch.setattr(reaction, "optimise", fake_optimise)
    reaction.relax_without_afir_bias(
        molecule, {"software": "xtb", "release_retry_attempt": 0}
    )
    reaction.relax_without_afir_bias(
        molecule, {"software": "xtb", "release_retry_attempt": 1}
    )

    assert job_names == ["relax", "relax_attempt_1"]


def test_release_probe_rejects_bond_heuristic_contact_outside_radius_sum(tmp_path, monkeypatch):
    from pyar.data import new_atomic_data

    class Geometry:
        name = "candidate"
        atoms_list = ["C", "N"]

        def __init__(self, distance):
            self.coordinates = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])

        def mol_to_xyz(self, path):
            with open(path, "w") as stream:
                stream.write("candidate geometry\n")

    radii_sum = (
        new_atomic_data.covalent_radius["C"]
        + new_atomic_data.covalent_radius["N"]
    )
    geometry = Geometry(radii_sum * 1.05)
    monkeypatch.chdir(tmp_path)

    survived = reaction._write_release_probe(
        geometry, geometry, {}, True,
        {"forming_pairs": [[0, 1]]},
    )

    assert survived is False


def test_release_probe_handles_failed_relaxation_without_coordinates(tmp_path, monkeypatch):
    class Geometry:
        name = "candidate"
        atoms_list = ["C", "N"]
        coordinates = None

        def mol_to_xyz(self, path):
            with open(path, "w") as stream:
                stream.write("candidate geometry\n")

    monkeypatch.chdir(tmp_path)
    survived = reaction._write_release_probe(
        Geometry(), Geometry(), {}, False, {"forming_pairs": [[0, 1]]}
    )

    assert survived is False
    probe = json.loads((tmp_path / "release_probe.json").read_text())
    assert probe["post_relax_distance"] == []
    assert probe["release_state"] == "RELEASE_FAILED"


def test_unbiased_relaxation_rejects_backend_without_native_optimization(monkeypatch, tmp_path):
    molecule = Molecule(["H"], np.zeros((1, 3)), name="candidate")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(reaction, "backend_supports_native_optimization", lambda _: False)

    with pytest.raises(ValueError, match="does not advertise native optimization"):
        reaction.relax_without_afir_bias(molecule, {"software": "unsupported"})


def test_release_retry_escalates_margin_once_and_preserves_ceiling(monkeypatch):
    class Candidate:
        coordinates = np.zeros((1, 3))
        name = "candidate"

        def copy(self):
            return Candidate()

        def is_bonded(self):
            return True

    evidence = {
        "state": "CANDIDATE",
        "forming_pairs": [[0, 1]],
        # High instantaneous resistance no longer disqualifies a failed
        # candidate from continuing under bias.
        "alpha_critical": 0.2,
        "persistence_counter": 3,
        "bond_order_available": False,
    }
    retry_params = {
        "software": "xtb",
        "bias_controller": "adaptive",
        "bias_alpha_margin": 0.002,
        "gamma": 300.0,
    }
    prepare = Mock(return_value=retry_params)
    monkeypatch.setattr(reaction, "_prepare_release_retry", prepare)
    monkeypatch.setattr(reaction, "_set_release_outcome", Mock())
    monkeypatch.setattr(reaction, "relax_without_afir_bias", Mock(return_value=True))
    monkeypatch.setattr(reaction, "_write_release_probe", Mock(side_effect=[False, True]))
    monkeypatch.setattr(reaction, "optimise", Mock(return_value=True))
    monkeypatch.setattr(reaction, "_release_evidence", Mock(return_value=evidence))

    result = reaction._adaptive_release_probe(
        Candidate(), Candidate(),
        {**retry_params, "bias_alpha_margin": 0.001, "bias_max": 300.0}, evidence,
    )

    assert result[4] is True
    assert result[5] == 1
    prepare.assert_called_once()
    assert prepare.call_args.args[2] == pytest.approx(0.002)


def test_retry_that_is_not_candidate_updates_summary_and_preserves_probe_history(
    monkeypatch, tmp_path
):
    class Candidate:
        coordinates = np.zeros((1, 3))
        name = "candidate"

        def copy(self):
            return Candidate()

        def is_bonded(self):
            return True

        def mol_to_xyz(self, path):
            (tmp_path / path).write_text("retry geometry\n")

    evidence = {
        "state": "CANDIDATE",
        "forming_pairs": [[0, 1]],
        "persistence_counter": 3,
    }
    retry_evidence = {"state": "DRIVING", "reason": "contact_not_persistent"}
    retry_params = {"release_retry_attempt": 1, "bias_alpha_margin": 0.002}
    monkeypatch.chdir(tmp_path)
    (tmp_path / "release_probe.json").write_text(json.dumps({
        "release_attempt": 0,
        "release_state": "RELEASE_FAILED",
    }))
    first_probe = json.dumps({"release_attempt": 0, "release_state": "RELEASE_FAILED"})
    (tmp_path / "release_attempts.jsonl").write_text(first_probe + "\n")
    monkeypatch.setattr(reaction, "_set_release_outcome", Mock())
    monkeypatch.setattr(reaction, "relax_without_afir_bias", Mock(return_value=True))
    monkeypatch.setattr(reaction, "_write_release_probe", Mock(return_value=False))
    monkeypatch.setattr(reaction, "_prepare_release_retry", Mock(return_value=retry_params))
    monkeypatch.setattr(reaction, "optimise", Mock(return_value=True))
    monkeypatch.setattr(reaction, "_release_evidence", Mock(return_value=retry_evidence))

    result = reaction._adaptive_release_probe(
        Candidate(), Candidate(),
        {"release_retry_limit": 1, "release_margin_factor": 2.0}, evidence,
    )

    summary = json.loads((tmp_path / "release_probe.json").read_text())
    history = (tmp_path / "release_attempts.jsonl").read_text().splitlines()
    assert result[4] is False
    assert summary["release_attempt"] == 0
    assert summary["retry_outcome"] == "retry_not_candidate"
    assert summary["retry_attempt"] == 1
    assert summary["retry_reason"] == "release_evidence_not_candidate"
    assert summary["retry_evidence"] == retry_evidence
    assert (tmp_path / summary["retry_geometry"]).exists()
    assert history == [first_probe]


def test_release_retry_resets_contact_persistence_but_keeps_controller_load(tmp_path, monkeypatch):
    class Candidate:
        name = "candidate"
        coordinates = np.asarray([[0.0, 0.0, 0.0]])

    monkeypatch.chdir(tmp_path)
    job_dir = tmp_path / "job_candidate"
    job_dir.mkdir()
    checkpoint = {
        "configuration": {"qc_params": {"bias_alpha_margin": 0.001}},
        "controller": {"configuration": {"safety_margin": 0.001}, "decision": {"alpha": 0.2}},
        "release_tracker": {"state": "CANDIDATE"},
        "positions_angstrom": [[1.0, 0.0, 0.0]],
    }
    path = job_dir / "pyar_bias_controller_state.json"
    path.write_text(json.dumps(checkpoint))

    params = reaction._prepare_release_retry(Candidate(), {"bias_alpha_margin": 0.001}, 0.002, 1)
    saved = json.loads(path.read_text())
    assert params["bias_controller_restart"] is True
    assert params["adaptive_max_segments"] == 8
    assert params["adaptive_suppress_release_candidate"] is True
    assert saved["controller"]["decision"]["alpha"] == pytest.approx(0.2)
    assert saved["controller"]["configuration"]["safety_margin"] == pytest.approx(0.002)
    assert "release_tracker" not in saved
