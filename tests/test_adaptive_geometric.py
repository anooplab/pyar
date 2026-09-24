"""Regression coverage for accepted-segment adaptive bias integration."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
from ase import Atoms
from ase.units import Bohr, Hartree

from pyar.backends.geometric import (
    Geometric, PyarGeometricCalculator, _CONTROLLER_STATE_FILE,
)
from pyar.energy_gradient_providers import EnergyGradientResult
from pyar.reaction_analysis import analyse_reaction_trace
from pyar.reaction_trace import load_trace_records
from pyar.workflows.reaction import without_afir_bias


class HarmonicProvider:
    def evaluate(self, atoms, coordinates_bohr):
        vector = coordinates_bohr[1] - coordinates_bohr[0]
        distance = np.linalg.norm(vector)
        derivative = 0.1 * (distance - 4.0)
        gradient = derivative * vector / distance
        return EnergyGradientResult(0.05 * (distance - 4.0)**2, np.array([-gradient, gradient]))


def pair(distance):
    return Atoms("CC", positions=[[0., 0., 0.], [distance * 0.52917726, 0., 0.]])


@pytest.fixture
def calculation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("pyar.backends.geometric._resolve_backend_evaluator",
                        lambda *args: HarmonicProvider())
    return dict(software="xtb", gamma=100., bias_controller="adaptive",
                bias_alpha_margin=0.001, bias_alpha_smoothing=0.5, trace_enabled=True)


@pytest.mark.parametrize("kind", ["afir", "softmin"])
def test_frozen_segment_force_matches_energy_difference(calculation, kind):
    calculator = PyarGeometricCalculator(dict(calculation, bias_potential=kind), [[0], [1]])
    calculator.calculate(pair(3.8))
    state = calculator.bias_controller.state_dict()
    force = calculator.results["forces"][1, 0] * Bohr / Hartree
    energies = []
    h = 1e-5
    for distance in (3.8+h, 3.8-h, 3.8+h):
        calculator.calculate(pair(distance))
        energies.append(calculator.results["energy"] / Hartree)
        assert calculator.bias_controller.state_dict() == state
    assert energies[0] == energies[2]
    assert force == pytest.approx(-(energies[0]-energies[1])/(2*h), rel=1e-6)


@pytest.mark.parametrize("kind", ["afir", "softmin"])
def test_acceptance_preserves_energy_and_restores_history(calculation, kind):
    params = dict(calculation, bias_potential=kind)
    calculator = PyarGeometricCalculator(params, [[0], [1]])
    calculator.calculate(pair(3.8))
    calculator.calculate(pair(3.79))
    previous_energy = calculator.results["energy"]
    assert calculator.accept_geometry(pair(3.79))
    calculator.calculate(pair(3.79))
    assert calculator.results["energy"] == pytest.approx(previous_energy, abs=1e-12)
    state = calculator.bias_controller.state_dict()
    checkpoint = Path(_CONTROLLER_STATE_FILE).read_text()
    # Rejected/unaccepted probes must not replace the accepted checkpoint.
    calculator.calculate(pair(3.7))
    assert Path(_CONTROLLER_STATE_FILE).read_text() == checkpoint
    # The original process is no longer a trace writer after restart.
    calculator._trace_recorder = None
    restored = PyarGeometricCalculator(dict(params, bias_controller_restart=True), [[0], [1]])
    restored.calculate(pair(3.79))
    assert restored.bias_controller.state_dict() == state
    assert restored.results["energy"] == pytest.approx(previous_energy, abs=1e-12)
    for active in (calculator, restored):
        active.calculate(pair(3.78))
        active.accept_geometry(pair(3.78))
        active.calculate(pair(3.78))
    assert calculator.bias_controller.state_dict() == restored.bias_controller.state_dict()
    np.testing.assert_allclose(calculator.results["forces"], restored.results["forces"])
    records = load_trace_records("reaction_trace")
    assert records[-1]["bias_controller"] == restored.bias_controller.state_dict()
    assert records[-1]["bias_controller"]["alpha_max"] == restored.alpha_max
    assert records[-1]["bias_controller"]["decision"]["alpha_critical"] is not None
    assert records[-1]["bias_controller"]["decision"]["alpha_target"] is not None
    assert records[-1]["bias_parameters"]["potential"] == kind
    assert records[-1]["bias_parameters"]["gamma_kj_mol"] == 100.0
    if kind == "softmin":
        assert records[-1]["bias_parameters"]["beta_per_bohr"] == 1.0
    assert json.loads(Path("pyar_geometric_state.json").read_text())["bias_controller"] == records[-1]["bias_controller"]
    analyse_reaction_trace(Path.cwd())
    with open("path_summary.csv") as handle:
        rows = list(csv.DictReader(handle))
    assert float(rows[-1]["bias_alpha"]) == restored.bias_controller.decision.alpha
    assert float(rows[-1]["bias_alpha_max_hartree_per_bohr"]) == restored.alpha_max
    assert rows[-1]["bias_potential"] == kind
    assert float(rows[-1]["bias_gamma_kj_mol"]) == 100.0
    assert float(rows[-1]["bias_energy_offset_hartree"]) == restored.bias_controller.energy_offset


def test_restart_rejects_mismatched_geometry_and_configuration(calculation):
    calculator = PyarGeometricCalculator(calculation, [[0], [1]])
    calculator.calculate(pair(3.8))
    restored = PyarGeometricCalculator(dict(calculation, bias_controller_restart=True), [[0], [1]])
    with pytest.raises(ValueError, match="geometry"):
        restored.calculate(pair(3.7))
    with pytest.raises(ValueError, match="configuration"):
        PyarGeometricCalculator(dict(calculation, gamma=10., bias_controller_restart=True), [[0], [1]])


def test_cli_controller_parameters_reach_bias_controller(calculation):
    calculator = PyarGeometricCalculator(
        dict(calculation, bias_alpha_min=0.02, bias_alpha_margin=0.03,
             bias_alpha_smoothing=0.4, bias_alpha_epsilon=1e-9), [[0], [1]]
    )
    assert calculator.bias_controller.configuration() == {
        "policy": "adaptive",
        "alpha_min": 0.02,
        "safety_margin": 0.03,
        "smoothing": 0.4,
        "smoothing_mode": "decrease_only",
        "epsilon": 1e-9,
        "scheduled_alpha": None,
    }


@pytest.mark.parametrize("policy", ["fixed", "scheduled", "adaptive"])
def test_unbiased_relaxation_bypasses_positive_controller_minimum(calculation, policy):
    params = without_afir_bias(dict(calculation, bias_controller=policy, bias_alpha_min=0.01))
    calculator = PyarGeometricCalculator(params, [[0], [1]])
    with patch.object(calculator.bias_controller, "select", side_effect=AssertionError("disabled")):
        calculator.calculate(pair(3.8))
    assert calculator.results["energy"]/Hartree == pytest.approx(0.002)
    assert calculator.bias_controller.decision is None
    assert not Path(_CONTROLLER_STATE_FILE).exists()


def test_adaptive_command_uses_accepted_step_driver(calculation):
    molecule = SimpleNamespace(name="pair", title="pair", atoms_list=["C", "C"],
        number_of_atoms=2, charge=0, multiplicity=1, scftype="rhf",
        coordinates=pair(3.8).positions, fragments=[[0], [1]])
    with patch("pyar.backends.geometric._find_geometric_executable", return_value="geometric-optimize"):
        wrapper = Geometric(molecule, calculation)
    command = wrapper._build_command()
    assert command[1:3] == ["-m", "pyar.backends.adaptive_geometric"]
    assert json.loads(command[-1])["qc_params"]["bias_controller"] == "adaptive"


def test_tiny_margin_reports_stall_instead_of_success(calculation):
    from pyar.backends.adaptive_geometric import run_adaptive_optimization

    pair(3.8).write("start.xyz")
    arguments = dict(qc_params=dict(calculation, bias_alpha_margin=1e-7,
                                   opt_cycles=10, opt_threshold="normal"),
                     fragment_indices=[[0], [1]])
    with pytest.raises(RuntimeError, match="stalled below its force ceiling"):
        run_adaptive_optimization("start.xyz", arguments)


def test_iteration_limit_preserves_endpoint_for_parent_relaxation(calculation):
    from pyar.backends.adaptive_geometric import run_adaptive_optimization

    pair(3.8).write("start.xyz")
    arguments = dict(qc_params=dict(calculation, opt_cycles=1,
                                   opt_threshold="normal"),
                     fragment_indices=[[0], [1]])
    optimizer = run_adaptive_optimization("start.xyz", arguments)
    state = json.loads(Path("pyar_geometric_state.json").read_text())
    assert state["optimization_status"] == "cycle_exceeded"
    assert Path("start_optim.xyz").exists()
    from geometric.optimize import OPT_STATE
    assert optimizer.state == OPT_STATE.FAILED


def test_physical_energy_is_independent_of_bias_continuity_offset(calculation):
    calculator = PyarGeometricCalculator(calculation, [[0], [1]])
    wrapper = Geometric.__new__(Geometric)
    calculator.calculate(pair(3.8))
    physical = wrapper._read_final_energy()
    calculator.bias_controller.energy_offset += 100.
    calculator.calculate(pair(3.8))
    state = json.loads(Path("pyar_geometric_state.json").read_text())
    assert wrapper._read_final_energy() == physical
    assert physical == pytest.approx(0.002)
    assert state["total_energy_hartree"] > 99.
    assert state["total_energy_hartree"] == pytest.approx(
        state["backend_energy_hartree"] + state["bias_energy_hartree"])


def test_real_optimizer_updates_only_accepted_steps(calculation, monkeypatch):
    pytest.importorskip("geometric")
    from pyar.backends.adaptive_geometric import AdaptiveOptimizer, run_adaptive_optimization
    from geometric.optimize import OPT_STATE

    class RejectedTrialProvider(HarmonicProvider):
        calls = 0

        def evaluate(self, atoms, coordinates_bohr):
            self.calls += 1
            result = super().evaluate(atoms, coordinates_bohr)
            # Deliberately bad first trial to exercise geomeTRIC's rejection.
            if self.calls == 2:
                return EnergyGradientResult(result.energy_hartree + 1., result.gradient_hartree_per_bohr)
            return result

    provider = RejectedTrialProvider()
    monkeypatch.setattr("pyar.backends.geometric._resolve_backend_evaluator", lambda *args: provider)
    events = []
    evaluate_step = AdaptiveOptimizer.evaluateStep

    def observe(self):
        trial = self.X.copy()
        before = self.engine.calculator.bias_controller.state_dict()
        evaluate_step(self)
        after = self.engine.calculator.bias_controller.state_dict()
        accepted = np.array_equal(trial, self.X)
        events.append(accepted)
        if not accepted:
            assert after == before
        elif self.state != OPT_STATE.FAILED:
            assert after["segment_index"] == before["segment_index"] + 1

    monkeypatch.setattr(AdaptiveOptimizer, "evaluateStep", observe)
    pair(3.8).write("start.xyz")
    arguments = dict(qc_params=dict(calculation, opt_cycles=100, opt_threshold="tight"),
                     fragment_indices=[[0], [1]])
    optimizer = run_adaptive_optimization("start.xyz", arguments)
    assert optimizer.state == OPT_STATE.CONVERGED
    assert False in events
    assert True in events
    assert optimizer.engine.calculator.bias_controller.segment_index == sum(events)
    saved_index = optimizer.engine.calculator.bias_controller.segment_index
    expected_distance = 4. - optimizer.engine.calculator.alpha_max / 0.1
    assert np.linalg.norm(optimizer.X.reshape(-1, 3)[1] - optimizer.X.reshape(-1, 3)[0]) == pytest.approx(expected_distance, abs=0.003)
    # Resume starts from the accepted checkpoint even with the original input.
    arguments["qc_params"]["bias_controller_restart"] = True
    resumed = run_adaptive_optimization("start.xyz", arguments)
    assert resumed.state == OPT_STATE.CONVERGED
    assert resumed.engine.calculator.bias_controller.segment_index > saved_index
