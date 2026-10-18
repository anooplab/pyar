"""Objective-aware microsolvation selection and persistence contracts."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.microsolvation.placement import generate_microsolvation_candidates
from pyar.microsolvation.request import MicrosolvationRequest
from pyar.persistence import read_json
from pyar.state.microsolvation import MicrosolvationRunState
from pyar.workflows.microsolvation import (
    _deduplicate_and_select, _energy_prescreen, _optimize_candidates, _quality_key, microsolvate,
)


def methane():
    return Molecule(["C", "H", "H", "H", "H"],
                    [[0, 0, 0], [0.629, 0.629, 0.629], [-0.629, -0.629, 0.629],
                     [-0.629, 0.629, -0.629], [0.629, -0.629, -0.629]],
                    name="methane", title="methane", charge=0, multiplicity=1)


def water():
    return Molecule(["O", "H", "H"], [[0, 0, 0], [0.9572, 0, 0], [-0.239, 0.927, 0]],
                    name="water", title="water", charge=0, multiplicity=1)


def request(*, objective="shell", score="geometric", backend=None, count=2, maximum=2):
    settings = {} if backend is None else {"software": backend}
    return MicrosolvationRequest(
        methane(), water(), count, 4, settings, maximum, None, 32, 1.4, 3.5,
        objective, score,
    )


def test_historical_default_request_serialization_is_unchanged():
    payload = request().to_state_dict()
    assert payload["selection_policy"]["name"] == "surface-coverage-energy-diversity-v1"
    assert "objective" not in payload["selection_policy"]
    assert "placement_policy" not in payload
    assert "site_contact_policy" not in payload


def test_default_shell_quality_key_retains_historical_energy_tiebreaker():
    req = request()
    tied = [
        SimpleNamespace(name="higher-energy", energy=-10.0,
                        microsolvation_coverage={"coverage_fraction": 0.5}),
        SimpleNamespace(name="lower-energy", energy=-10.01,
                        microsolvation_coverage={"coverage_fraction": 0.5}),
    ]
    assert [molecule.name for molecule in sorted(tied, key=lambda item: _quality_key(item, req))] == [
        "lower-energy", "higher-energy",
    ]


def test_auto_wall_selection_policy_is_separate_from_confinement():
    qc = {"software": "xtb", "geometry_optimizer": "geometric"}
    auto = MicrosolvationRequest(methane(), water(), 1, 4, qc, 2, None, 32, 1.4, 3.5,
                                 "shell", "geometric", (), None, "auto")
    explicit = MicrosolvationRequest(methane(), water(), 1, 4, qc, 2, None, 32, 1.4, 3.5,
                                     "shell", "geometric", (), None, "outer-sphere", 4.0, 0.5)
    auto_state = auto.to_state_dict()
    explicit_state = explicit.to_state_dict()
    assert auto_state["selection_policy"]["name"] == "objective-aware-v2"
    assert auto_state["selection_policy"]["ranking"] == "solute-surface-coverage-then-physical-energy"
    assert auto_state["confinement"]["mode"] == "case-dependent-outer-sphere-v1"
    # The legacy explicit-wall policy remains byte-compatible for restarts.
    assert explicit_state["selection_policy"] == {
        "name": "objective-aware-v1", "objective": "shell", "placement_score": "geometric",
        "energy_unit": "hartree", "ranking": "unbiased-post-relaxation-energy",
        "cluster_energy_frontier_hartree": None, "unbiased_energy_frontier_hartree": 0.005,
        "prescreen_pool_multiplier": 4,
    }
    tied_candidates = [
        SimpleNamespace(name="higher-energy", energy=-10.0,
                        microsolvation_coverage={"coverage_fraction": 0.5}),
        SimpleNamespace(name="lower-energy", energy=-10.01,
                        microsolvation_coverage={"coverage_fraction": 0.5}),
    ]
    assert [_quality_key(candidate, request()) for candidate in tied_candidates] == [
        _quality_key(candidate, auto) for candidate in tied_candidates
    ]


def test_auto_wall_without_candidate_wall_use_keeps_shell_selection_order():
    qc = {"software": "xtb", "geometry_optimizer": "geometric"}
    auto = MicrosolvationRequest(methane(), water(), 1, 4, qc, 1, None, 32, 1.4, 3.5,
                                 "shell", "geometric", (), None, "auto")
    selected_names = []
    for req in (request(), auto):
        candidates = []
        for name, energy, coverage, shift in (
            ("high-coverage", -10.0, 0.9, [3.0, 0.0, 0.0]),
            ("low-energy", -10.03, 0.2, [-3.0, 0.0, 0.0]),
        ):
            molecule = req.solute.merged_with(req.solvent)
            molecule.name = name
            molecule.solute_atom_indices = tuple(range(5))
            molecule.solvent_fragments = (tuple(range(5, 8)),)
            molecule.fragments = [list(molecule.solute_atom_indices), list(molecule.solvent_fragments[0])]
            molecule.fragments_coordinates = molecule.split_coordinates()
            molecule.fragments_atoms_list = molecule.split_atoms_lists()
            molecule.coordinates[5:] += np.asarray(shift)
            molecule.energy = energy
            molecule.microsolvation_coverage = {"coverage_fraction": coverage, "shell_valid": True}
            candidates.append(molecule)
        selected, diagnostics = _deduplicate_and_select(candidates, 1, req)
        selected_names.append(selected[0].name)
        assert diagnostics["selection_mode"] == "solute-surface-coverage-then-energy-and-structural-diversity"
    assert selected_names == ["high-coverage", "high-coverage"]


def test_energy_and_cluster_policy_require_a_real_energy_provider():
    with pytest.raises(ValueError, match="requires --backend"):
        request(score="energy")
    with pytest.raises(ValueError, match="requires --placement-score energy"):
        request(objective="cluster", backend="xtb")


def test_cluster_growth_targets_existing_solvent_after_first_step():
    req = request(objective="cluster", score="energy", backend="xtb")
    first = generate_microsolvation_candidates(req.solute, req, step=10000,
                                               solute_atom_count=len(req.solute))
    assert first and first[0].microsolvation_target_kind == "solute"
    second = generate_microsolvation_candidates(first[0], req, step=20000,
                                                solute_atom_count=len(req.solute))
    assert second
    assert {candidate.microsolvation_target_kind for candidate in second} == {"solute", "solvent"}
    assert any(candidate.microsolvation_target_fragment == 0 for candidate in second)


def test_cluster_energy_ranking_keeps_low_energy_low_coverage_candidate():
    req = request(objective="cluster", score="energy", backend="xtb", count=1, maximum=1)
    candidates = []
    for name, energy, coverage, oxygen in (
        ("covered-high-energy", -10.0, 0.9, [3.0, 0.0, 0.0]),
        ("droplet-low-energy", -10.02, 0.1, [-3.0, 0.0, 0.0]),
    ):
        molecule = req.solute.merged_with(req.solvent)
        molecule.name = name
        molecule.solute_atom_indices = tuple(range(5))
        molecule.solvent_fragments = (tuple(range(5, 8)),)
        molecule.fragments = [list(molecule.solute_atom_indices), list(molecule.solvent_fragments[0])]
        molecule.fragments_coordinates = molecule.split_coordinates()
        molecule.fragments_atoms_list = molecule.split_atoms_lists()
        molecule.coordinates[5:] += np.asarray(oxygen)
        molecule.energy = energy
        molecule.microsolvation_coverage = {"coverage_fraction": coverage, "shell_valid": coverage > 0.5}
        candidates.append(molecule)
    selected, diagnostics = _deduplicate_and_select(candidates, 1, req)
    assert selected[0].name == "droplet-low-energy"
    assert diagnostics["selection_mode"] == "fixed-composition-energy-then-diversity"


def test_wall_shell_ranking_uses_unbiased_energy_after_shell_filter():
    qc = {"software": "xtb", "geometry_optimizer": "geometric"}
    req = MicrosolvationRequest(methane(), water(), 1, 4, qc, 2, None, 32, 1.4, 3.5,
                                "shell", "geometric", (), None, "outer-sphere", 4.0, 0.5)
    candidates = []
    for name, energy, oxygen, valid in (
        ("covered-higher-energy", -10.0, [3.0, 0.0, 0.0], True),
        ("lower-energy", -10.01, [-3.0, 0.0, 0.0], True),
    ):
        molecule = req.solute.merged_with(req.solvent)
        molecule.name = name
        molecule.solute_atom_indices = tuple(range(5))
        molecule.solvent_fragments = (tuple(range(5, 8)),)
        molecule.fragments = [list(molecule.solute_atom_indices), list(molecule.solvent_fragments[0])]
        molecule.fragments_coordinates = molecule.split_coordinates()
        molecule.fragments_atoms_list = molecule.split_atoms_lists()
        molecule.coordinates[5:] += np.asarray(oxygen)
        molecule.energy = energy
        molecule.microsolvation_coverage = {"coverage_fraction": 0.5 if valid else 0.0,
                                            "shell_valid": valid}
        if valid:
            candidates.append(molecule)
    selected, diagnostics = _deduplicate_and_select(candidates, 1, req)
    assert selected[0].name == "lower-energy"
    assert diagnostics["selection_mode"] == "unbiased-post-relaxation-energy-then-diversity"
    assert diagnostics["energy_unit"] == "hartree"


def test_auto_wall_policy_derives_case_specific_radius_and_thermal_scale_strength():
    qc = {"software": "xtb", "geometry_optimizer": "geometric"}
    req = MicrosolvationRequest(methane(), water(), 1, 4, qc, 2, None, 32, 1.4, 3.5,
                                "shell", "geometric", (), None, "auto")
    confinement = req.to_state_dict()["confinement"]
    assert confinement["mode"] == "case-dependent-outer-sphere-v1"
    expected_radius = (
        confinement["radius_derivation"]["solute_surface_radius_angstrom"]
        + confinement["radius_derivation"]["solvent_fragment_com_extent_angstrom"]
        + 0.5
    )
    assert confinement["radius_angstrom"] == pytest.approx(expected_radius)
    strength = confinement["force_constant_hartree_per_angstrom2"]
    assert 0.01 < strength < 0.04
    energy_at_half_angstrom = 0.5 * strength * 0.5**2
    assert energy_at_half_angstrom == pytest.approx(confinement["strength_derivation"]["target_penalty_hartree"])
    assert confinement["trigger"] == "apply-only-after-unbiased-optimization-fails-shell-validity"


def test_auto_wall_radius_expands_with_solvent_fragment_extent():
    qc = {"software": "xtb", "geometry_optimizer": "geometric"}
    ordinary = MicrosolvationRequest(methane(), water(), 1, 4, qc, 2, None, 32, 1.4, 3.5,
                                    "shell", "geometric", (), None, "auto")
    extended_water = water()
    extended_water.coordinates[1:] *= 1.5
    larger_case = MicrosolvationRequest(methane(), extended_water, 1, 4, qc, 2, None, 32, 1.4, 3.5,
                                        "shell", "geometric", (), None, "auto")
    ordinary_policy = ordinary.to_state_dict()["confinement"]
    larger_policy = larger_case.to_state_dict()["confinement"]
    assert (larger_policy["radius_angstrom"] > ordinary_policy["radius_angstrom"])
    assert (larger_policy["force_constant_hartree_per_angstrom2"]
            == ordinary_policy["force_constant_hartree_per_angstrom2"])


@pytest.mark.parametrize("escaped,expected_calls", [(False, 1), (True, 3)])
def test_auto_wall_optimizes_unbiased_first_and_only_rescues_escaped_candidate(
        tmp_path, monkeypatch, escaped, expected_calls):
    from pathlib import Path
    from pyar.microsolvation.placement import generate_microsolvation_candidates

    qc = {"software": "xtb", "geometry_optimizer": "geometric"}
    req = MicrosolvationRequest(methane(), water(), 1, 3, qc, 1, None, 32, 1.4, 3.5,
                                "shell", "geometric", (), None, "auto")
    candidate = generate_microsolvation_candidates(
        req.solute, req, step=10000, solute_atom_count=len(req.solute))[0]
    candidate.name = "auto-escape" if escaped else "auto-contained"
    calls = []

    def fake_optimize(molecule, parameters):
        calls.append((Path.cwd().name, "microsolvation_wall" in parameters))
        if len(calls) == 1 and escaped:
            molecule.coordinates[5:] += [20.0, 0.0, 0.0]
        elif "microsolvation_wall" in parameters:
            assert calls[0][1] is False
            molecule.coordinates[5:] -= [20.0, 0.0, 0.0]
            state_dir = Path.cwd() / f"job_{molecule.name}"
            state_dir.mkdir(parents=True, exist_ok=True)
            wall = parameters["microsolvation_wall"]
            (state_dir / "pyar_geometric_state.json").write_text(
                '{"total_energy_hartree": -9.0, "wall_energy_hartree": 0.0, '
                '"wall_energy_max_hartree": 0.1}')
            assert wall["radius_angstrom"] == pytest.approx(
                req.to_state_dict()["confinement"]["radius_angstrom"])
        molecule.energy = -10.0 + len(calls) * 1e-4
        return "converged"

    monkeypatch.setattr("pyar.optimiser.optimise", fake_optimize)
    retained, rejected, optimized, failures = _optimize_candidates(
        [candidate], req, tmp_path / "candidates")
    assert not failures and rejected == 0 and optimized == 1
    assert len(retained) == 1 and retained[0].microsolvation_coverage["shell_valid"]
    assert len(calls) == expected_calls
    assert retained[0].microsolvation_wall_applied is escaped
    assert retained[0].microsolvation_wall_decision == (
        "applied-after-shell-escape" if escaped else "not-applied-unbiased-geometry-shell-valid")
    reports = retained[0].microsolvation_optimization_diagnostics
    assert len(reports) == expected_calls
    assert all(item["status"] == "INVALID_RESULT" for item in reports)
    assert all(not Path(item["diagnostic_path"]).is_absolute() for item in reports)
    assert {item["stage"] for item in reports} == (
        {"unbiased_screen", "wall_search", "unbiased_relaxation"} if escaped
        else {"unbiased_screen"})
    assert all((tmp_path / "candidates" / item["diagnostic_path"]).is_file() for item in reports)


def test_auto_wall_policy_is_inactive_for_cluster_objective():
    qc = {"software": "xtb"}
    req = MicrosolvationRequest(methane(), water(), 1, 4, qc, 1, None, 32, 1.4, 3.5,
                                "cluster", "energy", (), None, "auto")
    assert req.to_state_dict()["confinement"] == {
        "mode": "none; finite-cluster objective",
        "requested_policy": "auto-disabled-for-cluster-objective",
    }
    assert req.to_state_dict()["selection_policy"]["ranking"] == "fixed-composition-energy"


def test_energy_prescreen_keeps_low_energy_pool_and_prescreen_failures(monkeypatch):
    req = request(objective="cluster", score="energy", backend="xtb", maximum=1)
    candidates = []
    for index in range(8):
        molecule = req.solute.merged_with(req.solvent)
        molecule.name = f"candidate-{index}"
        molecule.coordinates[0, 0] = float(index)
        molecule.microsolvation_coverage = {"coverage_fraction": index / 8}
        candidates.append(molecule)

    class Provider:
        def evaluate(self, _molecule, coordinates_bohr):
            first_x_angstrom = float(coordinates_bohr[0, 0]) * 0.52917726
            if first_x_angstrom >= 7:
                raise RuntimeError("simulated single-point failure")
            return SimpleNamespace(energy_hartree=first_x_angstrom)

    monkeypatch.setattr("pyar.energy_gradient_providers.get_energy_gradient_provider",
                        lambda *_args, **_kwargs: Provider())
    retained, records = _energy_prescreen(candidates, req)
    assert {molecule.name for molecule in retained} == {
        "candidate-0", "candidate-1", "candidate-2", "candidate-3", "candidate-7",
    }
    failure = next(record for record in records if record["name"] == "candidate-7")
    assert failure["status"] == "failed"
    assert failure["selected_for_optimization"] is True


def test_cluster_v2_persists_non_shell_survivor_and_rejects_policy_change(tmp_path, monkeypatch):
    req = request(objective="cluster", score="energy", backend="xtb", count=2, maximum=1)

    class Provider:
        def evaluate(self, molecule, coordinates_bohr):
            return SimpleNamespace(energy_hartree=-10.0 - len(molecule) / 1000)

    monkeypatch.setattr("pyar.energy_gradient_providers.get_energy_gradient_provider",
                        lambda *_args, **_kwargs: Provider())

    def optimize(molecule, _qc):
        if molecule.name.endswith("_000"):
            return False
        molecule.coordinates[5:] += [45.0, 0.0, 0.0]
        molecule.energy = (-20.0 if getattr(molecule, "microsolvation_target_kind", "solute") == "solvent"
                           else -10.0) - len(molecule) / 1000
        return "converged"

    monkeypatch.setattr("pyar.optimiser.optimise", optimize)
    result = microsolvate(req, output=tmp_path / "cluster")
    assert result.status == "completed"
    state = MicrosolvationRunState.load(tmp_path / "cluster", req.to_state_dict())
    assert state is not None
    final = state.restore_pool()[0]
    assert not final.microsolvation_coverage["shell_valid"]
    assert final.microsolvation_target_kind == "solvent"
    committed = read_json(tmp_path / "cluster/steps/step_002/selection.json")
    assert committed["validation"]["shell_filtered"] is False
    assert committed["step_record"]["selection"]["selection_basis"] == "fixed-composition-energy"
    first_selection = read_json(tmp_path / "cluster/steps/step_001/selection.json")
    failures = first_selection["step_record"]["selection"]["optimization_failures"]
    assert failures and failures[0]["name"] == "micro_10000_000"
    assert (tmp_path / "cluster/steps/step_001/attempts/attempt_000/execution/candidates"
            / "micro_10000_000.xyz").is_file()

    shell = request(count=2, maximum=1)
    with pytest.raises(Exception, match="differs"):
        MicrosolvationRunState.load(tmp_path / "cluster", shell.to_state_dict())
