import json

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.evaluation_budget import EvaluationBudgetExhausted, reserve_backend_evaluation
from pyar.microsolvation.perturbation import (
    has_interfragment_clash, propose_rigid_move, preserves_cluster_identity,
)
from pyar.microsolvation.request import MicrosolvationRequest, MicrosolvationRequestError


def cluster():
    molecule = Molecule(["C", "H", "H", "H", "H", "O", "H", "H"],
        [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63],
         [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63],
         [4, 0, 0], [4.76, 0.59, 0], [3.24, 0.59, 0]], name="cluster")
    molecule.solute_atom_indices = (0, 1, 2, 3, 4)
    molecule.solvent_fragments = ((5, 6, 7),)
    molecule.charge = 0
    molecule.multiplicity = 1
    return molecule


def test_rigid_move_is_seeded_and_preserves_fragment_identity_and_composition():
    original = cluster()
    first, description = propose_rigid_move(original, np.random.Generator(np.random.PCG64(17)))
    second, second_description = propose_rigid_move(original, np.random.Generator(np.random.PCG64(17)))
    assert description == second_description
    assert np.array_equal(first.coordinates, second.coordinates)
    assert first.atoms_list == original.atoms_list
    assert first.solute_atom_indices == original.solute_atom_indices
    assert first.solvent_fragments == original.solvent_fragments
    assert first.charge == 0 and first.multiplicity == 1
    assert np.array_equal(first.coordinates[:5], original.coordinates[:5])


def test_solvent_site_exchange_swaps_local_geometries_at_fixed_sites():
    molecule = Molecule(["He", "O", "H", "H", "O", "H", "H"],
        [[-5, 0, 0], [0, 0, 0], [.96, 0, 0], [-.24, .93, 0],
         [5, 0, 0], [5, .96, 0], [4.07, -.24, 0]])
    molecule.solute_atom_indices = (0,)
    molecule.solvent_fragments = ((1, 2, 3), (4, 5, 6))
    molecule.charge = 0
    molecule.multiplicity = 1
    moved, description = propose_rigid_move(molecule, np.random.Generator(np.random.PCG64(3)))
    assert description["move"] == "solvent-site-exchange"
    masses = np.asarray(molecule.atomic_mass)
    for fragment in ((1, 2, 3), (4, 5, 6)):
        assert np.allclose(np.average(moved.coordinates[list(fragment)], axis=0,
                                      weights=masses[list(fragment)]),
                           np.average(molecule.coordinates[list(fragment)], axis=0,
                                      weights=masses[list(fragment)]))
    assert not np.array_equal(moved.coordinates, molecule.coordinates)
    assert moved.atoms_list == molecule.atoms_list
    assert moved.solvent_fragments == molecule.solvent_fragments


def test_multi_fragment_move_preserves_each_solvent_geometry_and_mapping():
    atoms = ["He"]
    coordinates = [[-5, 0, 0]]
    fragments = []
    for index, x in enumerate((0., 6., 12., 18.)):
        start = len(atoms)
        atoms.extend(["O", "H", "H"])
        coordinates.extend([[x, 0, 0], [x + .96, 0, 0], [x - .24, .93, 0]])
        fragments.append((start, start + 1, start + 2))
    molecule = Molecule(atoms, coordinates)
    molecule.solute_atom_indices = (0,)
    molecule.solvent_fragments = tuple(fragments)
    molecule.charge = 0
    molecule.multiplicity = 1
    moved, description = propose_rigid_move(molecule, np.random.Generator(np.random.PCG64(1)))
    assert description["move"] == "multi-fragment-rigid"
    assert len(description["fragments"]) == 3
    assert np.array_equal(moved.coordinates[0], molecule.coordinates[0])
    assert moved.atoms_list == molecule.atoms_list
    assert moved.solvent_fragments == molecule.solvent_fragments
    for fragment in fragments:
        indices = list(fragment)
        before = np.linalg.norm(molecule.coordinates[indices[0]] - molecule.coordinates[indices[1]])
        after = np.linalg.norm(moved.coordinates[indices[0]] - moved.coordinates[indices[1]])
        assert after == pytest.approx(before)


def test_interfragment_clash_detects_overlap_but_accepts_separated_fragments():
    molecule = cluster()
    assert not has_interfragment_clash(molecule)
    molecule.coordinates[5:] -= [4, 0, 0]
    assert has_interfragment_clash(molecule)


def test_optimized_cluster_identity_rejects_atom_reordering_and_solute_graph_change():
    molecule = cluster()
    solute = Molecule(molecule.atoms_list[:5], molecule.coordinates[:5], charge=0, multiplicity=1)
    solvent = Molecule(molecule.atoms_list[5:], molecule.coordinates[5:], charge=0, multiplicity=1)
    request = MicrosolvationRequest(solute, solvent, 1)
    assert preserves_cluster_identity(molecule, request, 1)
    reordered = molecule.copy()
    reordered.atoms_list[5], reordered.atoms_list[6] = reordered.atoms_list[6], reordered.atoms_list[5]
    assert not preserves_cluster_identity(reordered, request, 1)
    broken = molecule.copy()
    broken.solute_atom_indices = molecule.solute_atom_indices
    broken.solvent_fragments = molecule.solvent_fragments
    broken.coordinates[1] += [20, 0, 0]
    assert not preserves_cluster_identity(broken, request, 1)


def test_backend_evaluation_budget_is_strict_and_persistent(tmp_path):
    counter = tmp_path / "counter.json"
    counter.write_text(json.dumps({"used": 0, "limit": 2, "deadline_epoch": 4_000_000_000}))
    assert reserve_backend_evaluation(counter) == 1
    assert reserve_backend_evaluation(counter) == 2
    with pytest.raises(EvaluationBudgetExhausted):
        reserve_backend_evaluation(counter)
    assert json.loads(counter.read_text())["used"] == 2


def test_backend_evaluation_budget_checks_deadline_before_reserving(tmp_path):
    counter = tmp_path / "counter.json"
    counter.write_text(json.dumps({"used": 0, "limit": 2, "deadline_epoch": 1}))
    with pytest.raises(EvaluationBudgetExhausted, match="time limit"):
        reserve_backend_evaluation(counter)
    assert json.loads(counter.read_text())["used"] == 0


def test_ms4_request_requires_cluster_backend_and_explicit_hard_budgets():
    solute = Molecule(["He"], [[0, 0, 0]], charge=0, multiplicity=1)
    solvent = Molecule(["He"], [[0, 0, 0]], charge=0, multiplicity=1)
    with pytest.raises(MicrosolvationRequestError, match="cluster"):
        MicrosolvationRequest(solute, solvent, 1, resampling_moves=2,
                              resampling_evaluation_budget=20, resampling_time_limit=30)
    request = MicrosolvationRequest(solute, solvent, 1, objective="cluster", placement_score="energy",
        backend_parameters={"software": "xtb", "geometry_optimizer": "geometric"},
        resampling_moves=2, resampling_evaluation_budget=20, resampling_time_limit=30, resampling_seed=9)
    assert request.to_state_dict()["resampling_policy"]["rng"] == "numpy-PCG64"


def test_v2_resampling_stage_persists_proposals_and_is_idempotent(tmp_path, monkeypatch):
    import importlib
    workflow = importlib.import_module("pyar.workflows.microsolvation")
    from pyar.persistence import read_json

    solute = Molecule(["He"], [[0, 0, 0]], name="helium", charge=0, multiplicity=1)
    solvent = Molecule(["He"], [[0, 0, 0]], name="helium", charge=0, multiplicity=1)
    request = MicrosolvationRequest(solute, solvent, 1, 2, {"software": "xtb", "geometry_optimizer": "geometric"},
        2, surface_points_per_atom=8, objective="cluster", placement_score="energy",
        resampling_moves=2, resampling_evaluation_budget=20, resampling_time_limit=30, resampling_seed=12)

    monkeypatch.setattr(workflow, "_energy_prescreen", lambda candidates, req: (candidates, []))
    def optimize(candidates, req, directory, *, evaluation_budget_path=None):
        for candidate in candidates:
            candidate.energy = 0.0
        return list(candidates), 0, len(candidates), []
    monkeypatch.setattr(workflow, "_optimize_candidates", optimize)
    root = tmp_path / "ms4"
    result = workflow.microsolvate(request, output=root)
    assert result.status == "completed"
    state = read_json(root / "state.json")
    assert len(state["resampling"]["attempts"]) == 2
    proposal_files = sorted((root / "resampling/attempts").glob("attempt_*/proposal.json"))
    assert len(proposal_files) == 2
    for attempt in state["resampling"]["attempts"]:
        if "optimized_candidate" in attempt:
            assert attempt["optimized_candidate"]["role"] == "resampling-candidate"
            assert (root / attempt["optimized_candidate"]["path"]).is_file()
    before = {path: path.read_bytes() for path in (root / "resampling").rglob("*") if path.is_file()}
    workflow.microsolvate(request, output=root)
    assert before == {path: path.read_bytes() for path in before}


def test_interrupted_resampling_never_replays_committed_proposal(tmp_path, monkeypatch):
    import importlib
    workflow = importlib.import_module("pyar.workflows.microsolvation")
    from pyar.persistence import read_json

    solute = Molecule(["He"], [[0, 0, 0]], name="helium", charge=0, multiplicity=1)
    solvent = Molecule(["He"], [[0, 0, 0]], name="helium", charge=0, multiplicity=1)
    request = MicrosolvationRequest(solute, solvent, 1, 2, {"software": "xtb", "geometry_optimizer": "geometric"},
        2, surface_points_per_atom=8, objective="cluster", placement_score="energy",
        resampling_moves=1, resampling_evaluation_budget=20, resampling_time_limit=30, resampling_seed=12)
    monkeypatch.setattr(workflow, "_energy_prescreen", lambda candidates, req: (candidates, []))
    calls = []
    def crash_during_resampling(candidates, req, directory, *, evaluation_budget_path=None):
        calls.append(directory)
        if evaluation_budget_path is not None:
            raise RuntimeError("simulated process interruption")
        for candidate in candidates:
            candidate.energy = 0.0
        return list(candidates), 0, len(candidates), []
    monkeypatch.setattr(workflow, "_optimize_candidates", crash_during_resampling)
    root = tmp_path / "interrupted"
    with pytest.raises(RuntimeError, match="interruption"):
        workflow.microsolvate(request, output=root)
    proposal = (root / "resampling/attempts/attempt_00000/proposal.json").read_bytes()
    def optimize(candidates, req, directory, *, evaluation_budget_path=None):
        for candidate in candidates:
            candidate.energy = 0.0
        return list(candidates), 0, len(candidates), []
    monkeypatch.setattr(workflow, "_optimize_candidates", optimize)
    result = workflow.microsolvate(request, output=root)
    state = read_json(root / "state.json")
    assert result.status == "completed"
    assert len(state["resampling"]["attempts"]) == 1
    assert state["resampling"]["attempts"][0]["status"] == "interrupted"
    assert (root / "resampling/attempts/attempt_00000/proposal.json").read_bytes() == proposal
