"""Tests for original-solute-centred microsolvation."""

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.microsolvation.placement import (
    _has_clash,
    generate_microsolvation_candidates,
)
from pyar.microsolvation.request import MicrosolvationRequest
from pyar.microsolvation.surface import build_solute_surface
from pyar.scripts import modern_microsolvate
from pyar.state.microsolvation import MicrosolvationRunState, MicrosolvationStateError
from pyar.workflows.microsolvation import microsolvate


def methane():
    return Molecule(
        ["C", "H", "H", "H", "H"],
        [[0, 0, 0], [0.629, 0.629, 0.629], [-0.629, -0.629, 0.629],
         [-0.629, 0.629, -0.629], [0.629, -0.629, -0.629]],
        name="methane", title="methane", charge=0, multiplicity=1,
    )


def water():
    return Molecule(
        ["O", "H", "H"], [[0, 0, 0], [0.9572, 0, 0], [-0.239, 0.927, 0]],
        name="water", title="water", charge=0, multiplicity=1,
    )


def request(*, count=2, site=None, orientations=4, maximum=2):
    return MicrosolvationRequest(
        methane(), water(), count, orientations, {}, maximum, site,
        surface_points_per_atom=32, probe_radius=1.4, shell_tolerance=3.5,
    )


def test_fibonacci_surface_is_deterministic_and_exposed_points_exist():
    first = build_solute_surface(methane(), points_per_atom=32, probe_radius=1.4)
    second = build_solute_surface(methane(), points_per_atom=32, probe_radius=1.4)
    assert len(first.points) > 0
    np.testing.assert_array_equal(first.points, second.points)
    assert set(first.parent_atoms).issubset(set(range(5)))


def test_second_water_targets_original_methane_surface_and_sees_full_cluster():
    req = request()
    first = generate_microsolvation_candidates(
        req.solute, req, step=10000, solute_atom_count=len(req.solute),
    )
    assert first
    second = generate_microsolvation_candidates(
        first[0], req, step=20000, solute_atom_count=len(req.solute),
    )
    assert second
    surface = build_solute_surface(req.solute, points_per_atom=32, probe_radius=1.4)
    original_sites = {tuple(np.round(point, 10)) for point in surface.points}
    for candidate in second:
        assert candidate.microsolvation_target_atom < len(req.solute)
        assert tuple(np.round(candidate.microsolvation_target_point, 10)) in original_sites
        incoming = Molecule(candidate.atoms_list[-3:], candidate.coordinates[-3:])
        full_cluster = Molecule(candidate.atoms_list[:-3], candidate.coordinates[:-3])
        assert len(full_cluster) == len(req.solute) + 3
        assert not _has_clash(full_cluster, incoming)
        assert candidate.solute_atom_indices == tuple(range(5))
        assert len(candidate.solvent_fragments) == 2


def test_directed_site_restricts_targets_to_original_solute_atom():
    req = request(site=(3,))
    candidates = generate_microsolvation_candidates(
        req.solute, req, step=10000, solute_atom_count=len(req.solute),
    )
    assert candidates
    assert {candidate.microsolvation_target_atom for candidate in candidates} == {3}


def test_multiple_sites_use_union_of_original_solute_regions():
    req = request(site=(3, 4))
    candidates = generate_microsolvation_candidates(
        req.solute, req, step=10000, solute_atom_count=len(req.solute),
    )
    assert candidates
    assert {candidate.microsolvation_target_atom for candidate in candidates} <= {3, 4}


def test_site_validation_happens_in_request():
    with pytest.raises(ValueError, match="out of range"):
        request(site=(5,))


def test_geometry_only_workflow_is_bounded_and_records_coverage_and_identity(tmp_path):
    req = request(count=2, orientations=3, maximum=2)
    result = microsolvate(req, output=tmp_path / "microsolvation")
    assert result.status == "completed"
    assert len(result.selected_paths) <= 2
    assert result.metadata["completed_solvents"] == 2
    assert len(result.metadata["coverage_history"]) == 2
    assert (tmp_path / "microsolvation" / "step_001" / "coverage.json").is_file()
    state = MicrosolvationRunState.load(tmp_path / "microsolvation", req.to_state_dict())
    assert state is not None
    restored = state.restore_pool()
    assert all(molecule.energy is None for molecule in restored)
    assert all(molecule.solute_atom_indices == tuple(range(5)) for molecule in restored)
    assert all(len(molecule.solvent_fragments) == 2 for molecule in restored)


def test_state_rejects_changed_request(tmp_path):
    first = request(count=1)
    microsolvate(first, output=tmp_path / "microsolvation")
    with pytest.raises(MicrosolvationStateError, match="differs"):
        MicrosolvationRunState.load(tmp_path / "microsolvation", request(count=2).to_state_dict())


def test_restart_resumes_at_next_unfinished_solvent_step(tmp_path, monkeypatch):
    req = request(count=2, orientations=3, maximum=2)
    import pyar.workflows.microsolvation as workflow
    original = workflow.generate_microsolvation_candidates
    calls = []

    def interrupt_second_step(seed, request, *, step, solute_atom_count):
        calls.append(step)
        if step >= 20000:
            raise RuntimeError("simulated interruption")
        return original(seed, request, step=step, solute_atom_count=solute_atom_count)

    monkeypatch.setattr(workflow, "generate_microsolvation_candidates", interrupt_second_step)
    with pytest.raises(RuntimeError, match="simulated interruption"):
        workflow.microsolvate(req, output=tmp_path / "microsolvation")
    state = MicrosolvationRunState.load(tmp_path / "microsolvation", req.to_state_dict())
    assert state.data["next_step"] == 2
    assert len(state.data["completed_steps"]) == 1

    calls.clear()

    def record_resume(seed, request, *, step, solute_atom_count):
        calls.append(step)
        return original(seed, request, step=step, solute_atom_count=solute_atom_count)

    monkeypatch.setattr(workflow, "generate_microsolvation_candidates", record_resume)
    result = workflow.microsolvate(req, output=tmp_path / "microsolvation")
    assert result.status == "completed"
    assert len(result.metadata["coverage_history"]) == 2
    assert calls and all(step >= 20000 for step in calls)


def test_modified_completed_snapshot_blocks_restart(tmp_path):
    req = request(count=1, orientations=3, maximum=2)
    microsolvate(req, output=tmp_path / "microsolvation")
    state = MicrosolvationRunState.load(tmp_path / "microsolvation", req.to_state_dict())
    snapshot = state.directory / state.data["completed_steps"][0]["selected_paths"][0]
    snapshot.write_text("tampered\n")
    with pytest.raises(MicrosolvationStateError, match="modified"):
        MicrosolvationRunState.load(tmp_path / "microsolvation", req.to_state_dict())


def _write_xyz(path, molecule):
    molecule.mol_to_xyz(str(path))


def test_check_validates_without_creating_output_or_running_workflow(tmp_path, monkeypatch, capsys):
    solute_path, solvent_path = tmp_path / "methane.xyz", tmp_path / "water.xyz"
    _write_xyz(solute_path, methane())
    _write_xyz(solvent_path, water())
    output = tmp_path / "microsolvation"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(modern_microsolvate, "microsolvate", lambda *a, **k: pytest.fail("workflow called"))
    modern_microsolvate.main([str(solute_path), str(solvent_path), "--count", "2", "--check", "--output", str(output)])
    assert not output.exists()
    assert "no calculations were performed" in capsys.readouterr().out


def test_xtb_resolution_uses_modern_model_without_dft_defaults(tmp_path):
    solute_path, solvent_path = tmp_path / "methane.xyz", tmp_path / "water.xyz"
    _write_xyz(solute_path, methane())
    _write_xyz(solvent_path, water())
    args = modern_microsolvate.build_parser().parse_args(
        [str(solute_path), str(solvent_path), "--count", "1", "--backend", "xtb"]
    )
    resolved = modern_microsolvate.resolve_request(args)
    assert resolved.backend_parameters["software"] == "xtb"
    assert resolved.backend_parameters["xtb_model"] == "gfn2"
    assert resolved.backend_parameters["method"] is None
    assert resolved.backend_parameters["basis"] is None
