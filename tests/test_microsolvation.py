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


def test_outer_sphere_wall_requires_geometric_backend_and_changes_fingerprint():
    base = request(count=1)
    assert "wall" not in base.to_state_dict()["confinement"]
    backend = {"software": "xtb", "geometry_optimizer": "geometric"}
    walled = MicrosolvationRequest(
        methane(), water(), 1, 4, backend, 2, None, 32, 1.4, 3.5,
        "shell", "geometric", (), None, "outer-sphere", 8.0, 0.5,
    )
    assert walled.to_state_dict()["confinement"]["radius_angstrom"] == 8.0
    with pytest.raises(ValueError, match="geometry-optimizer geometric"):
        MicrosolvationRequest(
            methane(), water(), 1, backend_parameters={"software": "xtb"},
            wall="outer-sphere", wall_radius=8.0, wall_force_constant=0.5,
        )
    with pytest.raises(ValueError, match="only supported for shell"):
        MicrosolvationRequest(
            methane(), water(), 1, backend_parameters=backend, objective="cluster",
            wall="outer-sphere", wall_radius=8.0, wall_force_constant=0.5,
        )


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (
            {"wall": "outer-sphere", "wall_radius": 8.0, "wall_force_constant": 0.5},
            "outer-sphere wall requires an energy-gradient backend",
        ),
        (
            {
                "wall": "outer-sphere",
                "wall_radius": 8.0,
                "wall_force_constant": 0.5,
                "backend_parameters": {"software": "mopac", "geometry_optimizer": "geometric"},
            },
            "outer-sphere wall requires a registered energy-gradient provider",
        ),
        (
            {"wall": "auto"},
            "automatic wall policy requires an energy-gradient backend",
        ),
        (
            {"wall": "auto", "backend_parameters": {"software": "mopac", "geometry_optimizer": "geometric"}},
            "automatic wall policy requires a registered energy-gradient provider",
        ),
        (
            {"placement_score": "energy", "backend_parameters": {"software": "mopac"}},
            "energy placement scoring requires a registered energy-gradient provider",
        ),
    ],
)
def test_energy_gradient_backend_requirement_messages(options, message):
    with pytest.raises(ValueError, match=message):
        MicrosolvationRequest(methane(), water(), 1, **options)


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (
            {
                "wall": "outer-sphere",
                "wall_radius": 0,
                "wall_force_constant": 0,
                "objective": "cluster",
            },
            "outer-sphere wall is only supported for shell objective",
        ),
        (
            {"wall": "outer-sphere", "wall_radius": 0, "wall_force_constant": 0},
            "outer-sphere wall requires a positive --wall-radius",
        ),
        (
            {
                "wall": "outer-sphere",
                "wall_radius": 8,
                "wall_force_constant": 0.5,
                "backend_parameters": {"geometry_optimizer": "native"},
            },
            "outer-sphere wall requires an energy-gradient backend",
        ),
    ],
)
def test_wall_policy_preserves_first_error_order(options, message):
    with pytest.raises(ValueError, match=message):
        MicrosolvationRequest(methane(), water(), 1, **options)


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"wall": "unknown", "placement_score": "unknown"},
         "wall must be 'none', 'outer-sphere', or 'auto'"),
        ({"wall": "none", "wall_radius": 8.0},
         "--wall-radius and --wall-force-constant require --wall outer-sphere"),
        ({"wall": "auto", "wall_force_constant": 0.5},
         "automatic wall policy derives its own radius and strength"),
    ],
)
def test_wall_mode_validation_preserves_error_order(options, message):
    with pytest.raises(ValueError, match=message):
        MicrosolvationRequest(methane(), water(), 1, **options)


def test_outer_sphere_wall_force_is_analytic_and_translation_invariant():
    from pyar.microsolvation.wall import outer_sphere_energy_gradient

    xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [4.0, 0.0, 0.0]])
    masses = np.array([12.0, 1.0, 16.0])
    kwargs = dict(radius_angstrom=2.0, force_constant_hartree_per_angstrom2=0.5,
                  solvent_fragments=[[2]])
    energy, gradient = outer_sphere_energy_gradient(xyz, masses, [0, 1], **kwargs)
    assert energy > 0
    np.testing.assert_allclose(gradient.sum(axis=0), 0.0, atol=1e-14)
    delta = 1e-6
    shifted = xyz.copy()
    shifted[2, 0] += delta
    plus = outer_sphere_energy_gradient(shifted, masses, [0, 1], **kwargs)[0]
    shifted[2, 0] -= 2 * delta
    minus = outer_sphere_energy_gradient(shifted, masses, [0, 1], **kwargs)[0]
    assert gradient[2, 0] == pytest.approx((plus - minus) / (2 * delta), rel=1e-6)
    translated = xyz + [12.3, -4.0, 9.1]
    moved_energy, moved_gradient = outer_sphere_energy_gradient(translated, masses, [0, 1], **kwargs)
    assert moved_energy == pytest.approx(energy)
    np.testing.assert_allclose(moved_gradient, gradient, atol=1e-14)
    angle = 0.73
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    rotated_energy, rotated_gradient = outer_sphere_energy_gradient(xyz @ rotation.T, masses, [0, 1], **kwargs)
    assert rotated_energy == pytest.approx(energy)
    np.testing.assert_allclose(rotated_gradient, gradient @ rotation.T, atol=1e-14)
    zero_energy, zero_gradient = outer_sphere_energy_gradient(
        xyz, masses, [0, 1], radius_angstrom=2.0, force_constant_hartree_per_angstrom2=0.0)
    assert zero_energy == 0.0
    np.testing.assert_array_equal(zero_gradient, np.zeros_like(xyz))
    boundary_x = 2.0 + masses[1] / masses[:2].sum()
    inside = xyz.copy(); inside[2, 0] = boundary_x - 0.01
    boundary = xyz.copy(); boundary[2, 0] = boundary_x
    outside = xyz.copy(); outside[2, 0] = boundary_x + 0.01
    e_inside, g_inside = outer_sphere_energy_gradient(inside, masses, [0, 1], **kwargs)
    e_boundary, g_boundary = outer_sphere_energy_gradient(boundary, masses, [0, 1], **kwargs)
    e_outside, g_outside = outer_sphere_energy_gradient(outside, masses, [0, 1], **kwargs)
    assert e_inside == e_boundary == 0.0
    assert np.linalg.norm(g_inside) == np.linalg.norm(g_boundary) == 0.0
    assert e_outside < 3e-5 and np.linalg.norm(g_outside) < 0.01


def test_outer_wall_applies_translation_force_to_solvent_fragment_com():
    from pyar.microsolvation.wall import outer_sphere_energy_gradient

    xyz = np.array([[0., 0., 0.], [1., 0., 0.], [4., 0., 0.], [4.5, 0.8, 0.]])
    masses = np.array([12., 1., 16., 1.])
    kwargs = dict(radius_angstrom=2.0, force_constant_hartree_per_angstrom2=0.5,
                  solvent_fragments=[[2, 3]])
    energy, gradient = outer_sphere_energy_gradient(xyz, masses, [0, 1], **kwargs)
    delta = 1e-6
    moved = xyz.copy(); moved[3, 1] += delta
    plus = outer_sphere_energy_gradient(moved, masses, [0, 1], **kwargs)[0]
    moved[3, 1] -= 2 * delta
    minus = outer_sphere_energy_gradient(moved, masses, [0, 1], **kwargs)[0]
    assert gradient[3, 1] == pytest.approx((plus - minus) / (2 * delta), rel=1e-6)
    # The wall translates the entire solvent fragment without applying torque.
    assert gradient[2, 0] / masses[2] == pytest.approx(gradient[3, 0] / masses[3])
    assert gradient[2, 1] / masses[2] == pytest.approx(gradient[3, 1] / masses[3])
    assert np.isfinite(energy)


def test_outer_wall_gradient_matches_central_differences_for_solute_and_solvent_atoms():
    from pyar.microsolvation.wall import outer_sphere_energy_gradient

    xyz = np.array([
        [0.0, 0.0, 0.0], [1.1, 0.1, 0.0], [-0.4, 0.9, 0.2],
        [4.0, 0.2, 0.1], [4.4, 0.8, -0.1], [3.8, -0.7, 0.2],
        [-3.6, 0.1, 0.3], [-4.1, 0.7, 0.1], [-3.2, -0.5, -0.2],
    ])
    masses = np.array([12.0, 1.0, 1.0, 16.0, 1.0, 1.0, 16.0, 1.0, 1.0])
    kwargs = dict(radius_angstrom=2.2, force_constant_hartree_per_angstrom2=0.37,
                  solvent_fragments=[[3, 4, 5], [6, 7, 8]])
    _, analytic = outer_sphere_energy_gradient(xyz, masses, [0, 1, 2], **kwargs)
    numerical = np.zeros_like(xyz)
    delta = 1e-6
    for atom in range(len(xyz)):
        for axis in range(3):
            plus = xyz.copy(); plus[atom, axis] += delta
            minus = xyz.copy(); minus[atom, axis] -= delta
            e_plus = outer_sphere_energy_gradient(plus, masses, [0, 1, 2], **kwargs)[0]
            e_minus = outer_sphere_energy_gradient(minus, masses, [0, 1, 2], **kwargs)[0]
            numerical[atom, axis] = (e_plus - e_minus) / (2 * delta)
    np.testing.assert_allclose(analytic, numerical, rtol=2e-6, atol=2e-8)


def test_outer_wall_rejects_nonfinite_coordinates_and_bad_fragment_partition():
    from pyar.microsolvation.wall import outer_sphere_energy_gradient

    xyz = np.array([[0., 0., 0.], [2., 0., 0.]])
    masses = np.array([12., 16.])
    with pytest.raises(ValueError, match="finite"):
        outer_sphere_energy_gradient([[0., 0., 0.], [float("nan"), 0., 0.]], masses, [0],
                                     radius_angstrom=1., force_constant_hartree_per_angstrom2=1.)
    with pytest.raises(ValueError, match="partition"):
        outer_sphere_energy_gradient(xyz, masses, [0], radius_angstrom=1.,
                                     force_constant_hartree_per_angstrom2=1., solvent_fragments=[[0]])


def test_wall_request_fingerprint_rejects_changed_parameters(tmp_path):
    from pyar.microsolvation.layout_v2 import MicrosolvationRunStateV2
    from pyar.state.microsolvation import MicrosolvationStateError

    qc = {"software": "xtb", "geometry_optimizer": "geometric"}
    first = MicrosolvationRequest(methane(), water(), 1, backend_parameters=qc,
                                  wall="outer-sphere", wall_radius=4.0, wall_force_constant=0.5)
    changed = MicrosolvationRequest(methane(), water(), 1, backend_parameters=qc,
                                    wall="outer-sphere", wall_radius=3.0, wall_force_constant=0.5)
    root = tmp_path / "wall-run"
    initial = first.solute.copy()
    initial.solute_atom_indices = tuple(range(len(initial)))
    initial.solvent_fragments = ()
    initial.microsolvation_coverage = {}
    MicrosolvationRunStateV2.create(root, first.to_state_dict(), initial, first.solvent)
    with pytest.raises(MicrosolvationStateError, match="differs"):
        MicrosolvationRunStateV2.load(root, changed.to_state_dict())


def test_geometry_only_workflow_is_bounded_and_records_coverage_and_identity(tmp_path):
    req = request(count=2, orientations=3, maximum=2)
    result = microsolvate(req, output=tmp_path / "microsolvation")
    assert result.status == "completed"
    assert len(result.selected_paths) <= 2
    assert result.metadata["completed_solvents"] == 2
    assert len(result.metadata["coverage_history"]) == 2
    assert (tmp_path / "microsolvation" / "steps" / "step_001" / "coverage.json").is_file()
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
