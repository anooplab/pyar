"""Real constrained optimizer and scan continuation boundary regressions."""

import json
from pathlib import Path

import numpy as np
import pytest

from pyar.backends.bond_scan import BondScanRequest, load_bond_scan_result, run_bond_scan
from pyar.core.molecule import Molecule
from pyar.energy_gradient_providers import EnergyGradientResult


def test_real_geometric_scan_constrains_distance_and_relaxes_other_bonds(tmp_path, monkeypatch):
    """A nonzero constrained force must not be mistaken for nonconvergence."""
    from pyar import energy_gradient_providers as providers
    caller = tmp_path / "caller"
    caller.mkdir()
    sentinel = caller / "pyar_geometric_state.json"
    sentinel.write_text("caller state")
    monkeypatch.chdir(caller)

    class PairPotential:
        def evaluate(self, molecule, coordinates_bohr):
            xyz = np.asarray(coordinates_bohr)
            gradient = np.zeros_like(xyz)
            energy = 0.
            for i, j in ((0, 1), (0, 2), (1, 2)):
                delta = xyz[i] - xyz[j]
                distance = np.linalg.norm(delta)
                displacement = distance - 2.5
                energy += 0.5 * displacement**2
                force = displacement * delta / distance
                gradient[i] += force
                gradient[j] -= force
            return EnergyGradientResult(energy, gradient)

    monkeypatch.setitem(providers.ENERGY_GRADIENT_PROVIDERS, "xtb", lambda params: PairPotential())
    molecule = Molecule(["H", "H", "H"], [[0., 0., 0.], [1.6, 0., 0.], [0.8, 1.1, 0.]])
    molecule.fragments = [[0], [1, 2]]
    request = BondScanRequest(0, 1, 1.6, 1.3, 3)
    result = run_bond_scan(molecule, request, tmp_path, {"software": "xtb", "opt_cycles": 100})
    assert result.success, result.status
    assert Path.cwd() == caller
    assert sentinel.read_text() == "caller state"
    assert (tmp_path / "pyar_geometric_state.json").is_file()
    for (_, xyz, _), target in zip(result.frames, np.linspace(1.6, 1.3, 3)):
        assert np.linalg.norm(xyz[0] - xyz[1]) == pytest.approx(target, abs=1e-3)
        assert np.linalg.norm(xyz[0] - xyz[2]) == pytest.approx(2.5 * 0.529177210903, abs=1e-3)
    assert load_bond_scan_result(molecule, request, tmp_path).success
    original_charge = molecule.charge
    molecule.charge = original_charge + 1
    assert not load_bond_scan_result(molecule, request, tmp_path).success
    molecule.charge = original_charge
    assert not load_bond_scan_result(molecule, BondScanRequest(0, 1, 1.6, 1.2, 3), tmp_path).success
    result.trajectory_path.write_text(result.trajectory_path.read_text() + "tampering")
    assert not load_bond_scan_result(molecule, request, tmp_path).success


@pytest.fixture
def scan_continuation(tmp_path):
    from pyar.neb import _write_xyz_trajectory
    root = tmp_path / "scan"
    directory = root / "orientation_000"
    scan = directory / "scan"
    scan.mkdir(parents=True)
    frames = [np.array([[0., 0., 0.], [d, 0., 0.]]) for d in (2., 1.5, 1.)]
    _write_xyz_trajectory(directory / "start.xyz", ["H", "H"], frames[:1], [-1.])
    _write_xyz_trajectory(scan / "final_scan.xyz", ["H", "H"], frames[-1:], [-1.])
    _write_xyz_trajectory(scan / "scan_trajectory.xyz", ["H", "H"], frames, [-1., 0., -1.])
    profile = scan / "scan_profile.json"
    profile.write_text(json.dumps({"points": [{"energy_hartree": energy} for energy in (-1., 0., -1.)]}))
    (root / "request.json").write_text(json.dumps({"orientation_definitions": [{
        "orientation": 0, "molecule": {"atoms": ["H", "H"], "charge": -1, "multiplicity": 2}}]}))
    result = {"orientation": 0, "scan_status": "success", "scan_profile_json": str(profile),
              "trajectory_path": str(scan / "scan_trajectory.xyz"), "final_scan_path": str(scan / "final_scan.xyz")}
    return root, result


@pytest.mark.parametrize("through,expected", [
    ("neb", ["relax", "neb"]), ("ts", ["relax", "neb", "ts"]),
    ("frequency", ["relax", "neb", "ts", "frequency"]),
    ("irc", ["relax", "neb", "ts", "frequency", "irc"]),
    ("endpoints", ["relax", "neb", "ts", "frequency", "irc", "endpoint-relax"]),
    ("endpoint-frequency", ["relax", "neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency"]),
    ("all", ["relax", "neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency"]),
])
@pytest.mark.parametrize("software", ["xtb", "orca", "gaussian", "aimnet_2"])
def test_continuation_stage_order_and_backend_identity(scan_continuation, monkeypatch, through, expected, software):
    from pyar.workflows.scan_path import continue_scan_paths
    from pyar.neb import _STAGE_GATES
    root, result = scan_continuation
    calls = []
    def run(**kwargs):
        calls.append(kwargs)
        return {**{key: True for key in _STAGE_GATES[kwargs["stage"]]}, "interior_maximum": True}
    monkeypatch.setattr("pyar.neb.run_neb", run)
    continue_scan_paths(root, [result], {"software": software, "xtb_model": "gfn2", "nprocs": 3},
                        through, {"ts_optimizer": "sella", "sella_internal_coordinates": True})
    assert [call["stage"] for call in calls] == expected
    assert all((call["software"], call["charge"], call["multiplicity"], call["nprocs"])
               == (software, -1, 2, 3) for call in calls)
    assert all(call["reuse"] for call in calls)
    assert calls[1]["ts_guess"].name == "scan_waypoint.xyz"
    assert all(call["ts_guess"] is None for call in calls if call["stage"] == "ts")
    assert result["continuation_status"] == "complete"
    assert result["reactant_product_connection_confirmed"] == (through in {"all", "endpoint-frequency"})


@pytest.mark.parametrize("failure", ["neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency"])
def test_continuation_stops_at_failed_scientific_gate(scan_continuation, monkeypatch, failure):
    from pyar.workflows.scan_path import continue_scan_paths
    from pyar.neb import _STAGE_GATES
    root, result = scan_continuation
    calls = []
    def run(**kwargs):
        stage = kwargs["stage"]
        calls.append(stage)
        return {**{key: stage != failure for key in _STAGE_GATES[stage]}, "interior_maximum": True}
    monkeypatch.setattr("pyar.neb.run_neb", run)
    continue_scan_paths(root, [result], {"software": "xtb"}, "all", {})
    assert calls[-1] == failure
    assert result["continuation_status"] == "scientific_gate_failed"
    assert not result["reactant_product_connection_confirmed"]


def test_cli_xtb_full_cycle_does_not_require_basis(tmp_path, monkeypatch):
    from pyar.scripts import scan_bond
    path = tmp_path / "H.xyz"
    path.write_text("1\nH\nH 0 0 0\n")
    calls = []
    monkeypatch.setattr(scan_bond, "run_scan_bond", lambda *args: calls.append(args))
    scan_bond.main([str(path), str(path), "--atoms", "0", "0", "-N", "1", "--software", "xtb",
                    "--xtb-model", "gfn2", "--through", "all", "--ts-optimizer", "sella"])
    assert calls[0][4]["xtb_model"] == "gfn2"
    assert calls[0][9] == "all"
    assert calls[0][10]["ts_optimizer"] == "sella"


def test_orca_builtin_xtb_gradient_omits_dft_keywords():
    from ase import Atoms
    from pyar.energy_gradient_providers import OrcaEnergyGradientProvider
    keywords = OrcaEnergyGradientProvider({"method": "XTB2", "basis": None})._build_keyword(Atoms("H2"))
    assert keywords.startswith("! GFN2-xTB ENGRAD")
    assert all(word not in keywords for word in ("def2", "RI ", "D3BJ", "KDIIS"))


def test_scan_aborts_at_unconverged_point_without_final_artifact(tmp_path, monkeypatch):
    from geometric.optimize import OPT_STATE
    captured = []
    class FailingOptimizer:
        state = OPT_STATE.FAILED
        def __init__(self, coordinates, molecule, internal, engine, scratch, params, **kwargs):
            captured.append((internal, params))
            self.progress = SimpleNamespace(xyzs=[molecule.xyzs[0]], qm_energies=[-1.])
        def optimizeGeometry(self):
            return self.progress
    from types import SimpleNamespace
    monkeypatch.setattr('geometric.optimize.Optimizer', FailingOptimizer)
    molecule = Molecule(['H', 'H', 'H'], [[0., 0., 0.], [1.6, 0., 0.], [0.8, 1.1, 0.]])
    result = run_bond_scan(molecule, BondScanRequest(0, 1, 1.6, 1.3, 3), tmp_path, {'software': 'xtb'})
    assert not result.success
    assert result.status == 'point_not_converged'
    assert len(captured) == 1
    assert captured[0][0].conmethod == 1
    assert not (tmp_path / 'final_scan.xyz').exists()


def test_corrupt_scan_orientation_does_not_abort_other_orientations(scan_continuation, monkeypatch):
    from pyar.workflows.scan_path import continue_scan_paths
    from pyar.neb import _STAGE_GATES
    root, result = scan_continuation
    corrupt = dict(result, trajectory_path=str(root / 'missing.xyz'))
    monkeypatch.setattr('pyar.neb.run_neb', lambda **kwargs: {
        **{key: True for key in _STAGE_GATES[kwargs['stage']]}, 'interior_maximum': True})
    continue_scan_paths(root, [corrupt, result], {'software': 'xtb'}, 'neb', {})
    assert corrupt['continuation_status'] == 'failed'
    assert corrupt['failed_stage'] == 'initialization'
    assert result['continuation_status'] == 'complete'


@pytest.mark.parametrize('options', [{'images': 4}, {'max_cycles': 0}, {'ts_fmax': float('nan')},
                                    {'climb': -1}, {'interpolation': 'unknown'}])
def test_invalid_continuation_settings_rejected_before_scanning(options):
    from pyar.workflows.scan_path import validate_continuation
    with pytest.raises(ValueError):
        validate_continuation('all', options)


def test_orca_open_shell_singlet_keeps_unrestricted_gradient_method():
    from ase import Atoms
    from pyar.energy_gradient_providers import OrcaEnergyGradientProvider
    provider = OrcaEnergyGradientProvider({'method': 'BP86', 'basis': 'def2-SVP',
                                          'multiplicity': 1, 'scftype': 'uhf'})
    assert ' UKS' in provider._build_keyword(Atoms('H2'))


def test_cli_reports_scientific_failure_with_nonzero_exit(tmp_path, monkeypatch, capsys):
    from pyar.scripts import scan_bond
    path = tmp_path / 'H.xyz'
    path.write_text('1\nH\nH 0 0 0\n')
    monkeypatch.setattr(scan_bond, 'run_scan_bond', lambda *args: {
        'status': 'failed', 'output_dir': str(tmp_path),
        'results': [{'continuation_status': 'scientific_gate_failed', 'failed_stage': 'frequency'}]})
    with pytest.raises(SystemExit) as error:
        scan_bond.main([str(path), str(path), '--atoms', '0', '0', '-N', '1', '--software', 'xtb', '--through', 'all'])
    assert error.value.code == 1
    assert 'failed through all' in capsys.readouterr().out


def test_scan_restores_caller_directory_after_failure(tmp_path, monkeypatch):
    from pyar.backends import bond_scan
    caller = Path.cwd()
    def fail(*args):
        assert Path.cwd() == tmp_path
        raise RuntimeError('fixture failure')
    monkeypatch.setattr(bond_scan, '_run_bond_scan', fail)
    with pytest.raises(RuntimeError, match='fixture failure'):
        bond_scan.run_bond_scan(None, None, tmp_path, {})
    assert Path.cwd() == caller
