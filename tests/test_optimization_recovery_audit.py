"""Regression coverage for backend optimization and recovery safety.

Fake executables exercise native subprocess/file boundaries without chemistry tools.
"""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from pyar import optimiser
from pyar.backend_errors import BackendExecutionError
from pyar.core.molecule import Molecule


@pytest.fixture
def molecule():
    return Molecule(['H', 'H'], [[0., 0., 0.], [0., 0., .74]], name='audit')


@pytest.fixture
def qc():
    return dict(software='orca', method='HF', basis='def2-SVP', nprocs=1,
                scf_cycles=100, opt_cycles=10, opt_threshold='normal', gamma=None)


@pytest.fixture
def fake_program(tmp_path, monkeypatch):
    """A real child process writes only specified files under the test cwd."""
    monkeypatch.chdir(tmp_path)
    executable = tmp_path / 'fake-calculator'
    executable.write_text(
        f'#!{sys.executable}\n'
        'import json, pathlib, sys\n'
        'case = json.loads(pathlib.Path("case.json").read_text())\n'
        'for name, content in case.get("files", {}).items():\n'
        '    pathlib.Path(name).write_text(content)\n'
        'print(case.get("stdout", ""))\n'
        'sys.exit(case.get("returncode", 0))\n'
    )
    executable.chmod(0o755)
    def configure(*, stdout='', files=None, returncode=0):
        Path('case.json').write_text(json.dumps(dict(stdout=stdout, files=files or {}, returncode=returncode)))
        return str(executable)
    return configure


def xyz(energy=-1., symbols=('H', 'H'), distance=.74):
    rows = [f'{symbol} 0 0 {index * distance}' for index, symbol in enumerate(symbols)]
    return f'{len(rows)}\nenergy:{energy}\n' + '\n'.join(rows) + '\n'


def orca_runner(monkeypatch, molecule, qc, executable):
    from pyar.backends import orca
    monkeypatch.setattr(orca, 'require_executable', lambda *args: executable)
    return orca.Orca(molecule, qc)


def test_orca_normal_termination_without_geometry_convergence_is_rejected(fake_program, monkeypatch, molecule, qc):
    executable = fake_program(stdout='FINAL SINGLE POINT ENERGY -1.0\n****ORCA TERMINATED NORMALLY****',
                              files={'trial_audit.xyz': xyz()})
    assert orca_runner(monkeypatch, molecule, qc, executable).optimize() is not True


@pytest.mark.parametrize('energy_line', ['', 'FINAL SINGLE POINT ENERGY nan\n'])
def test_orca_missing_or_nonfinite_energy_is_rejected(fake_program, monkeypatch, molecule, qc, energy_line):
    executable = fake_program(stdout=energy_line + 'THE OPTIMIZATION HAS CONVERGED\n****ORCA TERMINATED NORMALLY****',
                              files={'trial_audit.xyz': xyz()})
    assert orca_runner(monkeypatch, molecule, qc, executable).optimize() is not True


def test_orca_does_not_reuse_old_wrong_symbol_xyz(fake_program, monkeypatch, molecule, qc):
    executable = fake_program(stdout='THE OPTIMIZATION HAS CONVERGED\nFINAL SINGLE POINT ENERGY -1\n****ORCA TERMINATED NORMALLY****')
    runner = orca_runner(monkeypatch, molecule, qc, executable)
    Path('trial_audit.xyz').write_text(xyz(symbols=('He', 'He')))
    assert runner.optimize() is not True


def test_orca_rejects_fresh_geometry_with_wrong_symbols(fake_program, monkeypatch, molecule, qc):
    executable = fake_program(
        stdout='THE OPTIMIZATION HAS CONVERGED\nFINAL SINGLE POINT ENERGY -1\n'
               '****ORCA TERMINATED NORMALLY****',
        files={'trial_audit.xyz': xyz(symbols=('He', 'He'))},
    )
    assert orca_runner(monkeypatch, molecule, qc, executable).optimize() is False


def test_orca_nonzero_exit_does_not_promote_existing_frame(fake_program, monkeypatch, molecule, qc):
    executable = fake_program(returncode=7, files={'trial_audit.xyz': xyz()},
                              stdout='FINAL SINGLE POINT ENERGY -1\n****ORCA TERMINATED NORMALLY****')
    assert orca_runner(monkeypatch, molecule, qc, executable).optimize() is False
    assert Path('trial_audit.xyz').is_file()  # retained on disk, not a validated accepted frame


def test_orca_one_atom_has_matrix_coordinates(fake_program, monkeypatch, qc):
    molecule = Molecule(['He'], [[0., 0., 0.]], name='audit')
    executable = fake_program(stdout='THE OPTIMIZATION HAS CONVERGED\nFINAL SINGLE POINT ENERGY -2\n****ORCA TERMINATED NORMALLY****',
                              files={'trial_audit.xyz': xyz(symbols=('He',))})
    runner = orca_runner(monkeypatch, molecule, qc, executable)
    assert runner.optimize() is True
    assert runner.optimized_coordinates.shape == (1, 3)


@pytest.mark.parametrize('scc_marker', [False, True])
def test_xtb_stale_success_cannot_authorize_current_run(fake_program, monkeypatch, molecule, scc_marker):
    from pyar.backends import xtb
    executable = fake_program(stdout='SCC failed to converge')
    monkeypatch.setattr(xtb, 'require_executable', lambda *args: executable)
    runner = xtb.Xtb(molecule, {'opt_threshold': 'normal'})
    Path('.xtboptok').touch()
    Path('xtbopt.xyz').write_text(xyz())
    Path('energy').write_text('$energy\n 1 -1.0\n$end\n')
    if scc_marker:
        Path('.sccnotconverged').touch()
    assert runner.optimize() is not True


def test_xtb_scc_marker_without_success_returns_failure(fake_program, monkeypatch, molecule):
    from pyar.backends import xtb
    executable = fake_program(stdout='SCC failed', files={'.sccnotconverged': ''})
    monkeypatch.setattr(xtb, 'require_executable', lambda *args: executable)
    assert xtb.Xtb(molecule, {'opt_threshold': 'normal'}).optimize() == 'SCFFailed'


def test_xtb_zero_exit_unknown_option_is_rejected(fake_program, monkeypatch, molecule):
    from pyar.backends import xtb
    executable = fake_program(stdout="Unknown option '--gxtb' provided")
    monkeypatch.setattr(xtb, 'require_executable', lambda *args: executable)
    assert xtb.Xtb(molecule, {'opt_threshold': 'normal'}).optimize() is False


@pytest.mark.parametrize('problem', ['missing', 'wrong_method', 'wrong_symbols'])
def test_generic_cache_needs_valid_convergence_and_request_identity(tmp_path, monkeypatch, molecule, qc, problem):
    monkeypatch.chdir(tmp_path)
    directory = tmp_path / 'job_audit'
    directory.mkdir()
    content = xyz(symbols=('He', 'He') if problem == 'wrong_symbols' else ('H', 'H'))
    (directory / 'result_audit.xyz').write_text(content)
    if problem == 'wrong_method':
        (directory / 'request.json').write_text(json.dumps({'software': 'xtb', 'xtb_model': 'gxtb'}))
    # If cache validation rejects it, an unavailable calculator must fail closed.
    monkeypatch.setattr(optimiser, 'build_geometry', lambda *args: (_ for _ in ()).throw(BackendExecutionError('cache rejected')))
    with pytest.raises(BackendExecutionError):
        optimiser.optimise(molecule, qc)


def test_matching_cache_is_reused_without_running_backend(tmp_path, monkeypatch, molecule, qc):
    monkeypatch.chdir(tmp_path)
    provenance = {"backend": "orca", "method": "HF", "version": "test-1"}
    monkeypatch.setattr(optimiser, '_backend_cache_provenance', lambda _qc: provenance)
    coordinates = np.asarray([[0., 0., 0.], [0., 0., .75]])
    calls = []

    class Geometry:
        atoms_list = ['H', 'H']
        optimized_coordinates = coordinates
        energy = -1.1

        def optimize(self):
            calls.append('run')
            from pyar.backends import write_xyz
            write_xyz(self.atoms_list, self.optimized_coordinates,
                      'result_audit.xyz', energy=self.energy)
            return True

    monkeypatch.setattr(optimiser, 'build_geometry', lambda *_args: Geometry())
    assert optimiser.optimise(molecule, qc) is True
    same_request = Molecule(['H', 'H'], [[0., 0., 0.], [0., 0., .74]], name='audit')
    assert optimiser.optimise(same_request, qc) is True
    assert calls == ['run']
    assert same_request.energy == -1.1


def test_unknown_backend_version_disables_cache_reuse_and_write(tmp_path, monkeypatch, molecule, qc):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(optimiser, '_backend_cache_provenance', lambda _qc: None)
    calls = []

    class Geometry:
        atoms_list = ['H', 'H']
        optimized_coordinates = np.asarray([[0., 0., 0.], [0., 0., .75]])
        energy = -1.1

        def optimize(self):
            calls.append(1)
            from pyar.backends import write_xyz
            write_xyz(self.atoms_list, self.optimized_coordinates,
                      'result_audit.xyz', energy=self.energy)
            return True

    monkeypatch.setattr(optimiser, 'build_geometry', lambda *_args: Geometry())
    assert optimiser.optimise(molecule, qc) is True
    assert optimiser.optimise(molecule, qc) is True
    assert calls == [1, 1]
    assert not (tmp_path / 'job_audit' / 'optimization_cache.json').exists()


def test_backend_version_change_invalidates_cache(tmp_path, monkeypatch, molecule, qc):
    monkeypatch.chdir(tmp_path)
    version = {'value': 'xTB 6.6.1'}
    monkeypatch.setattr(optimiser, '_backend_cache_provenance',
                        lambda _qc: {'backend': 'xtb', 'version': version['value']})
    calls = []

    class Geometry:
        atoms_list = ['H', 'H']
        optimized_coordinates = np.asarray([[0., 0., 0.], [0., 0., .75]])
        energy = -1.1

        def optimize(self):
            calls.append(version['value'])
            from pyar.backends import write_xyz
            write_xyz(self.atoms_list, self.optimized_coordinates,
                      'result_audit.xyz', energy=self.energy)
            return True

    monkeypatch.setattr(optimiser, 'build_geometry', lambda *_args: Geometry())
    assert optimiser.optimise(molecule, dict(qc, software='xtb')) is True
    version['value'] = 'xTB 6.7.1'
    assert optimiser.optimise(molecule, dict(qc, software='xtb')) is True
    assert calls == ['xTB 6.6.1', 'xTB 6.7.1']
    record = json.loads((tmp_path / 'job_audit' / 'optimization_cache.json').read_text())
    assert record['backend_provenance']['version'] == 'xTB 6.7.1'


def test_cache_with_changed_method_is_not_reused(tmp_path, monkeypatch, molecule, qc):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'job_audit').mkdir()
    result = tmp_path / 'job_audit' / 'result_audit.xyz'
    result.write_text(xyz())
    cache = tmp_path / 'job_audit' / 'optimization_cache.json'
    cache.write_text(json.dumps({
        'schema_version': 1, 'request_sha256': 'old-method',
        'result_sha256': 'valid-looking', 'status': 'converged',
    }))
    monkeypatch.setattr(
        optimiser, 'build_geometry',
        lambda *_args: (_ for _ in ()).throw(BackendExecutionError('new method executed')),
    )
    with pytest.raises(BackendExecutionError, match='new method executed'):
        optimiser.optimise(molecule, dict(qc, method='PBE0'))
    assert list((tmp_path / 'job_audit' / 'stale_cache').glob('attempt_*/result_audit.xyz'))


def geometric_case(fake_program, monkeypatch, *, status='cycle_exceeded', returncode=1,
                   energy=-1., state_distance=.74, write_geometry=True, force=0.):
    from pyar.backends import geometric
    state = dict(backend_energy_hartree=energy, positions_angstrom=[[0., 0., 0.], [0., 0., state_distance]],
                 forces_ev_per_angstrom=[[force, 0., 0.], [0., 0., 0.]],
                 optimization_status=status)
    files = {'pyar_geometric_state.json': json.dumps(state)}
    if write_geometry:
        files['trial_audit_optim.xyz'] = xyz()
    executable = fake_program(stdout='Maximum iterations reached' if returncode else '', files=files, returncode=returncode)
    monkeypatch.setattr(geometric, '_find_geometric_executable', lambda: executable)
    return geometric


def test_geometric_iteration_limit_retains_provisional_geometry(fake_program, monkeypatch, molecule, qc):
    geometric = geometric_case(fake_program, monkeypatch)
    runner = geometric.Geometric(molecule, dict(qc, software='xtb'))
    status = runner.optimize()
    assert status == 'CycleExceeded'
    assert not optimiser.is_success(status) and optimiser.is_usable(status)
    np.testing.assert_allclose(runner.optimized_coordinates, molecule.coordinates)
    assert runner.energy == -1.
    assert Path('result_audit.xyz').is_file()


def test_geometric_provisional_cache_is_not_promoted(fake_program, monkeypatch, molecule, qc):
    geometric_case(fake_program, monkeypatch)
    # optimise changes cwd; make the child configuration visible in that native job.
    Path('job_audit').mkdir()
    Path('job_audit/case.json').write_bytes(Path('case.json').read_bytes())
    params = dict(qc, software='xtb', geometry_optimizer='geometric')
    assert optimiser.optimise(molecule, params) == 'CycleExceeded'
    assert optimiser.optimise(molecule, params) != True


@pytest.mark.parametrize('problem', ['nan_energy', 'frame_mismatch', 'no_convergence_evidence', 'fallback_xyz', 'nonfinite_force'])
def test_geometric_invalid_result_is_not_success(fake_program, monkeypatch, molecule, qc, problem):
    state_status = 'converged' if problem == 'nonfinite_force' else 'unknown'
    geometric = geometric_case(fake_program, monkeypatch, status=state_status, returncode=0,
                               energy=float('nan') if problem == 'nan_energy' else -1.,
                               state_distance=3. if problem == 'frame_mismatch' else .74,
                               write_geometry=problem != 'fallback_xyz',
                               force=float('inf') if problem == 'nonfinite_force' else 0.)
    if problem == 'fallback_xyz':
        Path('unrelated_stale.xyz').write_text(xyz())
    with pytest.raises(BackendExecutionError):
        geometric.Geometric(molecule, dict(qc, software='xtb')).optimize()


def test_geometric_release_candidate_is_not_reported_as_converged(fake_program, monkeypatch, molecule, qc):
    from pyar.backends import geometric
    geometric_case(fake_program, monkeypatch, status='release_candidate', returncode=0)
    runner = geometric.Geometric(molecule, dict(qc, software='xtb'))
    assert runner.optimize() == 'ReleaseCandidate'
    assert not optimiser.is_success('ReleaseCandidate')
    assert optimiser.is_usable('ReleaseCandidate')


def test_geometric_missing_final_state_is_explicit(fake_program, monkeypatch, molecule, qc):
    from pyar.backends import geometric
    executable = fake_program()
    monkeypatch.setattr(geometric, '_find_geometric_executable', lambda: executable)
    with pytest.raises(BackendExecutionError, match='no usable final state'):
        geometric.Geometric(molecule, dict(qc, software='xtb')).optimize()


def test_sella_exception_chain_is_preserved(monkeypatch):
    from pyar import neb
    failure = ArithmeticError('singular internal coordinates')
    class BrokenSella:
        def __init__(self, *args, **kwargs):
            raise failure
    monkeypatch.setattr(neb, '_sella_api', lambda: (BrokenSella, 'fixture'))
    with pytest.raises(RuntimeError, match='Sella TS optimization failed') as captured:
        neb._optimize_sella(['H', 'H'], [[0., 0., 0.], [0., 0., .74]], SimpleNamespace(), 2, .01)
    assert captured.value.__cause__ is failure


@pytest.mark.parametrize('field,value', [('energy_hartree', float('nan')), ('gradient_hartree_per_bohr', [[float('inf'), 0., 0.]])])
def test_provider_contract_rejects_nonfinite_values(field, value):
    from pyar.energy_gradient_providers import EnergyGradientResult
    values = dict(energy_hartree=-1., gradient_hartree_per_bohr=[[0., 0., 0.]])
    values[field] = value
    with pytest.raises(ValueError):
        EnergyGradientResult(**values)


def test_incompatible_adaptive_checkpoint_is_rejected(tmp_path, monkeypatch):
    from pyar.backends.geometric import PyarGeometricCalculator
    monkeypatch.chdir(tmp_path)
    Path('pyar_bias_controller_state.json').write_text(json.dumps({'configuration': {}}))
    with pytest.raises(ValueError, match='Incompatible bias-controller restart configuration'):
        PyarGeometricCalculator(dict(software='xtb', gamma=100., bias_controller='adaptive',
                                     bias_controller_restart=True), fragment_indices=[[0], [1]])


@pytest.mark.parametrize('tamper', ['method', 'artifact', 'incomplete_summary'])
def test_path_stage_reuse_rejects_wrong_method_and_invalid_artifact(tmp_path, tamper):
    from pyar.neb import _save_stage, _load_stage
    calculator = SimpleNamespace(qc_params={'software': 'xtb', 'xtb_model': 'gfn2', 'charge': 0})
    geometry = tmp_path / 'ts.xyz'
    geometry.write_text(xyz())
    _save_stage(tmp_path, 'ts', {'ts_optimization_converged': True}, calculator, [], [geometry])
    if tamper == 'method':
        calculator.qc_params['xtb_model'] = 'gxtb'
    elif tamper == 'artifact':
        geometry.write_text(xyz(distance=3.))
    else:
        (tmp_path / 'ts_summary.json').write_text(json.dumps({'stage': 'ts', 'status': 'running'}))
    with pytest.raises(ValueError):
        _load_stage(tmp_path, 'ts', calculator)


@pytest.mark.parametrize('had_optimized_coordinates', [True, False])
def test_backend_failure_does_not_destroy_original_input(
        tmp_path, monkeypatch, molecule, qc, had_optimized_coordinates):
    monkeypatch.chdir(tmp_path)
    coordinates = molecule.coordinates.copy()
    molecule.energy = -0.5
    if had_optimized_coordinates:
        molecule.optimized_coordinates = coordinates.copy()
    class Failure:
        def optimize(self):
            molecule.optimized_coordinates = np.full_like(coordinates, 99.)
            raise BackendExecutionError('SCF failed at trial step')
    monkeypatch.setattr(optimiser, 'build_geometry', lambda *args: Failure())
    with pytest.raises(BackendExecutionError, match='SCF failed'):
        optimiser.optimise(molecule, qc)
    np.testing.assert_array_equal(molecule.coordinates, coordinates)
    if had_optimized_coordinates:
        np.testing.assert_array_equal(molecule.optimized_coordinates, coordinates)
    else:
        assert not hasattr(molecule, 'optimized_coordinates')
    assert molecule.energy == -0.5


def test_truncated_cache_raises_backend_error(tmp_path, monkeypatch, molecule, qc):
    monkeypatch.chdir(tmp_path)
    Path('job_audit').mkdir()
    Path('job_audit/result_audit.xyz').write_text('2\nenergy:-1\nH 0 0 0\n')
    with pytest.raises(BackendExecutionError):
        optimiser.optimise(molecule, qc)
