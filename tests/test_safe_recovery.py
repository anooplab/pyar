"""SAFE recovery invariants, synthetic native evidence and interruption boundaries."""
from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pyar.persistence import read_json, write_json, canonical_json_bytes
from pyar.recovery import controller
from pyar.recovery.models import (OptimizationRequest, ExecutionObservation, ProcessRecord, Outcome)
from pyar.recovery.observations import snapshot_artifacts
from pyar.recovery.policy import RecoveryPolicy
from pyar.optimization_service import inspect_optimization


@pytest.fixture
def molecule(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return SimpleNamespace(name='h2', atoms_list=['H', 'H'], coordinates=[[0.,0.,0.],[0.,0.,1.]],
        charge=0, multiplicity=1, scftype='rhf', fragments=[])


@pytest.fixture
def settings():
    return dict(software='xtb', xtb_model='gfn2', geometry_optimizer='geometric',
                opt_target='minimum', opt_cycles=3, opt_threshold='normal', nprocs=1)


def xyz(coordinates, energy=-1.):
    return '2\nenergy: ' + str(energy) + '\n' + ''.join(
        f'H {x} {y} {z}\n' for x,y,z in coordinates)


def produce(directory, *, converged=False, progress=True, unknown=False):
    value = read_json(directory / 'input.json')
    req = controller._request(value['request'])
    root = directory / ('job_' + req.name)
    root.mkdir()
    (directory / 'worker.out').write_text('synthetic native execution\n')
    (root / ('trial_' + req.name + '.xyz')).write_text(xyz(req.coordinates_angstrom))
    start = req.coordinates_angstrom
    end = [list(row) for row in start]
    end[-1][-1] -= .02 if progress else 0
    gradient = .00001 if converged else .005
    energy = -1.1 if progress else -1.
    manifest = []
    for index, coords, e, g in [(0, start, -1., .01), (1, end, energy, gradient if progress else .01)]:
        filename = f'ore_accepted_{index:06d}.json'
        write_json(root / filename, dict(schema_version=1, version='1.1.1', step=index,
            objective_fingerprint=req.objective_fingerprint, symbols=['H','H'],
            coordinates_angstrom=coords, energy_hartree=e,
            gradient_hartree_per_bohr=[[0,0,g],[0,0,-g]], gradient_rms=g, gradient_max=g))
        manifest.append(filename)
    write_json(root / 'ore_accepted.json', dict(schema_version=1, files=manifest))
    text = ('geomeTRIC started. Version: 1.1.1\n'
        f'Step 0 : Gradient = .01/.01 Energy = -1.0\n'
        f'Step 1 : Grad = {gradient}/{gradient} E (change) = {energy}\n')
    text += 'Converged!\n' if converged else 'Maximum iterations reached\n'
    if unknown:
        text = text.replace('1.1.1', '99.0.0')
    (root / 'geometric.out').write_text(text)
    if converged:
        forces = [[0,0,-gradient*27.211386245988/.529177210903],
                  [0,0,gradient*27.211386245988/.529177210903]]
        write_json(root / 'pyar_geometric_state.json', dict(schema_version=1, software='xtb',
            opt_target='minimum', positions_angstrom=end, backend_energy_hartree=energy,
            total_energy_hartree=energy, forces_ev_per_angstrom=forces))
        (root / ('trial_' + req.name + '_optim.xyz')).write_text(xyz(end, energy))
    write_json(root / 'ore_evaluations.json', {'used': 2})
    for number in (1, 2):
        (root / f'ore_backend_{number:06d}.out').write_text(
            'xtb version 6.7.1\nHamiltonian GFN2-xTB\nnormal termination of xtb\n')
    execution = ExecutionObservation(req.request_fingerprint, snapshot_artifacts(root),
        (ProcessRecord('python', 0 if converged else 1),), True,
        True if converged else 'CycleExceeded', elapsed_seconds=.1)
    result = inspect_optimization(req, root, execution=execution)
    write_json(directory / 'result.json', result.to_dict())
    return 0


def fake_worker(monkeypatch, *, always_limit=False, progress=True, unknown=False):
    calls = []
    def run(directory, timeout):
        calls.append(directory)
        return produce(directory, converged=len(calls)>1 and not always_limit,
                       progress=progress, unknown=unknown)
    monkeypatch.setattr(controller, '_execute_attempt', run)
    return calls


def test_failure_then_success_and_replay(molecule, settings, monkeypatch):
    calls = fake_worker(monkeypatch)
    policy = RecoveryPolicy('safe')
    report = controller.run_recovery(molecule, settings, policy)
    assert report.result.status is Outcome.RECOVERED_CONVERGED
    assert report.reserved_cycles == 9 and report.backend_evaluations == 4
    assert [a.action for a in report.result.attempts] == ['initial', 'cycle_extension']
    first, second = report.result.attempts
    assert first.request.objective_fingerprint == second.request.objective_fingerprint
    assert first.original_scientific_fingerprint == second.original_scientific_fingerprint
    assert first.request.request_fingerprint != second.request.request_fingerprint
    assert second.request.coordinates_angstrom != first.request.coordinates_angstrom
    assert molecule.coordinates[-1][-1] == 1.
    # Derived outcome can be absent or corrupt; committed native evidence wins.
    outcome = Path('job_h2/recovery/outcome.json')
    outcome.write_text('{bad json')
    before = snapshot_artifacts(Path('job_h2/recovery/attempts'))
    replay = controller.run_recovery(molecule, settings, policy, resume=True)
    assert replay.result.converged and len(calls) == 2
    assert snapshot_artifacts(Path('job_h2/recovery/attempts')) == before


def test_ordinary_limit_never_completed(molecule, settings, monkeypatch):
    calls = fake_worker(monkeypatch, always_limit=True)
    report = controller.run_recovery(molecule, settings, RecoveryPolicy('safe', max_attempts=2))
    assert len(calls) == 2 and not report.result.converged
    assert report.result.status is Outcome.INCOMPLETE_USABLE
    assert not Path('job_h2/result_h2.xyz').exists()


@pytest.mark.parametrize('problem', ['no_progress', 'unknown'])
def test_no_blind_retry(molecule, settings, monkeypatch, problem):
    calls = fake_worker(monkeypatch, progress=problem!='no_progress', unknown=problem=='unknown')
    report = controller.run_recovery(molecule, settings, RecoveryPolicy('safe'))
    assert not report.result.converged and len(calls) == 1
    assert report.stop_reason == 'no_verified_safe_action'


@pytest.mark.parametrize('policy,expected', [
    (RecoveryPolicy('safe', max_cycles=3), 1),
    (RecoveryPolicy('safe', max_evaluations=2), 1),
    (RecoveryPolicy('safe', max_attempts=1), 1),
])
def test_finite_budgets(molecule, settings, monkeypatch, policy, expected):
    calls = fake_worker(monkeypatch, always_limit=True)
    report = controller.run_recovery(molecule, settings, policy)
    assert len(calls) == expected and not report.result.converged
    assert report.reserved_cycles <= policy.max_cycles
    assert report.backend_evaluations <= policy.max_evaluations


def test_crash_after_receipt_before_commit(molecule, settings, monkeypatch):
    calls = fake_worker(monkeypatch)
    real_write = controller.write_json
    def crash(path, value, **kwargs):
        if Path(path).name == 'commit.json':
            raise OSError('simulated power loss before commit')
        return real_write(path, value, **kwargs)
    monkeypatch.setattr(controller, 'write_json', crash)
    with pytest.raises(OSError, match='power loss'):
        controller.run_recovery(molecule, settings, RecoveryPolicy('safe'))
    assert len(calls) == 1
    monkeypatch.setattr(controller, 'write_json', real_write)
    report = controller.run_recovery(molecule, settings, RecoveryPolicy('safe'), resume=True)
    assert report.result.converged and len(calls) == 2  # first execution adopted, not repeated


def test_crash_after_commit_before_outcome(molecule, settings, monkeypatch):
    calls = fake_worker(monkeypatch)
    real_write = controller.write_json
    def crash(path, value, **kwargs):
        if Path(path).name == 'outcome.json':
            raise OSError('simulated power loss after commit')
        return real_write(path, value, **kwargs)
    monkeypatch.setattr(controller, 'write_json', crash)
    with pytest.raises(OSError):
        controller.run_recovery(molecule, settings, RecoveryPolicy('safe'))
    monkeypatch.setattr(controller, 'write_json', real_write)
    report = controller.run_recovery(molecule, settings, RecoveryPolicy('safe'), resume=True)
    assert report.result.converged and len(calls) == 2


def test_partial_worker_never_reissued(molecule, settings, monkeypatch):
    calls = []
    def fail(directory, timeout):
        calls.append(directory)
        (directory / 'worker.out').write_text('partial')
        return 1
    monkeypatch.setattr(controller, '_execute_attempt', fail)
    controller.run_recovery(molecule, settings, RecoveryPolicy('safe'))
    report = controller.run_recovery(molecule, settings, RecoveryPolicy('safe'), resume=True)
    assert len(calls) == 1 and not report.result.converged
    assert report.stop_reason == 'interrupted_attempt_without_receipt'


@pytest.mark.parametrize('problem', ['native_file', 'checkpoint', 'request', 'intent', 'symlink'])
def test_corrupt_or_incompatible_resume(molecule, settings, monkeypatch, problem):
    fake_worker(monkeypatch)
    controller.run_recovery(molecule, settings, RecoveryPolicy('safe'))
    first = Path('job_h2/recovery/attempts/attempt_001')
    if problem == 'request':
        settings = dict(settings, opt_threshold='tight')
    elif problem == 'intent':
        (first / 'input.json').write_text('{}')
    elif problem == 'symlink':
        path = first / 'job_h2/geometric.out'
        path.unlink()
        path.symlink_to('/etc/passwd')
    else:
        path = first / 'job_h2' / ('geometric.out' if problem=='native_file' else 'ore_accepted_000001.json')
        path.write_text('{}')
    with pytest.raises(ValueError):
        controller.run_recovery(molecule, settings, RecoveryPolicy('safe'), resume=True)


@pytest.mark.parametrize('changed', [dict(geometry_optimizer='native'), dict(software='orca'),
    dict(multiplicity=3), dict(gamma=10), dict(constraints=[1]), dict(custom_keywords='foo'),
    dict(xtb_model='gxtb'), dict(microsolvation_wall={'x':1}), dict(unknown_control=True)])
def test_unsupported_preflight_no_writes(molecule, settings, changed):
    if 'multiplicity' in changed:
        molecule.multiplicity = changed['multiplicity']
    before = set(Path('.').iterdir())
    with pytest.raises(ValueError):
        controller.preflight_recovery([molecule], dict(settings, **changed), RecoveryPolicy('safe'))
    assert set(Path('.').iterdir()) == before


def test_stale_cache_not_adopted(molecule, settings):
    Path('job_h2').mkdir()
    Path('job_h2/result_h2.xyz').write_text(xyz(molecule.coordinates))
    before = snapshot_artifacts(Path('job_h2'))
    with pytest.raises(ValueError, match='Existing job'):
        controller.run_recovery(molecule, settings, RecoveryPolicy('safe'))
    assert snapshot_artifacts(Path('job_h2')) == before


@pytest.mark.parametrize('values', [dict(max_attempts=0), dict(max_attempts=9), dict(max_cycles=-1),
    dict(max_evaluations=10001), dict(wall_seconds=float('nan')), dict(wall_seconds=0)])
def test_invalid_policy(values):
    with pytest.raises(ValueError):
        RecoveryPolicy('safe', **values)


def test_check_read_only_and_legacy_default(molecule, settings, monkeypatch):
    from pyar.scripts import optimize
    from pyar import optimiser
    Path('h2.xyz').write_text(xyz(molecule.coordinates))
    monkeypatch.setattr(optimize, 'preflight', lambda *a, **k: ['fake backend'])
    before = snapshot_artifacts(Path('.'))
    optimize.main(['h2.xyz', '--backend', 'xtb', '--geometry-optimizer', 'geometric',
                   '--recovery', 'safe', '--check'])
    assert snapshot_artifacts(Path('.')) == before
    called = []
    monkeypatch.setattr(optimiser, 'optimise', lambda *a: called.append(a) or True)
    monkeypatch.setattr(controller, 'run_recovery', lambda *a, **k: pytest.fail('recovery default changed'))
    optimize.main(['h2.xyz', '--backend', 'xtb'])
    assert len(called) == 1


def test_off_api_runs_exactly_once_without_recovery_directory(molecule, settings, monkeypatch):
    import pyar.optimization_service as service
    calls = []
    expected = controller.OptimizationResult(Outcome.INVALID_RESULT, (), (), ())
    monkeypatch.setattr(service, 'optimize_once', lambda *a: calls.append(a) or expected)
    report = controller.run_recovery(molecule, settings)
    assert report.result is expected and len(calls) == 1
    assert not Path('job_h2').exists()


def test_walltime_kills_backend_process_group(tmp_path):
    import os
    import subprocess
    import sys
    import time
    if os.name != 'posix':
        pytest.skip('Process group test requires POSIX')
    from unittest.mock import patch
    real_popen = subprocess.Popen
    script = "import subprocess,sys,time; subprocess.Popen([sys.executable,'-c',\"import time,pathlib; time.sleep(.5); pathlib.Path('leaked-child').write_text('bad')\"]); time.sleep(10)"
    def spawn(command, **kwargs):
        return real_popen([sys.executable, '-c', script], **kwargs)
    with patch.object(controller.subprocess, 'Popen', spawn):
        with pytest.raises(subprocess.TimeoutExpired):
            controller._execute_attempt(tmp_path, .2)
    time.sleep(.6)
    assert not (tmp_path / 'leaked-child').exists()


def test_fake_executable_end_to_end(molecule, settings, monkeypatch, tmp_path):
    """Real worker/driver/optimizer, synthetic harmonic backend; no chemistry executable."""
    import os
    import sys
    geometric = pytest.importorskip('geometric')
    if geometric.__version__ != '1.1.1':
        pytest.skip('Accepted-step pilot deliberately pins geomeTRIC 1.1.1')
    binary = tmp_path / 'bin'
    binary.mkdir()
    fake = binary / 'xtb'
    fake.write_text(f'#!{sys.executable}\n' + '''import math, pathlib, sys
rows = pathlib.Path(sys.argv[1]).read_text().splitlines()[2:]
xyz = [[float(x) / .529177210903 for x in row.split()[1:]] for row in rows]
v = [b-a for a,b in zip(*xyz)]
r = math.sqrt(sum(x*x for x in v))
d = r-1.4
g = [.1*d*x/r for x in v]
pathlib.Path('gradient').write_text('\\n'.join(' '.join(str(sign*x) for x in g) for sign in (-1,1)))
print('xtb version 6.7.1')
print('Hamiltonian GFN2-xTB')
print('| TOTAL ENERGY', -.99 + .05*d*d, 'Eh |')
print('normal termination of xtb')
''')
    fake.chmod(0o755)
    fake_source = fake.read_text()
    monkeypatch.setenv('PATH', str(binary)+os.pathsep+os.environ['PATH'])
    monkeypatch.setenv('PYTHONPATH', str(Path(__file__).resolve().parents[1]))
    report = controller.run_recovery(molecule, dict(settings, opt_cycles=2),
        RecoveryPolicy('safe', max_attempts=4, max_cycles=40, max_evaluations=50, wall_seconds=30))
    assert report.result.status is Outcome.RECOVERED_CONVERGED, report.to_dict()
    assert len(report.result.attempts) >= 2
    assert report.backend_evaluations <= 50
    # This actual worker must also replay without executing the fake backend.
    fake.write_text(f'#!{sys.executable}\nraise SystemExit(99)\n')
    replay = controller.run_recovery(molecule, dict(settings, opt_cycles=2), report.policy, resume=True)
    assert replay.result.converged
    # Put the synthetic backend back and exercise a real physical hard cap.
    # Two evaluations permit initial+one trial, never a third provider call.
    fake.write_text(fake_source)
    molecule.name = 'limited'
    limited = controller.run_recovery(molecule, settings,
        RecoveryPolicy('safe', max_evaluations=2, wall_seconds=30))
    assert not limited.result.converged and limited.backend_evaluations == 2
    assert limited.stop_reason == 'evaluation_budget_exhausted'
    backend = Path('job_limited/recovery/attempts/attempt_001/job_limited')
    assert len(list(backend.glob('ore_backend_*.out'))) == 2
    assert not Path('job_limited/result_limited.xyz').exists()


@pytest.mark.parametrize('problem', ['bad_residual', 'wrong_symbols', 'nan_energy', 'bad_initial',
    'path_traversal', 'wrong_objective', 'unknown_xtb', 'scc_failure', 'empty_manifest', 'negative_count'])
def test_fresh_invalid_checkpoint_not_authority(molecule, settings, tmp_path, problem):
    original = controller.prepare_request(molecule, settings, RecoveryPolicy('safe'))
    current = replace(original, settings_json=canonical_json_bytes(dict(original.settings, _ore_evaluation_limit=10)).decode())
    directory = tmp_path / 'attempt'
    directory.mkdir()
    write_json(directory / 'input.json', {'request': controller._encode(current)})
    produce(directory)
    root = directory / 'job_h2'
    last = root / 'ore_accepted_000001.json'
    value = read_json(last)
    if problem == 'bad_residual':
        value['gradient_rms'] = 0.
    elif problem == 'wrong_symbols':
        value['symbols'] = ['He', 'He']
    elif problem == 'nan_energy':
        last.write_text(last.read_text().replace('"energy_hartree":-1.1', '"energy_hartree":NaN'))
    elif problem == 'bad_initial':
        initial = read_json(root / 'ore_accepted_000000.json')
        initial['coordinates_angstrom'][-1][-1] = 7.
        write_json(root / 'ore_accepted_000000.json', initial)
    elif problem == 'path_traversal':
        write_json(root / 'ore_accepted.json', {'schema_version':1, 'files':['../../foreign.json']})
    elif problem == 'wrong_objective':
        value['objective_fingerprint'] = '0'*64
    elif problem in {'unknown_xtb', 'scc_failure'}:
        (root / 'ore_backend_000001.out').write_text('xtb version 99.0.0\n' if problem=='unknown_xtb'
            else 'xtb version 6.7.1\nHamiltonian GFN2-xTB\nSCC failed\nnormal termination of xtb\n')
    elif problem == 'empty_manifest':
        write_json(root / 'ore_accepted.json', {'schema_version':1, 'files':[]})
    elif problem == 'negative_count':
        write_json(root / 'ore_evaluations.json', {'used':-1})
    if problem in {'bad_residual', 'wrong_symbols', 'wrong_objective'}:
        write_json(last, value)
    execution = ExecutionObservation(current.request_fingerprint, snapshot_artifacts(root),
        (ProcessRecord('python', 1),), True, 'CycleExceeded')
    result = inspect_optimization(current, root, execution=execution)
    assert result.status is Outcome.INVALID_RESULT
    assert result.accepted_frame is None


def test_rejected_or_undecided_last_trial_not_adopted(molecule, settings, monkeypatch):
    calls = fake_worker(monkeypatch, always_limit=True)
    report = controller.run_recovery(molecule, settings, RecoveryPolicy('safe', max_attempts=1))
    attempt = report.result.attempts[0]
    assert attempt.steps[-1].index < attempt.request.settings['opt_cycles']
    assert report.result.accepted_frame.coordinates_angstrom[-1][-1] == .98


def test_check_resume_rejects_corrupt_native_artifact_without_writes(molecule, settings, monkeypatch):
    fake_worker(monkeypatch)
    policy = RecoveryPolicy('safe')
    controller.run_recovery(molecule, settings, policy)
    Path('job_h2/recovery/attempts/attempt_001/job_h2/geometric.out').write_text('corrupt')
    before = snapshot_artifacts(Path('job_h2'))
    with pytest.raises(ValueError):
        controller.preflight_recovery([molecule], settings, policy, resume=True)
    assert snapshot_artifacts(Path('job_h2')) == before


def test_physical_provider_is_not_called_after_hard_cap(tmp_path, monkeypatch):
    pytest.importorskip('ase')
    import numpy as np
    from ase import Atoms
    from pyar.backends.geometric import PyarGeometricCalculator
    from pyar.energy_gradient_providers import EnergyGradientResult
    monkeypatch.chdir(tmp_path)
    calls = []
    provider = SimpleNamespace(evaluate=lambda atoms, coords: calls.append(coords) or
        EnergyGradientResult(-1., np.array([[0.,0.,.01],[0.,0.,-.01]])))
    monkeypatch.setattr('pyar.backends.geometric._resolve_backend_evaluator', lambda *a: provider)
    calc = PyarGeometricCalculator(dict(software='xtb', gamma=None, _ore_safe=True, _ore_evaluation_limit=2))
    atoms = Atoms('H2', positions=[[0,0,0],[0,0,1.]])
    atoms.calc = calc
    atoms.get_forces()
    atoms.positions[-1,-1] += .1
    atoms.get_forces()
    atoms.positions[-1,-1] += .1
    with pytest.raises(RuntimeError, match='evaluation budget exhausted'):
        atoms.get_forces()
    assert len(calls) == 2
    assert read_json('ore_evaluations.json')['used'] == 2
