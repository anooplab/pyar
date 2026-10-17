"""R1 contracts: diagnosis cannot promote stale or incomplete optimizer output."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pyar.optimization_service import inspect_optimization, optimize_once
from pyar.recovery.models import OptimizationRequest, Outcome, FailureCode, ExecutionObservation, ProcessRecord
from pyar.recovery.observations import snapshot_artifacts


@pytest.fixture
def molecule():
    return SimpleNamespace(name='h2', atoms_list=['H', 'H'], coordinates=[[0., 0., 0.], [0., 0., .74]],
                           charge=0, multiplicity=1, scftype='rhf', fragments=[])


@pytest.fixture
def opt_request(molecule):
    return OptimizationRequest.from_molecule(molecule, {'software': 'xtb', 'xtb_model': 'gfn2',
                                                       'opt_threshold': 'normal', 'nprocs': 1})


def xyz(energy=-1., symbols=('H', 'H'), coordinates=None):
    positions = coordinates or [[0., 0., i * .74] for i in range(len(symbols))]
    return f'{len(symbols)}\nenergy: {energy}\n' + ''.join(
        f'{s} {x} {y} {z}\n' for s, (x,y,z) in zip(symbols, positions))


def receipt(opt_request, root, *, returncode=0, legacy_status=True):
    return ExecutionObservation(opt_request.request_fingerprint, snapshot_artifacts(root),
                                (ProcessRecord(opt_request.backend, returncode),), True,
                                legacy_status=legacy_status)


def xtb_output(root, *, energy=-1., converged=True):
    (root / 'trial_h2.xyz').write_text(xyz())
    (root / 'xtb.out').write_text('xtb version 6.7.1\nHamiltonian GFN2-xTB\n' +
        ('*** GEOMETRY OPTIMIZATION CONVERGED AFTER 2 ITERATIONS ***\n' if converged else '') +
        f'| TOTAL ENERGY {energy} Eh |\n| GRADIENT NORM 0.00001 Eh/a0 |\nnormal termination of xtb\n')
    (root / 'xtbopt.xyz').write_text(xyz(energy=energy))


def codes(result):
    return {d.code for d in result.diagnostics}


def test_valid_xtb_is_read_only_and_deterministic(tmp_path, opt_request):
    xtb_output(tmp_path)
    execution = receipt(opt_request, tmp_path)
    before = snapshot_artifacts(tmp_path)
    first = inspect_optimization(opt_request, tmp_path, execution=execution)
    second = inspect_optimization(opt_request, tmp_path, execution=execution)
    assert first.status is Outcome.CONVERGED and first.converged
    assert first.accepted_frame.physical_energy_hartree == -1.
    assert first.accepted_frame.coordinates_angstrom == opt_request.coordinates_angstrom
    assert not first.workflow_accepted
    assert first.to_json() == second.to_json()
    assert snapshot_artifacts(tmp_path) == before


@pytest.mark.parametrize('problem', ['missing_receipt', 'changed_method', 'changed_artifact', 'stale_marker', 'nonzero_exit'])
def test_xtb_false_success_is_rejected(tmp_path, opt_request, problem):
    xtb_output(tmp_path, converged=problem != 'stale_marker')
    (tmp_path / '.xtboptok').touch()
    execution = receipt(opt_request, tmp_path, returncode=1 if problem == 'nonzero_exit' else 0)
    if problem == 'missing_receipt':
        execution = None
    elif problem == 'changed_method':
        execution = ExecutionObservation('0' * 64, execution.artifacts, execution.processes, True, legacy_status=True)
    elif problem == 'changed_artifact':
        (tmp_path / 'xtbopt.xyz').write_text(xyz(energy=-2.))
    result = inspect_optimization(opt_request, tmp_path, execution=execution)
    assert not result.converged
    assert result.diagnostics


def test_xtb_scc_failure_wins_over_success_marker(tmp_path, opt_request):
    xtb_output(tmp_path)
    (tmp_path / '.xtboptok').touch()
    (tmp_path / '.sccnotconverged').touch()
    result = inspect_optimization(opt_request, tmp_path, execution=receipt(opt_request, tmp_path))
    assert not result.converged
    assert FailureCode.CALC_SCF_MAX_STEPS in codes(result)


@pytest.mark.parametrize('value', [0., -1.])
def test_real_zero_energy_is_valid(tmp_path, opt_request, value):
    xtb_output(tmp_path, energy=value)
    assert inspect_optimization(opt_request, tmp_path, execution=receipt(opt_request, tmp_path)).converged


@pytest.mark.parametrize('problem', ['missing_energy', 'nan_energy', 'nan_coordinates', 'wrong_symbols', 'truncated', 'unknown_version'])
def test_invalid_native_artifacts_fail_closed(tmp_path, opt_request, problem):
    xtb_output(tmp_path)
    if problem == 'missing_energy':
        path = tmp_path / 'xtb.out'
        path.write_text(path.read_text().replace('| TOTAL ENERGY -1.0 Eh |\n', ''))
    elif problem == 'nan_energy':
        xtb_output(tmp_path, energy=float('nan'))
    elif problem == 'nan_coordinates':
        (tmp_path / 'xtbopt.xyz').write_text(xyz(coordinates=[[float('nan'),0,0],[0,0,.74]]))
    elif problem == 'wrong_symbols':
        (tmp_path / 'xtbopt.xyz').write_text(xyz(symbols=('He', 'He')))
    elif problem == 'truncated':
        (tmp_path / 'xtbopt.xyz').write_text('2\nenergy:-1\nH 0 0 0\n')
    else:
        path = tmp_path / 'xtb.out'
        path.write_text(path.read_text().replace('6.7.1', '99.0.0'))
    assert not inspect_optimization(opt_request, tmp_path, execution=receipt(opt_request, tmp_path)).converged


def orca_request(molecule):
    return OptimizationRequest.from_molecule(molecule, dict(software='orca', method='HF', basis='def2-SVP',
                                                           opt_threshold='normal'))


def orca_output(root, *, converged=True, energy='-1.0'):
    (root / 'trial_h2.inp').write_text('! HF def2-SVP Opt\n*xyz 0 1\nH 0 0 0\nH 0 0 .74\n*\n')
    log = 'Program Version 6.1.1\n'
    if converged:
        log += 'THE OPTIMIZATION HAS CONVERGED\n'
    log += 'RMS gradient 0.00001 0.0001 YES\nMAX gradient 0.00002 0.0003 YES\n'
    log += 'CARTESIAN COORDINATES (ANGSTROEM)\n------------------\nH 0 0 0\nH 0 0 .74\n\n'
    if energy is not None:
        log += f'FINAL SINGLE POINT ENERGY {energy}\n'
    log += '****ORCA TERMINATED NORMALLY****\n'
    (root / 'trial_h2.out').write_text(log)
    (root / 'trial_h2.xyz').write_text(xyz())


@pytest.mark.parametrize('converged,energy,expected', [(False,'-1',False),(True,None,False),(True,'nan',False),(True,'0.0',True),(True,'-1',True)])
def test_orca_requires_convergence_and_real_energy(tmp_path, molecule, converged, energy, expected):
    req = orca_request(molecule)
    orca_output(tmp_path, converged=converged, energy=energy)
    assert inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path)).converged is expected


def geometric_request(molecule, **settings):
    return OptimizationRequest.from_molecule(molecule, dict(software='xtb', xtb_model='gfn2',
        geometry_optimizer='geometric', opt_threshold='normal', **settings))


def geometric_output(root, *, converged=True, checkpoint=False):
    (root / 'trial_h2.xyz').write_text(xyz())
    (root / 'trial_h2_optim.xyz').write_text(xyz())
    (root / 'geometric.out').write_text('geomeTRIC started. Version: 1.1.1\n'
        'Step 0 : Gradient = 0.0/0.0 (rms/max) Energy = -1.0\n' +
        ('Converged! =D\n' if converged else 'Maximum iterations reached (1)\n'))
    state = dict(software='xtb', gamma=0., opt_target='minimum', positions_angstrom=[[0,0,0],[0,0,.74]],
                 backend_energy_hartree=-1., total_energy_hartree=-1., forces_ev_per_angstrom=[[0,0,0],[0,0,0]])
    (root / 'pyar_geometric_state.json').write_text(json.dumps(state))
    if checkpoint:
        (root / 'pyar_bias_controller_state.json').write_text(json.dumps(dict(symbols=['H','H'],
            positions_angstrom=state['positions_angstrom'], configuration={})))


def test_geometric_convergence_requires_matching_finite_state(tmp_path, molecule):
    req = geometric_request(molecule)
    geometric_output(tmp_path)
    assert inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path)).converged
    state = json.loads((tmp_path / 'pyar_geometric_state.json').read_text())
    state['positions_angstrom'][1][2] = 3.
    (tmp_path / 'pyar_geometric_state.json').write_text(json.dumps(state))
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path))
    assert result.status is Outcome.INVALID_RESULT


def test_unverified_geometric_limit_is_never_converged_or_accepted(tmp_path, molecule):
    req = geometric_request(molecule)
    geometric_output(tmp_path, converged=False)
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path, returncode=1, legacy_status='CycleExceeded'))
    assert not result.converged and result.accepted_frame is None
    assert FailureCode.ARTIFACT_INCOMPLETE in codes(result)


def test_input_identity_and_numerical_fingerprints_are_separate(molecule):
    a = geometric_request(molecule, opt_cycles=10)
    b = geometric_request(molecule, opt_cycles=20)
    assert a.scientific_fingerprint == b.scientific_fingerprint
    assert a.request_fingerprint != b.request_fingerprint
    molecule.charge = 1
    c = geometric_request(molecule, opt_cycles=10)
    assert c.scientific_fingerprint != a.scientific_fingerprint


def test_settings_and_input_are_immutable_snapshots(molecule):
    settings = dict(software='xtb', custom_keywords=['a'])
    req = OptimizationRequest.from_molecule(molecule, settings)
    settings['custom_keywords'].append('b')
    molecule.coordinates[0][0] = 9.
    copy = req.settings
    copy['custom_keywords'].append('c')
    assert req.settings['custom_keywords'] == ['a']
    assert req.coordinates_angstrom[0][0] == 0.


def test_inspection_missing_directory_does_not_create_it(tmp_path, opt_request):
    root = tmp_path / 'absent'
    assert not inspect_optimization(opt_request, root).converged
    assert not root.exists()


def test_inspection_rejects_symlink_escape(tmp_path, opt_request):
    outside = tmp_path / 'outside'
    outside.write_text('private')
    root = tmp_path / 'run'
    root.mkdir()
    (root / 'xtb.out').symlink_to(outside)
    assert not inspect_optimization(opt_request, root).converged


def test_legacy_adapter_never_certifies_a_boolean():
    from pyar.recovery.models import legacy_diagnostic
    assert legacy_diagnostic(True).code is FailureCode.ARTIFACT_INCOMPLETE
    assert legacy_diagnostic('CycleExceeded').code is FailureCode.OPT_MAX_STEPS
    assert legacy_diagnostic('SCFFailed').code is FailureCode.CALC_SCF_MAX_STEPS
    assert legacy_diagnostic(None).code is FailureCode.UNKNOWN


def test_geometric_accepted_checkpoint_remains_provisional(tmp_path, molecule):
    from pyar.biases.controller import BiasController
    from pyar.biases.afir import alpha_from_gamma
    molecule.fragments = [[0], [1]]
    req = geometric_request(molecule, gamma=100., bias_controller='adaptive')
    geometric_output(tmp_path, converged=False, checkpoint=True)
    checkpoint = tmp_path / 'pyar_bias_controller_state.json'
    value = json.loads(checkpoint.read_text())
    settings = req.settings
    settings.pop('fragments')
    controller = BiasController('adaptive')
    controller.start_segment([[0.,0.,0.],[0.,0.,0.]], [[0.,0.,1.],[0.,0.,-1.]], 1., alpha_from_gamma(100.))
    value['controller'] = controller.state_dict()
    value['configuration'] = dict(qc_params=settings, fragment_indices=[[0],[1]], opt_target='minimum')
    checkpoint.write_text(json.dumps(value))
    state_path = tmp_path / 'pyar_geometric_state.json'
    state = json.loads(state_path.read_text())
    state.update(gamma=100., bias_controller=value['controller'])
    state_path.write_text(json.dumps(state))
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path, returncode=1, legacy_status='CycleExceeded'))
    assert result.status is Outcome.INCOMPLETE_USABLE
    assert not result.converged and result.accepted_frame is not None
    assert FailureCode.OPT_MAX_STEPS in codes(result)


@pytest.mark.parametrize('problem', ['forces_shape', 'nan_forces', 'force_mismatch', 'bias_mismatch', 'wall_mismatch', 'objective_mismatch'])
def test_geometric_forces_and_objective_are_validated(tmp_path, molecule, problem):
    req = geometric_request(molecule)
    geometric_output(tmp_path)
    path = tmp_path / 'pyar_geometric_state.json'
    state = json.loads(path.read_text())
    if problem == 'forces_shape':
        state['forces_ev_per_angstrom'] = [[0,0,0]]
    elif problem == 'nan_forces':
        state['forces_ev_per_angstrom'][0][0] = float('nan')
    elif problem == 'force_mismatch':
        state['forces_ev_per_angstrom'][0][0] = 1.
    elif problem == 'bias_mismatch':
        state['gamma'] = 100.
    elif problem == 'wall_mismatch':
        state['wall_configuration'] = {'radius_angstrom': 3.}
    else:
        state['backend_energy_hartree'] = -2.
    path.write_text(json.dumps(state))
    assert not inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path)).converged


def test_native_model_mismatch_is_rejected(tmp_path, opt_request):
    xtb_output(tmp_path)
    path = tmp_path / 'xtb.out'
    path.write_text(path.read_text().replace('GFN2-xTB', 'GFN1-xTB'))
    assert not inspect_optimization(opt_request, tmp_path, execution=receipt(opt_request, tmp_path)).converged


def test_original_exception_and_cause_survive_observation(tmp_path, opt_request):
    xtb_output(tmp_path)
    inner = ArithmeticError('physical calculator failed')
    outer = RuntimeError('optimizer boundary')
    outer.__cause__ = inner
    result = inspect_optimization(opt_request, tmp_path, execution=receipt(opt_request, tmp_path), exception=outer)
    assert not result.converged and result.exception is outer
    assert result.to_dict()['exception_chain'] == [
        {'type': 'RuntimeError', 'message': 'optimizer boundary'},
        {'type': 'ArithmeticError', 'message': 'physical calculator failed'},
    ]


def test_asset_identity_is_frozen_and_portable(tmp_path, molecule):
    first = tmp_path / 'first.model'
    second = tmp_path / 'second.model'
    first.write_bytes(b'model bytes')
    second.write_bytes(first.read_bytes())
    a = OptimizationRequest.from_molecule(molecule, dict(software='aimnet_2', model=str(first)))
    b = OptimizationRequest.from_molecule(molecule, dict(software='aimnet_2', model=str(second)))
    old = a.scientific_fingerprint
    assert old == b.scientific_fingerprint
    first.write_bytes(b'changed model')
    assert old == a.scientific_fingerprint


def test_accepted_trends_ignore_rejected_and_unknown_frames():
    from pyar.recovery.models import StepObservation
    from pyar.recovery.trajectory import analyze_accepted_steps
    steps = (StepObservation(0, -1., .1, .2, True),
             StepObservation(1, -100., .001, .002, False),
             StepObservation(2, -2., .05, .1, True),
             StepObservation(3, -99., .0001, .0002, None))
    trend = analyze_accepted_steps(steps)
    assert (trend.accepted_steps, trend.rejected_steps, trend.unknown_steps) == (2,1,1)
    assert trend.effective_energy_change_hartree == -1.
    assert trend.gradient_rms_change_hartree_per_bohr == -.05


def test_existing_output_never_triggers_execution(tmp_path, monkeypatch, molecule):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'job_h2').mkdir()
    sentinel = tmp_path / 'job_h2/result_h2.xyz'
    sentinel.write_text(xyz())
    from pyar import optimiser
    def forbidden(*args):
        pytest.fail('existing cache must not execute')
    monkeypatch.setattr(optimiser, 'optimise', forbidden)
    result = optimize_once(molecule, {'software': 'xtb'})
    assert not result.converged
    assert sentinel.read_text() == xyz()
    assert not (tmp_path / '.job_h2.pyar-v2.lock').exists()


def test_unsupported_target_is_rejected_before_execution(tmp_path, monkeypatch, molecule):
    monkeypatch.chdir(tmp_path)
    from pyar import optimiser
    monkeypatch.setattr(optimiser, 'optimise', lambda *a: pytest.fail('unsupported target executed'))
    result = optimize_once(molecule, dict(software='xtb', opt_target='saddle'))
    assert not result.converged
    assert list(tmp_path.iterdir()) == []


def test_facade_executes_once_and_preserves_input(tmp_path, monkeypatch, molecule):
    import copy
    from pyar import optimiser
    from pyar.recovery.observations import record_process
    monkeypatch.chdir(tmp_path)
    before = copy.deepcopy(molecule)
    settings = dict(software='xtb', opt_threshold='normal')
    calls = []
    def execute(candidate, params):
        calls.append((candidate, params))
        root = Path('job_h2')
        root.mkdir()
        xtb_output(root)
        candidate.coordinates = None
        params['extra'] = 'mutated'
        record_process(['xtb'], 0)
        return True
    monkeypatch.setattr(optimiser, 'optimise', execute)
    result = optimize_once(molecule, settings)
    assert len(calls) == 1 and result.converged
    assert molecule.coordinates == before.coordinates
    assert settings == dict(software='xtb', opt_threshold='normal')
    assert result.to_dict()['attempts'][0]['request']['backend'] == 'xtb'


def test_diagnostic_import_has_no_optional_scientific_dependencies(tmp_path):
    import subprocess, sys
    script = '''
import sys, importlib.abc
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, *args):
        if fullname.split('.')[0] in {'ase','geometric','sella','torch','mlatom','rdkit'}:
            raise ImportError('blocked optional import: '+fullname)
sys.meta_path.insert(0, Block())
import pyar.optimization_service
import pyar.recovery.models
'''
    completed = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize('backend', ['xtb', 'orca'])
def test_native_bias_is_not_silently_ignored(tmp_path, molecule, backend):
    req = OptimizationRequest.from_molecule(molecule, dict(software=backend, gamma=100.))
    if backend == 'xtb':
        xtb_output(tmp_path)
    else:
        orca_output(tmp_path)
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path))
    assert not result.converged
    assert FailureCode.BACKEND_UNSUPPORTED in codes(result)


def test_missing_process_status_is_never_certified(tmp_path, opt_request):
    xtb_output(tmp_path)
    result = inspect_optimization(opt_request, tmp_path, execution=receipt(opt_request, tmp_path, returncode=None))
    assert not result.converged


def test_observer_is_context_local_and_nested():
    from pyar.recovery.observations import capture_processes, record_process
    with capture_processes() as first:
        record_process(['xtb'], 7)
        with capture_processes() as second:
            record_process(['orca'], 0)
        record_process(['geometric-optimize'], 1)
    record_process(['unobserved'], 0)
    assert [(p.program,p.returncode) for p in first] == [('xtb',7),('geometric-optimize',1)]
    assert [(p.program,p.returncode) for p in second] == [('orca',0)]


def test_inspection_detects_artifact_changed_during_parse(tmp_path, opt_request, monkeypatch):
    import pyar.optimization_service as service
    xtb_output(tmp_path)
    execution = receipt(opt_request, tmp_path)
    original = service.parse_observation
    def changing(req, path):
        observed = original(req, path)
        path('xtbopt.xyz').write_text(xyz(energy=-9.))
        return observed
    monkeypatch.setattr(service, 'parse_observation', changing)
    assert not service.inspect_optimization(opt_request, tmp_path, execution=execution).converged


def test_single_atom_geometry_diagnostics_keep_matrix_shape(tmp_path):
    atom = SimpleNamespace(name='h2', atoms_list=['He'], coordinates=[[0.,0.,0.]], charge=0,
                           multiplicity=1, scftype='rhf', fragments=[])
    req = OptimizationRequest.from_molecule(atom, dict(software='xtb', opt_threshold='normal'))
    xtb_output(tmp_path)
    (tmp_path / 'trial_h2.xyz').write_text(xyz(symbols=('He',)))
    (tmp_path / 'xtbopt.xyz').write_text(xyz(symbols=('He',)))
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path))
    assert result.converged
    assert result.accepted_frame.coordinates_angstrom == ((0.,0.,0.),)


def test_fake_xtb_process_runs_through_the_opt_in_facade(tmp_path, monkeypatch):
    import copy, sys
    from pyar.backends import xtb
    from pyar.core.molecule import Molecule
    molecule = Molecule(['H','H'], [[0,0,0],[0,0,.74]], name='h2')
    original = copy.deepcopy(molecule.coordinates)
    executable = tmp_path / 'fake-xtb'
    executable.write_text(f'#!{sys.executable}\nfrom pathlib import Path\n' +
        'Path("xtbopt.xyz").write_text(' + repr(xyz()) + ')\n' +
        'Path(".xtboptok").touch()\n' +
        'print("xtb version 6.7.1\\nHamiltonian GFN2-xTB\\n*** GEOMETRY OPTIMIZATION CONVERGED AFTER 2 ITERATIONS ***\\n| TOTAL ENERGY -1.0 Eh |\\n| GRADIENT NORM 0.00001 Eh/a0 |\\nnormal termination of xtb")\n')
    executable.chmod(0o755)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(xtb, 'require_executable', lambda *args: str(executable))
    result = optimize_once(molecule, dict(software='xtb', opt_threshold='normal'))
    assert result.converged
    assert result.attempts[0].execution.processes[0].returncode == 0
    assert result.accepted_frame.physical_energy_hartree == -1.
    assert (molecule.coordinates == original).all()
    assert not (tmp_path / 'job_h2/.xtboptok').exists()  # Native adapter owns marker consumption.


def test_historical_output_can_be_diagnosed_without_certifying_current_run(tmp_path, opt_request):
    (tmp_path / 'xtb.out').write_text('xtb version 6.7.1\nSCC failed to converge\n')
    before = snapshot_artifacts(tmp_path)
    result = inspect_optimization(opt_request, tmp_path)
    assert not result.converged and result.accepted_frame is None
    assert FailureCode.CALC_SCF_MAX_STEPS in codes(result)
    assert FailureCode.CHECKPOINT_INCOMPATIBLE in codes(result)
    assert result.attempts[0].artifacts == before
    assert result.to_json() == inspect_optimization(opt_request, tmp_path).to_json()
    assert snapshot_artifacts(tmp_path) == before


def test_unknown_backend_preflight_does_not_create_output(tmp_path, monkeypatch, molecule):
    monkeypatch.chdir(tmp_path)
    result = optimize_once(molecule, {'software': 'unavailable-backend'})
    assert FailureCode.BACKEND_UNSUPPORTED in codes(result)
    assert list(tmp_path.iterdir()) == []


def test_geometric_limit_never_accepts_an_incompatible_checkpoint(tmp_path, molecule):
    req = geometric_request(molecule)
    geometric_output(tmp_path, converged=False, checkpoint=True)
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path, returncode=1))
    assert result.accepted_frame is None and not result.converged


def test_unknown_state_version_and_output_patch_are_not_adopted(tmp_path, molecule):
    req = geometric_request(molecule)
    geometric_output(tmp_path)
    log = tmp_path / 'geometric.out'
    log.write_text(log.read_text().replace('1.1.1', '1.1.99'))
    assert not inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path)).converged
    geometric_output(tmp_path)
    state_path = tmp_path / 'pyar_geometric_state.json'
    state = json.loads(state_path.read_text())
    state['schema_version'] = 99
    state_path.write_text(json.dumps(state))
    assert not inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path)).converged


def test_zero_exit_unknown_option_is_input_error_not_scf_diagnosis(tmp_path, opt_request):
    (tmp_path / 'xtb.out').write_text("xtb version 6.7.1\nUnknown option '--gxtb' provided\nnormal termination of xtb\n")
    result = inspect_optimization(opt_request, tmp_path)
    assert not result.converged
    assert FailureCode.USER_INPUT_INVALID in codes(result)
    assert FailureCode.CALC_SCF_MAX_STEPS not in codes(result)


def test_caught_native_process_error_is_retained(tmp_path, monkeypatch):
    import subprocess, sys
    from pyar.backends import xtb
    from pyar.core.molecule import Molecule
    executable = tmp_path / 'failing-xtb'
    executable.write_text(f'#!{sys.executable}\nimport sys\nprint("xtb version 6.7.1\\nSCC failed to converge")\nsys.exit(9)\n')
    executable.chmod(0o755)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(xtb, 'require_executable', lambda *args: str(executable))
    molecule = Molecule(['H','H'], [[0,0,0],[0,0,.74]], name='h2')
    result = optimize_once(molecule, dict(software='xtb', opt_threshold='normal'))
    assert result.status is Outcome.FAILED_FINAL
    process = result.attempts[0].execution.processes[0]
    assert process.returncode == 9
    assert isinstance(process.exception, subprocess.CalledProcessError)
    assert result.exception is not None
    assert 'CalledProcessError' in result.to_json()


def test_supported_capabilities_are_observational_only():
    from pyar.backend_capabilities import optimization_diagnostic_capabilities
    assert optimization_diagnostic_capabilities('xtb') == frozenset({'parsed_convergence', 'frame_validation'})
    assert optimization_diagnostic_capabilities('orca')
    assert optimization_diagnostic_capabilities('xtb', 'geometric')
    assert not optimization_diagnostic_capabilities('aimnet_2', 'geometric')
    assert not optimization_diagnostic_capabilities('gaussian')
    assert 'compatible_checkpoint' not in optimization_diagnostic_capabilities('xtb')


def test_optimizer_limit_wins_over_success_text(tmp_path, molecule):
    req = orca_request(molecule)
    orca_output(tmp_path)
    path = tmp_path / 'trial_h2.out'
    path.write_text(path.read_text() + '\nTHE OPTIMIZATION DID NOT CONVERGE: maximum number of geometry optimization cycles reached\n')
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path))
    assert not result.converged
    assert FailureCode.OPT_MAX_STEPS in codes(result)


def test_explicit_constraints_never_execute_unconstrained(tmp_path, monkeypatch, molecule):
    monkeypatch.chdir(tmp_path)
    result = optimize_once(molecule, dict(software='xtb', constraints={'distance': [0,1,.74]}))
    assert not result.converged
    assert list(tmp_path.iterdir()) == []


def test_rejected_final_trial_cannot_be_an_accepted_checkpoint(tmp_path, molecule):
    req = geometric_request(molecule)
    geometric_output(tmp_path, converged=False, checkpoint=True)
    checkpoint = tmp_path / 'pyar_bias_controller_state.json'
    value = json.loads(checkpoint.read_text())
    settings = req.settings
    settings.pop('fragments')
    value['configuration'] = dict(qc_params=settings, fragment_indices=[], opt_target='minimum')
    value['positions_angstrom'][1][2] = .8  # last accepted frame differs from final trial
    checkpoint.write_text(json.dumps(value))
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path, returncode=1))
    assert result.accepted_frame is None and not result.converged


def test_physical_energy_is_separate_from_biased_objective(tmp_path, molecule):
    req = geometric_request(molecule, gamma=100.)
    geometric_output(tmp_path)
    path = tmp_path / 'pyar_geometric_state.json'
    value = json.loads(path.read_text())
    value.update(gamma=100., backend_energy_hartree=-2., total_energy_hartree=-1.)
    path.write_text(json.dumps(value))
    result = inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path))
    assert result.converged  # Synthetic stationary effective-objective evidence, not a physical minimum claim.
    assert result.accepted_frame.physical_energy_hartree == -2.
    assert result.accepted_frame.effective_energy_hartree == -1.
    assert not result.workflow_accepted


def test_active_run_writer_is_respected(tmp_path, monkeypatch, molecule):
    from pyar.persistence import RunLock
    monkeypatch.chdir(tmp_path)
    root = tmp_path / 'job_h2'
    with RunLock(root):
        metadata = (tmp_path / '.job_h2.pyar-v2.lock').read_bytes()
        result = optimize_once(molecule, dict(software='xtb'))
        assert not result.converged
        assert not root.exists()
        assert (tmp_path / '.job_h2.pyar-v2.lock').read_bytes() == metadata


def test_unsafe_name_is_rejected_before_output(tmp_path, monkeypatch, molecule):
    monkeypatch.chdir(tmp_path)
    molecule.name = '../escape'
    with pytest.raises(ValueError, match='filename stem'):
        optimize_once(molecule, dict(software='xtb'))
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('change', ['fractional_charge', 'zero_multiplicity', 'conflicting_charge', 'conflicting_spin'])
def test_invalid_or_conflicting_electronic_state_is_not_rewritten(tmp_path, monkeypatch, molecule, change):
    monkeypatch.chdir(tmp_path)
    settings = {'software': 'xtb'}
    if change == 'fractional_charge':
        molecule.charge = .5
    elif change == 'zero_multiplicity':
        molecule.multiplicity = 0
    elif change == 'conflicting_charge':
        settings['charge'] = 1
    else:
        settings['xtb_unpaired_electrons'] = 2
    with pytest.raises(ValueError):
        optimize_once(molecule, settings)
    assert list(tmp_path.iterdir()) == []


def test_rounded_gradient_log_cannot_loosen_original_tolerance(tmp_path, molecule):
    req = OptimizationRequest.from_molecule(molecule, dict(software='xtb', geometry_optimizer='geometric',
                                                         opt_threshold='tight', xtb_model='gfn2'))
    geometric_output(tmp_path)
    log = tmp_path / 'geometric.out'
    log.write_text(log.read_text().replace('Gradient = 0.0/0.0', 'Gradient = 1.0e-5/1.0e-5'))
    state_path = tmp_path / 'pyar_geometric_state.json'
    state = json.loads(state_path.read_text())
    force = 1.000001e-5 * 27.211386245988 / .529177210903
    state['forces_ev_per_angstrom'] = [[0,0,force],[0,0,-force]]
    state_path.write_text(json.dumps(state))
    assert not inspect_optimization(req, tmp_path, execution=receipt(req, tmp_path)).converged


@pytest.mark.parametrize('legacy_status', [False, 'CycleExceeded', 'SCFFailed'])
def test_legacy_rejection_is_not_promoted_by_success_text(tmp_path, opt_request, legacy_status):
    xtb_output(tmp_path)
    result = inspect_optimization(opt_request, tmp_path,
                                 execution=receipt(opt_request, tmp_path, legacy_status=legacy_status))
    assert not result.converged
    assert result.diagnostics
