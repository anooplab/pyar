import json
from pathlib import Path
import pytest

from pyar.core.molecule import Molecule
from pyar.optimization_service import optimize_workflow_once
from pyar.recovery.models import Outcome


@pytest.fixture
def workflow_case(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = Molecule(['H', 'H'], [[0., 0., 0.], [0., 0., .74]], name='h2')
    return molecule, dict(software='xtb', geometry_optimizer='geometric', opt_cycles=10,
                         opt_threshold='normal')


def test_sidecar_failure_preserves_success_and_returns_no_report(workflow_case, monkeypatch):
    molecule, settings = workflow_case
    def broken_write(*args, **kwargs):
        raise OSError('disk full')
    monkeypatch.setattr('pyar.optimization_service.write_json', broken_write)
    status, result, report = optimize_workflow_once(
        molecule, settings, workflow_role='test', optimizer=lambda *args: True)
    assert status is True
    assert report is None
    assert not result.converged


def test_real_sidecar_namespace_collision_is_nonfatal(workflow_case):
    molecule, settings = workflow_case
    Path('.pyar_optimization_diagnostics').write_text('existing unrelated file')
    status, _, report = optimize_workflow_once(
        molecule, settings, workflow_role='test', optimizer=lambda *args: True)
    assert status is True and report is None
    assert Path('.pyar_optimization_diagnostics').read_text() == 'existing unrelated file'


def test_sidecar_failure_preserves_original_exception(workflow_case, monkeypatch):
    molecule, settings = workflow_case
    original = RuntimeError('backend failure')
    def fail(*args):
        raise original
    def broken_write(*args, **kwargs):
        raise OSError('disk full')
    monkeypatch.setattr('pyar.optimization_service.write_json', broken_write)
    with pytest.raises(RuntimeError) as caught:
        optimize_workflow_once(molecule, settings, workflow_role='test', optimizer=fail)
    assert caught.value is original
    assert original.ore_diagnostic_path is None


@pytest.mark.parametrize('boundary', ['snapshot_artifacts', 'write_json'])
def test_observer_cancellation_propagates(workflow_case, monkeypatch, boundary):
    molecule, settings = workflow_case
    def cancel(*args, **kwargs):
        raise KeyboardInterrupt('cancel')
    monkeypatch.setattr('pyar.optimization_service.' + boundary, cancel)
    with pytest.raises(KeyboardInterrupt, match='cancel'):
        optimize_workflow_once(molecule, settings, workflow_role='test', optimizer=lambda *args: True)


def test_owner_permission_distinguishes_provisional_from_validated(workflow_case):
    molecule, settings = workflow_case
    for allow in (False, True):
        status, result, report = optimize_workflow_once(molecule, settings, workflow_role='test',
            optimizer=lambda *args: 'CycleExceeded', allow_incomplete=allow)
        decision = json.loads(report.read_text())['owner_decision']
        assert decision['accepted'] is allow
        assert decision['provisional'] is allow
        assert not decision['validated_converged']
        assert not result.workflow_accepted


@pytest.mark.parametrize('backend,marker,code', [
    ('gaussian', 'Convergence failure', 'CALC_SCF_MAX_STEPS'),
    ('psi4', 'SCFConvergenceError', 'CALC_SCF_MAX_STEPS'),
    ('turbomole', 'SCF NOT CONVERGED', 'CALC_SCF_MAX_STEPS'),
    ('mopac', 'SCF FIELD WAS NOT ACHIEVED', 'CALC_SCF_MAX_STEPS'),
    ('aimnet_2', 'CUDA out of memory', 'EXEC_FAILED'),
    ('ani', 'OutOfMemoryError', 'EXEC_FAILED'),
    ('aiqm1_mlatom', 'Convergence failure', 'CALC_SCF_MAX_STEPS'),
    ('xtb-aimnet2', 'CUDA out of memory', 'EXEC_FAILED'),
    ('xtb-aiqm1', 'Convergence failure', 'CALC_SCF_MAX_STEPS'),
])
def test_diagnose_only_routes_have_negative_evidence(workflow_case, backend, marker, code):
    molecule, settings = workflow_case
    settings.update(software=backend, geometry_optimizer='native')
    def execute(*args):
        Path('job_h2').mkdir()
        Path('job_h2/trial_h2.out').write_text(marker)
        return True  # A zero exit or legacy success must never certify convergence.
    status, result, report = optimize_workflow_once(
        molecule, settings, workflow_role='test', optimizer=execute)
    assert status is True
    assert not result.converged
    assert result.accepted_frame is None
    assert code in {d.code.value for d in result.diagnostics}
    assert result.attempts[0].artifacts


@pytest.mark.parametrize('backend', [
    'gaussian', 'psi4', 'turbomole', 'xtb_turbo', 'mopac', 'aimnet_2', 'ani',
    'aiqm1_mlatom', 'obabel', 'xtb-aimnet2', 'xtb-aiqm1',
])
def test_diagnose_only_success_marker_is_never_a_certificate(workflow_case, backend):
    from pyar.backend_capabilities import optimization_diagnostic_mode
    molecule, settings = workflow_case
    settings.update(software=backend, geometry_optimizer='native')
    def execute(*args):
        Path('job_h2').mkdir()
        Path('job_h2/trial_h2.out').write_text('Normal termination\nOptimization completed\n')
        return True
    _, result, _ = optimize_workflow_once(molecule, settings, workflow_role='test', optimizer=execute)
    assert optimization_diagnostic_mode(backend) == 'diagnose_only'
    assert not result.converged and result.accepted_frame is None


def test_workflow_diagnostic_preserves_owner_status_and_records_provisional_result(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = Molecule(['H', 'H'], [[0., 0., 0.], [0., 0., .74]], name='h2')
    settings = {'software': 'xtb', 'geometry_optimizer': 'geometric', 'opt_cycles': 10,
                'charge': 0, 'multiplicity': 1, 'scftype': 'rhf', 'opt_threshold': 'normal'}

    def incomplete(candidate, parameters):
        candidate.energy = -1.0
        return 'CycleExceeded'

    legacy_status, result, report = optimize_workflow_once(
        molecule, settings, workflow_role='growth_preselection', optimizer=incomplete)

    assert legacy_status == 'CycleExceeded'
    assert result.status is Outcome.INVALID_RESULT
    value = json.loads(Path(report).read_text())
    assert value['workflow_role'] == 'growth_preselection'
    assert value['legacy_status'] == 'cycle_exceeded_usable'
    assert value['result']['converged'] is False
    assert value['result']['workflow_accepted'] is False
    assert Path(report).is_relative_to(tmp_path / '.pyar_optimization_diagnostics')


def test_workflow_diagnostic_records_and_reraises_backend_error(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = Molecule(['H'], [[0., 0., 0.]], name='h')
    settings = {'software': 'xtb', 'geometry_optimizer': 'geometric', 'opt_cycles': 10,
                'charge': 0, 'multiplicity': 1, 'scftype': 'rhf', 'opt_threshold': 'normal'}

    def fail(candidate, parameters):
        raise RuntimeError('backend boundary failure')

    try:
        optimize_workflow_once(molecule, settings, workflow_role='conformer_refinement', optimizer=fail)
    except RuntimeError as exc:
        assert str(exc) == 'backend boundary failure'
    else:
        raise AssertionError('The owning workflow exception contract must be preserved')

    reports = list((tmp_path / '.pyar_optimization_diagnostics').rglob('*.json'))
    assert len(reports) == 1
    value = json.loads(reports[0].read_text())
    assert value['result']['converged'] is False
    assert value['result']['exception_chain'][0]['message'] == 'backend boundary failure'
