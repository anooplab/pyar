from types import SimpleNamespace

import numpy as np
import pytest

from pyar.recovery.stages import run_scan_once, scan_diagnostic
from pyar.recovery.workflow import stage_diagnostic


def scan_case():
    molecule = SimpleNamespace(atoms_list=['H', 'H'])
    request = SimpleNamespace(n_points=2, atom_i=0, atom_j=1,
                              start_distance_angstrom=1., end_distance_angstrom=2.)
    result = SimpleNamespace(success=True, status='success',
        frames=[(['H', 'H'], np.array([[0., 0., 0.], [0., 0., d]]), '') for d in (1., 2.)],
        profile=[dict(target_distance_angstrom=d, energy_hartree=-1.) for d in (1., 2.)])
    return molecule, request, result


@pytest.mark.parametrize('corruption', ['distance', 'mapping', 'energy', 'missing_point'])
def test_constrained_scan_rejects_false_success(tmp_path, corruption):
    molecule, request, result = scan_case()
    if corruption == 'distance':
        result.frames[1][1][1, 2] = 2.1
    elif corruption == 'mapping':
        result.frames[1][0][1] = 'He'
    elif corruption == 'energy':
        result.profile[1]['energy_hartree'] = float('nan')
    else:
        result.frames.pop()
    observed = run_scan_once(lambda *args: result, molecule, request, tmp_path, {})
    assert not observed.success
    assert observed.status == 'scan_validation_failed'
    assert (tmp_path / 'optimization_diagnostic.json').is_file()


def test_valid_constrained_scan_is_not_changed(tmp_path):
    molecule, request, result = scan_case()
    assert run_scan_once(lambda *args: result, molecule, request, tmp_path, {}) is result
    assert result.success
    diagnostic = scan_diagnostic(molecule, request, result)
    assert diagnostic['stage_gate'] == 'passed'
    assert not diagnostic['restart_authority']


def test_scan_exception_and_cancellation_preserved(tmp_path):
    molecule, request, _ = scan_case()
    for exception in (RuntimeError('calculator failed'), KeyboardInterrupt('cancel')):
        def fail(*args):
            raise exception
        with pytest.raises(type(exception)) as caught:
            run_scan_once(fail, molecule, request, tmp_path, {})
        assert caught.value is exception


def test_scan_rejects_unsupported_recovery_without_calling_backend(tmp_path):
    molecule, request, _ = scan_case()
    def forbidden(*args):
        pytest.fail('Unsupported recovery must fail before backend execution')
    with pytest.raises(ValueError, match='recovery_unsupported'):
        run_scan_once(forbidden, molecule, request, tmp_path, {'workflow_recovery': {'mode': 'safe'}})


def test_partially_initialized_optimizer_evidence_is_best_effort(tmp_path):
    from pyar.recovery.stages import preserve_band_observations, preserve_optimizer_observations
    preserve_optimizer_observations(tmp_path, 'ts', ['H'], SimpleNamespace())
    preserve_band_observations(tmp_path, ['H'], SimpleNamespace())
    assert not list(tmp_path.glob('*.xyz'))


def test_ts_boolean_is_not_a_saddle_or_connection_certificate():
    value = stage_diagnostic('neb_ts', dict(optimizer=True, saddle_hessian=None, connection=None))
    assert value['gates']['optimizer'] == 'passed'
    assert value['gates']['saddle_hessian'] == 'unknown'
    assert value['stage_gate'] == 'not_assessed'
    assert not value['recovery_supported']


def test_stage_retains_original_exception_chain():
    try:
        try:
            raise ValueError('bad gradient')
        except ValueError as exc:
            raise RuntimeError('Sella failed') from exc
    except RuntimeError as exc:
        value = stage_diagnostic('neb_ts', dict(optimizer=None), exception=exc)
    assert [item['type'] for item in value['exception_chain']] == ['RuntimeError', 'ValueError']
    assert value['reason'] == 'EXEC_FAILED'
