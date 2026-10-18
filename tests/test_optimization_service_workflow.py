import json
from pathlib import Path

from pyar.core.molecule import Molecule
from pyar.optimization_service import optimize_workflow_once
from pyar.recovery.models import Outcome


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
