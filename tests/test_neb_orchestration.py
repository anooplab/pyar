"""Public failure and validation contracts across NEB orchestration boundaries."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pyar import neb

pytestmark = pytest.mark.usefixtures("inline_directory_worker")


@pytest.mark.parametrize('stage', ['all', 'relax'])
@pytest.mark.parametrize('exception_type', [RuntimeError, KeyboardInterrupt])
def test_interrupted_stage_preserves_failure_receipt_and_cwd(tmp_path, monkeypatch, stage, exception_type):
    calculator = SimpleNamespace(backend_energy_gradient_evaluations=7, _pyar_optimizer_steps_total=4)
    monkeypatch.setattr('pyar.backends.geometric.PyarGeometricCalculator', lambda **kwargs: calculator)
    error = exception_type('interrupted native stage')

    def interrupt(*args):
        calculator.backend_energy_gradient_evaluations += 3
        calculator._pyar_optimizer_steps_total += 2
        raise error

    monkeypatch.setattr(neb, '_execute_stage', interrupt)
    previous = Path.cwd()
    with pytest.raises(exception_type) as caught:
        neb.run_neb('start.xyz', 'end.xyz', 'guess.xyz', software='xtb', output=tmp_path, stage=stage)
    assert caught.value is error
    assert Path.cwd() == previous
    receipt = json.loads((tmp_path / 'relax_summary.json').read_text())
    assert receipt['status'] == 'failed'
    assert receipt['backend_energy_gradient_evaluations'] == 3
    assert receipt['optimizer_steps'] == 2
    assert receipt['error'] == str(error)
    workflow = tmp_path / 'workflow_summary.json'
    if stage == 'all':
        assert json.loads(workflow.read_text())['failed_stage'] == 'relax'
    else:
        assert not workflow.exists()


def test_route_error_precedes_numerical_controls_and_calculator_creation(tmp_path, monkeypatch):
    def forbidden(**kwargs):
        pytest.fail('invalid request reached calculator construction')

    monkeypatch.setattr('pyar.backends.geometric.PyarGeometricCalculator', forbidden)
    with pytest.raises(ValueError, match='images must be an odd integer'):
        neb.run_neb('start.xyz', 'end.xyz', 'guess.xyz', software='xtb', output=tmp_path,
                    images=2, max_cycles=0)
    assert not list(tmp_path.iterdir())
