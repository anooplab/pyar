"""Stop signals cannot be mistaken for a pool or commit a partial addition."""
import importlib

import pytest

from pyar.core.molecule import Molecule
from pyar.growth.errors import GrowthStopped
from pyar.growth.request import GrowRequest
from pyar.growth import service
from pyar.persistence import read_json

from pyar.growth.policy import GrowthPolicy, SelectionOptions


def inputs():
    return Molecule(['He'], [[0,0,0]], name='seed'), Molecule(['Ne'], [[0,0,0]], name='monomer')


def test_addition_stop_marker_raises_without_starting_work(tmp_path):
    (tmp_path/'stop').touch()
    seed, monomer = inputs()
    with pytest.raises(GrowthStopped):
        service.add_one('addition', [seed], monomer, GrowthPolicy(1, {}, 1, None, None, SelectionOptions('auto', 'auto', 'euclidean', 'auto'), True), directory=tmp_path)
    assert sorted(path.name for path in tmp_path.iterdir()) == ['stop']


def test_aggregate_stop_keeps_checkpoint_pending_and_resumable(tmp_path, monkeypatch):
    workflow = importlib.import_module('pyar.workflows.aggregate')
    monkeypatch.chdir(tmp_path)
    def stop(*args, **kwargs):
        raise GrowthStopped('owner stop')
    monkeypatch.setattr(workflow, 'add_one', stop)
    seed, monomer = inputs()
    result = workflow.aggregate([seed,monomer], [1,1], 1, {}, 1, 0, 1, None)
    assert result.status == 'stopped'
    root = tmp_path/'aggregates'
    assert read_json(root/'state.json')['status'] == 'running'
    assert len(list(root.glob('pathways/pathway_*/steps/step_*/selection.json'))) == 1
    assert not (root/'final_selection.json').exists()
    monkeypatch.setattr(workflow, 'add_one', service.add_one)
    seed, monomer = inputs()
    resumed = workflow.aggregate([seed,monomer], [1,1], 1, {}, 1, 0, 1, None)
    assert resumed.status == 'completed'


def test_fixed_addend_stop_preserves_state_and_resumes(tmp_path, monkeypatch):
    workflow = importlib.import_module('pyar.workflows.grow')
    seed, monomer = inputs()
    request = GrowRequest(seed, monomer, 1, number_of_orientations=1, maximum_number_of_seeds=1,
                          connectivity_policy='off')
    def stop(*args, **kwargs):
        raise GrowthStopped('owner stop')
    monkeypatch.setattr(workflow, 'add_one', stop)
    result = workflow.grow(request, output=tmp_path/'grow')
    assert result.status == 'stopped'
    state = read_json(tmp_path/'grow/state.json')
    assert state['status'] == 'stopped' and not state['completed_steps']
    monkeypatch.setattr(workflow, 'add_one', service.add_one)
    assert workflow.grow(request, output=tmp_path/'grow').status == 'completed'
