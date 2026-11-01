"""Supported request-taking APIs preserve persisted scientific requests."""
from pyar.aggregation.request import AggregateRequest
from pyar.core.molecule import Molecule
from pyar.growth.request import GrowRequest
from pyar.persistence import read_json
from pyar.workflows.aggregate import run_aggregate
from pyar.workflows.grow import run_grow


def molecules():
    return Molecule(['He'], [[0,0,0]], name='first'), Molecule(['Ne'], [[0,0,0]], name='second')


def test_resolved_aggregate_api_preserves_request_and_selected_composition(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    request = AggregateRequest.from_options(list(molecules()), [1,1], 1, {}, 1, 0, 1, None, 'off')
    original = request.to_state_dict()
    result = run_aggregate(request)
    assert result.status == 'completed'
    assert read_json(result.state_path)['request'] == original
    assert len(result.selected_paths) == 1
    selected = Molecule.from_xyz(str(tmp_path/'aggregates'/result.selected_paths[0]))
    assert sorted(selected.atoms_list) == ['He','Ne']
    assert selected.charge == 0 and selected.multiplicity == 1


def test_resolved_grow_api_preserves_request_and_completed_restart(tmp_path):
    seed, monomer = molecules()
    request = GrowRequest(seed, monomer, 1, number_of_orientations=1,
                          maximum_number_of_seeds=1, connectivity_policy='off')
    result = run_grow(request, output=tmp_path/'grow')
    assert result.status == 'completed'
    assert read_json(result.state_path)['request'] == request.to_state_dict()
    repeated = run_grow(request, output=tmp_path/'grow')
    assert repeated.selected_paths == result.selected_paths
