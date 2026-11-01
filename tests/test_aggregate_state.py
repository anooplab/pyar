"""Scientific request and integrity guards for supported aggregation state."""
from copy import deepcopy

import pytest
from pyar.aggregation.request import AggregateRequest
from pyar.aggregation.layout_v2 import AggregateRunStateV2
from pyar.state.aggregate import AggregateStateError
from pyar.core.molecule import Molecule
from pyar.persistence import read_json, write_json


@pytest.fixture
def state(tmp_path):
    molecules = [Molecule(['C'], [[0, 0, 0]], name='a'), Molecule(['H'], [[0, 0, 0]], name='b')]
    request = AggregateRequest(tuple(molecules), (1, 2), 2,
                               backend_parameters={'software': 'xtb'}).to_state_dict()
    state = AggregateRunStateV2.create(tmp_path, request, molecules, ['abb', 'bab'], {})
    state.initialize_pathways()
    return state


def test_current_state_roundtrips_progress_read_only(state):
    state.complete_pathway(0, 'abb')
    before = state.state_file.read_bytes()
    loaded = AggregateRunStateV2.load(state.root_directory, state.data['request'])
    assert loaded.data['completed_pathways'] == [{'index': 0, 'label': 'abb', 'status': 'completed'}]
    assert state.state_file.read_bytes() == before


@pytest.mark.parametrize('key,value', [('maximum_number_of_seeds', 7), ('selection_policy_version', -1),
                                      ('backend_parameters', {'software': 'orca'})])
def test_scientific_request_change_rejected(state, key, value):
    changed = deepcopy(state.data['request'])
    changed[key] = value
    with pytest.raises(AggregateStateError, match='does not match'):
        AggregateRunStateV2.load(state.root_directory, changed)


def test_out_of_sequence_progress_rejected(state):
    state.complete_pathway(1, 'bab')
    with pytest.raises(AggregateStateError, match='out of order'):
        AggregateRunStateV2.load(state.root_directory, state.data['request'])


def test_modified_input_rejected(state):
    identity = read_json(state.aggregate_directory / 'request.json')
    path = state.aggregate_directory / identity['inputs'][0]['path']
    path.write_text('modified input')
    with pytest.raises(AggregateStateError, match='checksum'):
        AggregateRunStateV2.load(state.root_directory, state.data['request'])


def test_completed_state_cannot_be_reopened(state):
    record = read_json(state.state_file)
    record['status'] = 'completed'
    write_json(state.state_file, record)
    with pytest.raises(AggregateStateError, match='already'):
        AggregateRunStateV2.load(state.root_directory, state.data['request'])
