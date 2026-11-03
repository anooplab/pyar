"""The supported request API executes discovery without trusting restart metadata."""
from pyar.core.molecule import Molecule
from pyar.reaction_request import ReactionRequest
from pyar.workflows.reaction import run_reaction


def test_request_discovery_preserves_caller_settings(tmp_path, deterministic_workflow_backend):
    settings = dict(software='xtb', xtb_model='gfn2', gamma=None)
    request = ReactionRequest([
        Molecule(['He'], [[0, 0, 0]], name='A'),
        Molecule(['Ne'], [[0, 0, 0]], name='B')], settings, 100., 100., 1, None, 2.3,
        {'caller_metadata_is_not_restart_authority': True})
    result = run_reaction(request, directory=tmp_path)
    assert result.status == 'completed'
    assert settings == dict(software='xtb', xtb_model='gfn2', gamma=None)
    import json
    state = json.loads((tmp_path / 'reaction/state.json').read_text())
    assert 'caller_metadata_is_not_restart_authority' not in state['request']
    assert state['completed_jobs']
