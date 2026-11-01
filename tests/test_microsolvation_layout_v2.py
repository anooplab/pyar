"""Microsolvation storage changes preserve original-solute shell science."""
import importlib
from pathlib import Path
import shutil
from unittest.mock import patch

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.microsolvation.request import MicrosolvationRequest
from pyar.microsolvation.layout_v2 import MicrosolvationRunStateV2
from pyar.persistence import read_json, sha256_file
from pyar.microsolvation.layout_v2 import MicrosolvationRunStateV2 as MicrosolvationRunState
from pyar.state.microsolvation import MicrosolvationStateError

workflow = importlib.import_module('pyar.workflows.microsolvation')
FIXTURE = Path(__file__).parent / 'data/run_layout_v1/microsolvation_interrupted'


def request():
    solute = Molecule(['C','H','H','H','H'],
        [[0,0,0],[.629,.629,.629],[-.629,-.629,.629],[-.629,.629,-.629],[.629,-.629,-.629]],
        name='methane', title='methane', charge=0, multiplicity=1)
    solvent = Molecule(['O','H','H'], [[0,0,0],[.9572,0,0],[-.239,.927,0]],
                       name='water', title='water', charge=0, multiplicity=1)
    return MicrosolvationRequest(solute, solvent, 2, 3, {}, 2, None,
                                  surface_points_per_atom=32, probe_radius=1.4, shell_tolerance=3.5)


def test_real_geometry_only_science_matches_v1_and_results_stay_immutable(tmp_path):
    new = workflow.microsolvate(request(), output=tmp_path / 'new')
    assert new.status == 'completed'
    assert new.selected_paths
    assert new.metadata['coverage_history']
    assert new.metadata['solute_atom_indices']
    assert new.metadata['solvent_fragments']
    before = {p: p.read_bytes() for p in (tmp_path / 'new').rglob('*') if p.is_file()}
    workflow.microsolvate(request(), output=tmp_path / 'new')
    assert before == {p: p.read_bytes() for p in (tmp_path / 'new').rglob('*') if p.is_file()}


def test_old_fixture_is_rejected_without_mutation(tmp_path):
    from pyar.persistence import LayoutError
    root = tmp_path / 'old'
    shutil.copytree(FIXTURE, root)
    before = {p.relative_to(root): p.read_bytes() for p in root.rglob('*') if p.is_file()}
    with pytest.raises(LayoutError, match='original PyAR environment'):
        workflow.microsolvate(request(), output=root)
    assert before == {p.relative_to(root): p.read_bytes() for p in root.rglob('*') if p.is_file()}


def test_crash_after_selection_publication_reuses_committed_step(tmp_path, monkeypatch):
    original = MicrosolvationRunStateV2.save
    failed = []
    def crash(state):
        if len(state.data['completed_steps']) == 1 and not failed:
            failed.append(True)
            raise OSError('crash after selection publish')
        return original(state)
    monkeypatch.setattr(MicrosolvationRunStateV2, 'save', crash)
    with pytest.raises(OSError, match='selection publish'):
        workflow.microsolvate(request(), output=tmp_path / 'run')
    record = tmp_path / 'run/steps/step_001/selection.json'
    digest = sha256_file(record)
    monkeypatch.setattr(MicrosolvationRunStateV2, 'save', original)
    generate = workflow.generate_microsolvation_candidates
    calls = []
    def record_call(*args, **kwargs):
        calls.append(kwargs['step'])
        return generate(*args, **kwargs)
    with patch.object(workflow, 'generate_microsolvation_candidates', side_effect=record_call):
        result = workflow.microsolvate(request(), output=tmp_path / 'run')
    assert result.status == 'completed' and calls and all(step >= 20000 for step in calls)
    assert sha256_file(record) == digest


def test_surface_tamper_rejected_and_manifest_can_be_rebuilt(tmp_path):
    root = tmp_path / 'run'
    workflow.microsolvate(request(), output=root)
    hashes = {p: sha256_file(p) for p in root.glob('steps/step_*/selection.json')}
    (root / 'results/manifest.json').unlink()
    workflow.microsolvate(request(), output=root)
    assert hashes == {p: sha256_file(p) for p in hashes}
    surface = root / 'surface/initial_surface.json'
    surface.write_text('{}')
    with pytest.raises(MicrosolvationStateError, match='surface.*modified'):
        MicrosolvationRunState.load(root, request().to_state_dict())


def test_preflight_on_completed_run_is_read_only(tmp_path, monkeypatch):
    from pyar.scripts import microsolvate
    from pyar.cli import main
    monkeypatch.chdir(tmp_path)
    req = request()
    req.solute.mol_to_xyz('solute.xyz')
    req.solvent.mol_to_xyz('solvent.xyz')
    vector = ['solute.xyz','solvent.xyz','--count','2','--orientations','3',
              '--maximum-number-of-seeds','2','--surface-points-per-atom','32','--shell-tolerance','3.5']
    args = microsolvate.build_parser().parse_args(vector)
    workflow.microsolvate(microsolvate.resolve_request(args))
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    main(['microsolvate', *vector, '--check'])
    assert before == {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}


def test_final_provenance_tamper_fails_closed(tmp_path):
    import json
    root = tmp_path / 'run'
    workflow.microsolvate(request(), output=root)
    path = root / 'state.json'
    data = read_json(path)
    data['final_seeds'][0]['rank'] = 99
    path.write_text(json.dumps(data))
    with pytest.raises(MicrosolvationStateError, match='Final result'):
        MicrosolvationRunState.load(root, request().to_state_dict())


def test_uncommitted_stage_path_escape_is_rejected(tmp_path):
    req = request()
    root = tmp_path / 'run'
    initial = req.solute.copy()
    initial.solute_atom_indices = tuple(range(len(initial)))
    initial.solvent_fragments = ()
    initial.microsolvation_coverage = {}
    MicrosolvationRunStateV2.create(root, req.to_state_dict(), initial, req.solvent)
    outside = tmp_path / 'outside'
    outside.mkdir()
    (root / 'steps/step_001').symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='escapes'):
        workflow.microsolvate(req, output=root)
    assert not list(outside.iterdir())
