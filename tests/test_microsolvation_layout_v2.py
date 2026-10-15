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
from pyar.state.microsolvation import MicrosolvationRunState, MicrosolvationStateError

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
    old = workflow._microsolvate_v1(request(), output=tmp_path / 'old')
    new = workflow.microsolvate(request(), output=tmp_path / 'new')
    assert old.status == new.status == 'completed'
    assert [Path(p).read_bytes() for p in old.selected_paths] == [Path(p).read_bytes() for p in new.selected_paths]
    assert old.metadata['coverage_history'] == new.metadata['coverage_history']
    assert old.metadata['solute_atom_indices'] == new.metadata['solute_atom_indices']
    assert old.metadata['solvent_fragments'] == new.metadata['solvent_fragments']
    before = {p: p.read_bytes() for p in (tmp_path / 'new').rglob('*') if p.is_file()}
    workflow.microsolvate(request(), output=tmp_path / 'new')
    assert before == {p: p.read_bytes() for p in (tmp_path / 'new').rglob('*') if p.is_file()}


def test_durable_legacy_fixture_resumes_without_conversion(tmp_path):
    root = tmp_path / 'old'
    shutil.copytree(FIXTURE, root)
    state = MicrosolvationRunState.load(root, request().to_state_dict())
    assert state.data['next_step'] == 2
    before = {p: p.read_bytes() for p in (root / 'step_001/selected').glob('*.xyz')}
    calls = []
    generate = workflow.generate_microsolvation_candidates
    def record(*args, **kwargs):
        calls.append(kwargs['step'])
        return generate(*args, **kwargs)
    with patch.object(workflow, 'generate_microsolvation_candidates', side_effect=record):
        result = workflow.microsolvate(request(), output=root)
    assert calls and all(step >= 20000 for step in calls)
    assert result.status == 'completed' and not (root / 'layout.json').exists()
    assert before == {p: p.read_bytes() for p in (root / 'step_001/selected').glob('*.xyz')}


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
    from pyar.scripts import modern_microsolvate
    from pyar.modern_cli import main
    monkeypatch.chdir(tmp_path)
    req = request()
    req.solute.mol_to_xyz('solute.xyz')
    req.solvent.mol_to_xyz('solvent.xyz')
    vector = ['solute.xyz','solvent.xyz','--count','2','--orientations','3',
              '--maximum-number-of-seeds','2','--surface-points-per-atom','32','--shell-tolerance','3.5']
    args = modern_microsolvate.build_parser().parse_args(vector)
    workflow.microsolvate(modern_microsolvate.resolve_request(args))
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
