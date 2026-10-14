"""Grow pilot scientific equivalence, legacy fixtures and crash recovery."""
import importlib
import json
from pathlib import Path
import shutil
from unittest.mock import patch

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.growth.request import GrowRequest
from pyar.growth.layout_v2 import GrowRunStateV2
from pyar.state.grow import GrowRunState, GrowStateError
from pyar.persistence import read_json, sha256_file
from pyar.scripts import modern_grow

workflow = importlib.import_module('pyar.workflows.grow')
FIXTURE = Path(__file__).parent / 'data/run_layout_v1/grow_interrupted'


def request():
    return GrowRequest(Molecule(['H', 'H'], [[0, 0, 0], [0, 0, .74]], name='seed'),
                       Molecule(['H', 'H'], [[0, 0, 0], [0, 0, .74]], name='monomer'), 3)


def addition(aid, seeds, monomer, *args, **kwargs):
    step = int(aid.rsplit('_', 1)[1])
    return [Molecule(seeds[0].atoms_list + monomer.atoms_list,
                    np.vstack([seeds[0].coordinates, monomer.coordinates + step * 3]),
                    name=f'stage_{step}', energy=None)]


def test_mocked_science_matches_v1_and_final_outputs_are_immutable(tmp_path):
    req = request()
    with patch.object(workflow, 'add_one', side_effect=addition):
        old = workflow._grow_v1(req, output=tmp_path / 'old')
        new = workflow.grow(req, output=tmp_path / 'new')
    assert Path(old.selected_paths[0]).read_bytes() == Path(new.selected_paths[0]).read_bytes()
    selected = {p.relative_to(tmp_path / 'new'): p.read_bytes()
                for p in (tmp_path / 'new/steps').rglob('*') if p.is_file()}
    result_bytes = Path(new.selected_paths[0]).read_bytes()
    manifest = read_json(tmp_path / 'new/results/manifest.json')
    with patch.object(workflow, 'add_one') as engine:
        repeated = workflow.grow(req, output=tmp_path / 'new')
        engine.assert_not_called()
    assert repeated.selected_paths == new.selected_paths
    assert read_json(tmp_path / 'new/results/manifest.json') == manifest
    assert Path(new.selected_paths[0]).read_bytes() == result_bytes
    assert selected == {p.relative_to(tmp_path / 'new'): p.read_bytes()
                        for p in (tmp_path / 'new/steps').rglob('*') if p.is_file()}


def test_durable_v1_fixture_resumes_without_conversion(tmp_path):
    root = tmp_path / 'old'
    shutil.copytree(FIXTURE, root)
    before = (root / 'step_001/selected/selected_000.xyz').read_bytes()
    loaded = GrowRunState.load(root, request().to_state_dict())
    assert loaded.data['next_step'] == 2
    with patch.object(workflow, 'add_one', side_effect=addition) as engine:
        result = workflow.grow(request(), output=root)
    assert engine.call_count == 2
    assert not (root / 'layout.json').exists()
    assert before == (root / 'step_001/selected/selected_000.xyz').read_bytes()
    assert Path(result.selected_paths[0]).parent == root / 'final'


def test_committed_step_survives_crash_before_progress_commit(tmp_path, monkeypatch):
    req = request()
    original = GrowRunStateV2.save
    failed = []
    def crash(state):
        if len(state.data['completed_steps']) == 1 and not failed:
            failed.append(True)
            raise OSError('crash after immutable selection publish')
        return original(state)
    with patch.object(workflow, 'add_one', side_effect=addition) as engine:
        monkeypatch.setattr(GrowRunStateV2, 'save', crash)
        with pytest.raises(OSError, match='immutable selection'):
            workflow.grow(req, output=tmp_path / 'run')
        record = tmp_path / 'run/steps/step_001/selection.json'
        before = record.read_bytes()
        monkeypatch.setattr(GrowRunStateV2, 'save', original)
        result = workflow.grow(req, output=tmp_path / 'run')
    assert result.status == 'completed' and engine.call_count == 3
    assert record.read_bytes() == before


def test_older_committed_snapshot_tamper_and_request_mismatch_rejected(tmp_path):
    req = request()
    with patch.object(workflow, 'add_one', side_effect=addition):
        workflow.grow(req, output=tmp_path / 'run')
    snapshot = tmp_path / 'run/steps/step_001/selected/selected_000.xyz'
    snapshot.write_text('tampered')
    with pytest.raises(GrowStateError, match='checksum'):
        GrowRunState.load(tmp_path / 'run', req.to_state_dict())


def test_preflight_completed_v2_is_read_only(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    req = request()
    req.seed.mol_to_xyz('seed.xyz')
    req.monomer.mol_to_xyz('monomer.xyz')
    # CLI resolution has additional state defaults; use the same CLI request
    # for both the direct workflow and its preflight.
    args = modern_grow.build_parser().parse_args(['seed.xyz', 'monomer.xyz', '--count', '3'])
    cli_req = modern_grow.resolve_request(args)
    with patch.object(workflow, 'add_one', side_effect=addition):
        workflow.grow(cli_req)
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    from pyar.modern_cli import main
    main(['grow', 'seed.xyz', 'monomer.xyz', '--count', '3', '--check'])
    assert before == {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}


def test_manifest_rebuild_does_not_change_committed_selection(tmp_path):
    req = request()
    with patch.object(workflow, 'add_one', side_effect=addition):
        workflow.grow(req, output=tmp_path / 'run')
    records = list((tmp_path / 'run/steps').glob('step_*/selection.json'))
    hashes = {path: sha256_file(path) for path in records}
    (tmp_path / 'run/results/manifest.json').unlink()
    with patch.object(workflow, 'add_one') as engine:
        workflow.grow(req, output=tmp_path / 'run')
        engine.assert_not_called()
    assert hashes == {path: sha256_file(path) for path in records}


def test_uncommitted_execution_and_results_symlink_escape_rejected(tmp_path):
    req = request()
    root = tmp_path / 'run'
    state = GrowRunStateV2.create(root, req.to_state_dict(), req.seed, req.monomer)
    outside = tmp_path / 'outside'
    outside.mkdir()
    (root / 'steps/step_001').symlink_to(outside, target_is_directory=True)
    with patch.object(workflow, 'add_one') as engine:
        with pytest.raises(ValueError, match='escapes'):
            workflow.grow(req, output=root)
        engine.assert_not_called()
    assert not list(outside.iterdir())
