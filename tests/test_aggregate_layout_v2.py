"""Keep final cross-path selection ordering and legacy aggregate continuation."""
import importlib
import os
from pathlib import Path
import random
import shutil
from unittest.mock import patch

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.persistence import read_json, sha256_file
from pyar.aggregation.layout_v2 import AggregateRunStateV2
from pyar.state.aggregate import AggregateStateError

workflow = importlib.import_module('pyar.workflows.aggregate')
FIXTURE = Path(__file__).parent / 'data/run_layout_v1/aggregate_interrupted'


def inputs():
    return [Molecule(['H','H'], [[0,0,0],[0,0,.74]], name='A'),
            Molecule(['He'], [[0,0,0]], name='B')]


def run(function, directory, **extra):
    directory.mkdir()
    with workflow._working_directory(directory):
        random.seed(21); np.random.seed(21)
        return function(inputs(), [1,2], 4, {}, 2, 0, 2, None,
                        connectivity_policy='off', selection_feature='distance-histogram', **extra)


def test_real_science_preserves_source_order_names_groups_and_geometry(tmp_path, monkeypatch):
    # The legacy selector returns a set of tuples containing identity-hashed
    # Molecules. Separately allocated inputs can therefore visit the same
    # pathways in different orders despite identical RNG seeds. Compare both
    # layouts against the same ordered pathway plan without changing sampling.
    select_pathways = workflow.select_pathways
    monkeypatch.setattr(workflow, 'select_pathways', lambda *args:
                        sorted(select_pathways(*args), key=workflow._format_path))
    old = run(workflow._aggregate_v1, tmp_path / 'old')
    new = run(workflow.aggregate, tmp_path / 'new')
    assert old.status == new.status == 'completed'
    assert old.metadata['pathways'] == new.metadata['pathways']
    old_bytes = [Path(old.run_directory, path).read_bytes() for path in old.selected_paths]
    new_bytes = [Path(new.run_directory, path).read_bytes() for path in new.selected_paths]
    assert old_bytes == new_bytes
    manifest = read_json(Path(new.run_directory) / 'results/manifest.json')
    assert len(manifest['artifacts']) == len(old.selected_paths)
    for artifact, raw in zip(manifest['artifacts'], old_bytes):
        assert (Path(new.run_directory) / artifact['path']).read_bytes() == raw
    assert not (Path(old.run_directory) / 'layout.json').exists()
    with workflow._working_directory(tmp_path / 'new'):
        with pytest.raises(AggregateStateError, match='already'):
            workflow.aggregate(inputs(), [1,2], 4, {}, 2, 0, 2, None,
                connectivity_policy='off', selection_feature='distance-histogram')


def test_interrupted_path_reuses_verified_addition(tmp_path, monkeypatch):
    root = tmp_path / 'run'
    root.mkdir()
    original = workflow.add_one
    calls = []
    def interrupt(*args, **kwargs):
        calls.append(args[0])
        if len(calls) == 2:
            raise RuntimeError('interrupt second addition')
        return original(*args, **kwargs)
    with workflow._working_directory(root):
        monkeypatch.setattr(workflow, 'add_one', interrupt)
        with pytest.raises(RuntimeError, match='second addition'):
            workflow.aggregate(inputs(), [1,2], 4, {}, 2, 0, 2, None,
                connectivity_policy='off', selection_feature='distance-histogram')
        committed = list((root / 'aggregates/pathways').glob('pathway_*/steps/step_*/selection.json'))
        hashes = {path: sha256_file(path) for path in committed}
        monkeypatch.setattr(workflow, 'add_one', original)
        result = workflow.aggregate(inputs(), [1,2], 4, {}, 2, 0, 2, None,
            connectivity_policy='off', selection_feature='distance-histogram')
    assert result.status == 'completed'
    assert hashes == {path: sha256_file(path) for path in hashes}


def test_legacy_fixture_resumes_original_state_without_conversion(tmp_path):
    root = tmp_path / 'old'
    root.mkdir()
    shutil.copytree(FIXTURE, root / 'aggregates')
    molecules = [Molecule(['H'], [[0,0,0]], name=name) for name in ('input_a','input_b')]
    for molecule in molecules:
        molecule.fragments = [0]
    with workflow._working_directory(root), patch.object(workflow, 'add_one', return_value=[]) as engine:
        result = workflow.aggregate(molecules, [1,1], 2, {'software':None}, 2, 0, 2, None)
    assert result.status == 'completed' and engine.call_count == 1
    assert not (root / 'aggregates/layout.json').exists()


def test_request_and_selected_snapshot_tamper_fail_closed(tmp_path):
    result = run(workflow.aggregate, tmp_path / 'run')
    root = Path(result.run_directory)
    state = read_json(root / 'state.json')
    snapshot = next(root.glob('pathways/pathway_*/steps/step_001/selected/selected_*.xyz'))
    snapshot.write_text('tampered')
    # Model the valid crash boundary before final progress publication.
    # This exercises checksum validation rather than completed-run rejection.
    import json
    state['status'] = 'running'
    (root / 'state.json').write_text(json.dumps(state))
    with pytest.raises(AggregateStateError, match='checksum'):
        AggregateRunStateV2.load(root.parent, state['request'])


def test_final_commit_can_be_adopted_after_progress_crash_without_reselection(tmp_path, monkeypatch):
    root = tmp_path / 'run'
    root.mkdir()
    save = AggregateRunStateV2.save
    failed = []
    def crash(state):
        if state.data['status'] == 'completed' and not failed:
            failed.append(True)
            raise OSError('crash after final selection commit')
        save(state)
    with workflow._working_directory(root):
        monkeypatch.setattr(AggregateRunStateV2, 'save', crash)
        with pytest.raises(OSError, match='final selection'):
            workflow.aggregate(inputs(), [1,2], 4, {}, 2, 0, 2, None,
                connectivity_policy='off', selection_feature='distance-histogram')
        hashes = {p: sha256_file(p) for p in (root / 'aggregates').glob('pathways/pathway_*/steps/step_*/selection.json')}
        monkeypatch.setattr(AggregateRunStateV2, 'save', save)
        with patch.object(workflow, 'add_one') as add, patch.object(workflow, '_finalize_selected_geometries') as select:
            result = workflow.aggregate(inputs(), [1,2], 4, {}, 2, 0, 2, None,
                connectivity_policy='off', selection_feature='distance-histogram')
            add.assert_not_called()
            select.assert_not_called()
    assert result.status == 'completed'
    assert hashes == {p: sha256_file(p) for p in hashes}
