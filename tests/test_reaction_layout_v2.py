"""Immutable discovery boundaries and read-only completed evidence."""
from pathlib import Path

import pytest

from pyar.core.molecule import Molecule
from pyar.persistence import read_json, write_json, sha256_file
from pyar.state.reaction import ReactionRunState, ReactionStateError
from pyar.state.reaction_v2 import ReactionRunStateV2


def create(root):
    molecules = [Molecule(['H'], [[0, 0, 0]], name=name) for name in ('a', 'b')]
    request = dict(gamma_schedule=[100., 200.], backend_parameters={'software':'xtb'}, orientations=1)
    return ReactionRunStateV2.create(root, request, [molecules[0]], molecules), request


def test_geometry_hash_and_request_are_checked_before_resume(tmp_path):
    state, request = create(tmp_path)
    assert isinstance(ReactionRunState.load(tmp_path, request), ReactionRunStateV2)
    with pytest.raises(ReactionStateError, match='does not match'):
        ReactionRunState.load(tmp_path, dict(request, orientations=2))
    snapshot = state.reaction_directory / state.data['pending_orientations'][0]['path']
    snapshot.write_text('tampered')
    with pytest.raises(ReactionStateError, match='checksum'):
        ReactionRunState.load(tmp_path, request)


def test_committed_cycle_adopted_after_progress_publication_failure(tmp_path, monkeypatch):
    state, request = create(tmp_path)
    execution = state.cycle_execution('0100', 100.)
    execution.mkdir(parents=True)
    (execution / 'native.txt').write_text('native cache')
    def crash():
        raise OSError('progress crash')
    monkeypatch.setattr(state, 'save', crash)
    with pytest.raises(OSError, match='progress crash'):
        state.complete_cycle(100., [])
    digest = sha256_file(state.cycle_path(0))
    loaded = ReactionRunState.load(tmp_path, request)
    assert loaded.data['current_cycle'] == 1
    assert loaded.remaining_gamma_schedule() == [200.]
    assert sha256_file(state.cycle_path(0)) == digest


def test_completed_discovery_ignores_progress_reports_and_rejects_evidence_changes(tmp_path):
    state, request = create(tmp_path)
    execution = state.cycle_execution('0100', 100.)
    execution.mkdir(parents=True)
    native = execution / 'native.txt'
    native.write_text('accepted evidence')
    state.complete_cycle(100., [])
    state.finish('completed_no_candidates')
    commit = state.reaction_directory / 'discovery/commit.json'
    digest = sha256_file(commit)
    report = read_json(state.state_file)
    report['inspection_note'] = 'regenerated report'
    report['status'] = 'running'
    write_json(state.state_file, report)
    loaded = ReactionRunState.read_completed(tmp_path, request)
    assert loaded.data['status'] == 'completed_no_candidates'
    assert sha256_file(commit) == digest
    with pytest.raises(ReactionStateError, match='completed'):
        ReactionRunState.load(tmp_path, request)
    native.write_text('tampered')
    with pytest.raises(ReactionStateError, match='checksum'):
        ReactionRunState.read_completed(tmp_path, request)


def test_accepted_product_and_release_snapshot_recover_before_job_progress(tmp_path, monkeypatch):
    state, request = create(tmp_path)
    execution = state.cycle_execution('0100', 100.)
    candidates = execution / 'orientation_x/job_0100_x/candidate_ts'
    candidates.mkdir(parents=True)
    (candidates / 'metadata.json').write_text('{}')
    source = execution / 'orientation_x/result_relax.xyz'
    source.write_text('2\naccepted\nH 0 0 0\nH 0 0 .74\n')
    probe = execution / 'orientation_x/release_probe.json'
    probe.write_text('{"release_survived": true}')
    def crash():
        raise OSError('job progress crash')
    monkeypatch.setattr(state, 'save', crash)
    with pytest.raises(OSError, match='job progress crash'):
        state.record_product('0100_x', 100., 'h2', '[HH]', source,
                             trace_summary={'candidate_ts_directory': str(candidates)})
    # A subsequent native retry cannot alter the already accepted evidence.
    probe.write_text('overwritten native scratch')
    loaded = ReactionRunState.load(tmp_path, request)
    product = loaded.data['products'][0]
    assert product['job_name'] == '0100_x'
    recorded = next(ref for ref in product['evidence'] if ref['path'].endswith('release_probe.json'))
    assert read_json(loaded.reaction_directory / recorded['path'])['release_survived']
