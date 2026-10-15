"""Scan commit reuse across continuation controls, independent of progress reports."""
import importlib
from pathlib import Path

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.persistence import read_json, sha256_file

workflow = importlib.import_module('pyar.workflows.scan_bond')


@pytest.fixture
def scan_run(tmp_path, monkeypatch):
    a, b = tmp_path / 'a.xyz', tmp_path / 'b.xyz'
    for path in (a, b):
        path.write_text('1\ninput\nH 0 0 0\n')
    molecule = Molecule(['H', 'H'], [[0, 0, 0], [2, 0, 0]])
    monkeypatch.setattr(workflow, '_make_orientations', lambda *args: [(0, molecule, None)])
    monkeypatch.setattr(workflow, 'separated_reactant_identity', lambda *args: None)
    monkeypatch.setattr('pyar.backends.orca_scan.require_executable', lambda *args: 'orca')
    calls = []
    def execute(command, stdout_path=None, stderr_path=None):
        calls.append(command)
        Path(stdout_path).write_text('****ORCA TERMINATED NORMALLY****\n')
        Path('scan.allxyz').write_text('2\nfirst\nH 0 0 0\nH 2 0 0\n2\nlast\nH 0 0 0\nH 1 0 0\n')
        Path('scan.relaxscanact.dat').write_text('2.0 -1.0\n1.0 -0.5\n')
        return 0
    monkeypatch.setattr('pyar.backends.orca_scan.run_command', execute)
    root = tmp_path / 'scan'
    def run(through='scan', legacy=False, **options):
        function = workflow._run_scan_bond if legacy else workflow.run_scan_bond
        target = root.with_name('legacy') if legacy else root
        return function(a, b, (0, 0), 1, {'software':'orca', 'method':'BP86', 'basis':'def2-SVP'},
                                      target, scan_end=1., scan_points=2, through=through, **options)
    return root, run, calls


def test_scan_extension_preserves_commit_and_native_bytes(scan_run, monkeypatch):
    root, run, calls = scan_run
    assert run()['status'] == 'complete'
    commit = root / 'orientation_000/scan_commit.json'
    record = read_json(commit)
    hashes = {root / item['path']: (item['sha256'], (root / item['path']).stat().st_mtime_ns)
              for item in record['artifacts']}
    digest = sha256_file(commit)
    def continuation(root, results, *args):
        for result in results:
            result['continuation_status'] = 'complete'
    monkeypatch.setattr('pyar.workflows.scan_path.continue_scan_paths', continuation)
    (root / 'orientation_000/result.json').unlink()
    assert run('all')['status'] == 'complete'
    assert len(calls) == 1
    assert sha256_file(commit) == digest
    assert all((sha256_file(path), path.stat().st_mtime_ns) == evidence
               for path, evidence in hashes.items())
    assert all(not Path(ref['path']).is_absolute() for ref in record['artifacts'])
    assert read_json(root / 'layout.json')['layout_version'] == 2
    assert read_json(root / 'results/manifest.json')['artifacts']


def test_committed_scan_tampering_fails_without_recalculation(scan_run):
    root, run, calls = scan_run
    run()
    record = read_json(root / 'orientation_000/scan_commit.json')
    path = root / next(ref['path'] for ref in record['artifacts'] if ref['path'].endswith('scan_trajectory.xyz'))
    path.write_text('tampered')
    with pytest.raises(ValueError, match='checksum'):
        run()
    assert len(calls) == 1


def test_changed_scan_request_rejected_before_jobs(scan_run):
    root, run, calls = scan_run
    run()
    (root.parent / 'a.xyz').write_text('1\nchanged\nHe 0 0 0\n')
    with pytest.raises(FileExistsError, match='different request'):
        run()
    assert len(calls) == 1


def test_v2_stage_commits_index_cumulative_native_executor(tmp_path, monkeypatch):
    import pyar.neb as neb
    from pyar.workflows.scan_path import run_path_stages
    inputs = []
    for name in ('start', 'end', 'guess'):
        path = tmp_path / f'{name}.xyz'
        path.write_text('1\ninput\nH 0 0 0\n')
        inputs.append(path)
    calls = []
    def native(**kwargs):
        from pyar.persistence import write_json
        calls.append(kwargs['stage'])
        output = Path(kwargs['output'])
        output.mkdir(parents=True, exist_ok=True)
        result = dict(stage=kwargs['stage'], status='complete', inputs={}, artifacts={}, dependencies=[],
                      **{key: True for key in neb._STAGE_GATES[kwargs['stage']]})
        result['interior_maximum'] = True
        write_json(output / f"{kwargs['stage']}_summary.json", result)
        return result
    monkeypatch.setattr(neb, 'run_neb', native)
    result = run_path_stages(start=inputs[0], end=inputs[1], ts_guess=inputs[2],
        output=tmp_path / 'path', qc_params={'software':'xtb'}, through='neb', options={},
        layout_dependency='a' * 64)
    assert result['status'] == 'complete'
    assert calls == ['relax', 'neb']
    for stage in calls:
        record = read_json(tmp_path / f'path/stages/{stage}/attempts/attempt_000/commit.json')
        assert record['scientific_validation']
        assert all(not Path(ref['path']).is_absolute() for ref in record['artifacts'])


def test_native_scan_outputs_match_legacy_and_legacy_routing_is_preserved(scan_run):
    root, run, calls = scan_run
    old = run(legacy=True)
    new = run()
    for name in ('final_scan_path', 'trajectory_path', 'scan_profile_csv', 'scan_profile_json'):
        assert Path(old['results'][0][name]).read_bytes() == Path(new['results'][0][name]).read_bytes()
    args = read_json(root.with_name('legacy') / 'request.json')
    assert 'request_schema_version' not in args
    assert not (root.with_name('legacy') / 'layout.json').exists()
    # Feed the persisted v1 request to the public router without converting it.
    a, b = root.parent / 'a.xyz', root.parent / 'b.xyz'
    result = workflow.run_scan_bond(a, b, (0, 0), 1, args['qc_params'], root.with_name('legacy'),
                                    scan_end=1., scan_points=2, through='scan')
    assert result['status'] == 'complete'
    assert not (root.with_name('legacy') / 'layout.json').exists()
    assert len(calls) == 2


def test_free_relaxation_commit_survives_progress_loss(scan_run, monkeypatch):
    from pyar.backends import write_xyz
    root, run, calls = scan_run
    optimizations = []
    def optimize(molecule, parameters):
        optimizations.append(parameters)
        destination = Path('job_relaxed/result_relaxed.xyz')
        destination.parent.mkdir()
        write_xyz(molecule.atoms_list, molecule.coordinates, destination)
        return True
    monkeypatch.setattr(workflow, 'optimise', optimize)
    monkeypatch.setattr(workflow, 'molecule_identity_from_xyz', lambda *args: {})
    assert run(None)['status'] == 'complete'
    record = root / 'orientation_000/relax_commit.json'
    digest = sha256_file(record)
    (root / 'orientation_000/result.json').unlink()
    assert run(None)['status'] == 'complete'
    assert len(calls) == len(optimizations) == 1
    assert sha256_file(record) == digest
