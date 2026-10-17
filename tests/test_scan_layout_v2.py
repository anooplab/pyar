"""Scan commit reuse across continuation controls, independent of progress reports."""
import importlib
from pathlib import Path

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.persistence import read_json, sha256_file, write_json

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
    import os
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
    flushed = []
    original_fsync = os.fsync
    def record_fsync(descriptor):
        flushed.append(os.readlink(f'/proc/self/fd/{descriptor}'))
        original_fsync(descriptor)
    monkeypatch.setattr(os, 'fsync', record_fsync)
    result = run_path_stages(start=inputs[0], end=inputs[1], ts_guess=inputs[2],
        output=tmp_path / 'path', qc_params={'software':'xtb'}, through='neb', options={},
        layout_dependency='a' * 64)
    assert result['status'] == 'complete', result
    assert calls == ['relax', 'neb']
    for stage in calls:
        record = read_json(tmp_path / f'path/stages/{stage}/attempts/attempt_000/commit.json')
        assert record['scientific_validation']
        assert all(not Path(ref['path']).is_absolute() for ref in record['artifacts'])
        summary_index = next(index for index, path in enumerate(flushed)
                             if f'/execution/{stage}_summary.json' in path)
        commit_index = next(index for index, path in enumerate(flushed)
                            if f'/stages/{stage}/attempts/attempt_000/.pyar-json-' in path)
        assert summary_index < commit_index


def test_scan_path_seed_is_forwarded_only_to_neb_stage(tmp_path, monkeypatch):
    import pyar.neb as neb
    from pyar.workflows.scan_path import run_path_stages
    endpoints = []
    for name in ('start', 'end', 'seed'):
        path = tmp_path / f'{name}.xyz'
        path.write_text('2\nframe\nH 0 0 0\nH 1 0 0\n')
        endpoints.append(path)
    calls = []
    def native(**kwargs):
        calls.append(kwargs)
        result = dict(stage=kwargs['stage'], status='complete', inputs={}, artifacts={}, dependencies=[],
                      **{key: True for key in neb._STAGE_GATES[kwargs['stage']]})
        result['interior_maximum'] = True
        from pyar.persistence import write_json
        output = Path(kwargs['output'])
        output.mkdir(parents=True, exist_ok=True)
        write_json(output / f"{kwargs['stage']}_summary.json", result)
        return result
    monkeypatch.setattr(neb, 'run_neb', native)
    result = run_path_stages(start=endpoints[0], end=endpoints[1], path_seed=endpoints[2],
        output=tmp_path / 'path', qc_params={'software':'xtb'}, through='neb',
        options={'path_refinement':'none', 'path_seed_maximum_frame_index':1},
        layout_dependency='b' * 64)
    assert result['status'] == 'complete'
    assert [call['stage'] for call in calls] == ['relax', 'neb']
    assert calls[0].get('path_seed') is None
    assert 'path_seed_maximum_frame_index' not in calls[0]
    assert calls[1]['path_seed'].name == 'path_seed.xyz'
    assert calls[1]['path_seed_maximum_frame_index'] == 1
    assert (tmp_path / 'path/inputs/path_seed.xyz').is_file()


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


def test_scan_commit_recovers_and_rebuilds_a_damaged_result_report(scan_run):
    root, run, calls = scan_run
    assert run()['status'] == 'complete'
    commit = root / 'orientation_000/scan_commit.json'
    digest = sha256_file(commit)
    report = root / 'orientation_000/result.json'
    report.write_text('{"status": "running"}\n')
    assert run()['status'] == 'complete'
    assert len(calls) == 1
    assert sha256_file(commit) == digest
    rebuilt = read_json(report)
    assert rebuilt['scan_status'] == 'success'


def test_scan_rejects_symlinked_execution_before_backend_writes(scan_run, monkeypatch):
    root, run, calls = scan_run
    monkeypatch.setattr(workflow, 'run_orca_bond_scan',
                        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError('backend stopped')))
    assert run()['status'] == 'failed'
    execution = root / 'orientation_000/attempts/attempt_000/execution'
    execution.rename(execution.with_name('failed_execution'))
    outside = root.parent / 'outside_execution'
    outside.mkdir()
    execution.symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='escapes'):
        run()
    assert list(outside.iterdir()) == []
    assert not calls


def test_scan_creates_nested_output_parent(scan_run):
    root, run, calls = scan_run
    nested = root.parent / 'new_parent' / 'scan'
    result = workflow.run_scan_bond(
        root.parent / 'a.xyz', root.parent / 'b.xyz', (0, 0), 1,
        {'software':'orca', 'method':'BP86', 'basis':'def2-SVP'}, nested,
        scan_end=1., scan_points=2, through='scan')
    assert result['status'] == 'complete'
    assert nested.is_dir()
    assert len(calls) == 1


def test_scan_native_evidence_is_fsynced_before_its_commit(scan_run, monkeypatch):
    import os
    root, run, _calls = scan_run
    targets = []
    original = os.fsync
    def capture(descriptor):
        targets.append(os.readlink(f'/proc/self/fd/{descriptor}'))
        original(descriptor)
    monkeypatch.setattr(os, 'fsync', capture)
    assert run()['status'] == 'complete'
    assert any('/orientation_000/attempts/attempt_000/execution/scan/scan.allxyz' in path
               for path in targets)
    evidence_index = next(index for index, path in enumerate(targets)
                          if '/orientation_000/attempts/attempt_000/execution/scan/scan.allxyz' in path)
    commit_index = next(index for index, path in enumerate(targets)
                        if '/orientation_000/.pyar-json-' in path)
    assert evidence_index < commit_index


def test_v2_scan_products_are_grouped_deterministically_without_path_execution(tmp_path, monkeypatch):
    from pyar.backends import write_xyz
    from pyar.backends.orca_scan import OrcaBondScanResult

    fragment_a = tmp_path / 'a.xyz'
    fragment_b = tmp_path / 'b.xyz'
    write_xyz(['H'], [[0., 0., 0.]], fragment_a)
    write_xyz(['H'], [[1.4, 0., 0.]], fragment_b)
    workflow_module = importlib.import_module('pyar.workflows.scan_bond')
    orientations = []
    for index in range(5):
        molecule = Molecule(['H', 'H'], [[0., 0., 0.], [1.4, 0., 0.]], name=f'orientation_{index:03d}')
        molecule.fragments = [[0], [1]]
        orientations.append((index, molecule, None))
    monkeypatch.setattr(workflow_module, '_make_orientations', lambda *args: orientations)
    monkeypatch.setattr(workflow_module, 'separated_reactant_identity',
                        lambda *args: {'inchi': 'reactants', 'smiles': '[H].[H]'})
    scan_calls = []

    def fake_scan(molecule, request, directory, qc_params):
        scan_calls.append(molecule.name)
        directory.mkdir(parents=True, exist_ok=True)
        index = int(molecule.name.rsplit('_', 1)[1])
        distances = [1.4, 1.2, 1.0]
        frames = [(molecule.atoms_list, np.array([[0., 0., 0.], [distance, 0., 0.]]), f'frame {i}')
                  for i, distance in enumerate(distances)]
        energies = {
            0: [-1.0, -0.5, -0.8],
            1: [-1.0, -0.8, -0.95],
            2: [-1.0, -0.7, -0.9],
            3: [-1.0, -0.6, -0.9],
            4: [-1.0, -0.4, -0.9],
        }[index]
        trajectory = directory / 'scan_trajectory.xyz'
        final = directory / 'final_scan.xyz'
        scan_input = directory / 'scan.inp'
        scan_output = directory / 'scan.out'
        trajectory.write_text(''.join(
            f'2\nframe {i}\nH 0 0 0\nH {distance} 0 0\n'
            for i, distance in enumerate(distances)
        ))
        write_xyz(molecule.atoms_list, frames[-1][1], final)
        scan_input.write_text('scan request\n')
        scan_output.write_text('scan completed\n')
        profile = [{'target_distance_angstrom': distance, 'energy_hartree': energy}
                   for distance, energy in zip(distances, energies)]
        return OrcaBondScanResult(True, trajectory, final, frames[-1][1], scan_input,
                                  scan_output, 'success', frames=frames, profile=profile)

    def fake_optimize(molecule, qc_params):
        orientation = next(int(part.rsplit('_', 1)[1]) for part in Path.cwd().parts
                           if part.startswith('orientation_'))
        if orientation == 2:
            molecule.coordinates[1, 0] = 2.0  # A changed but off-target product.
        output = Path('job_relaxed/result_relaxed.xyz')
        output.parent.mkdir(parents=True)
        write_xyz(molecule.atoms_list, molecule.coordinates, output)
        return True

    def fake_identity(path):
        orientation = next(int(part.rsplit('_', 1)[1]) for part in Path(path).parts
                            if part.startswith('orientation_'))
        identities = {
            0: {'inchi': 'InChI=product-one', 'smiles': 'C1'},
            1: {'inchi': 'InChI=product-one', 'smiles': 'C1 alternate'},
            2: {'inchi': 'InChI=product-two', 'smiles': 'C2'},
            4: {'inchi': 'reactants', 'smiles': '[H].[H]'},
        }
        if orientation == 3:
            raise ValueError('identity provider failed')
        return identities[orientation]

    monkeypatch.setattr(workflow_module, 'run_orca_bond_scan', fake_scan)
    monkeypatch.setattr(workflow_module, 'optimise', fake_optimize)
    monkeypatch.setattr(workflow_module, 'molecule_identity_from_xyz', fake_identity)
    path_module = importlib.import_module('pyar.workflows.scan_path')
    monkeypatch.setattr(path_module, 'continue_scan_paths',
                        lambda *args: pytest.fail('M5 must not start NEB/path stages'))
    root = tmp_path / 'scan_bond'
    params = {'software': 'orca', 'method': 'BP86', 'basis': 'def2-SVP'}
    result = workflow.run_scan_bond(fragment_a, fragment_b, (0, 0), 5, params, root,
                                    scan_end=1.0, scan_points=3, through=None)
    assert result['status'] == 'complete'
    index = read_json(root / 'products/index.json')
    assert len(index['products']) == 2
    product_one = next(row for row in index['products'] if row['identity']['inchi'] == 'InChI=product-one')
    assert product_one['orientation_ids'] == [0, 1]
    assert product_one['representative_orientation'] == 1
    assert product_one['alternative_orientations'] == [0]
    assert product_one['identity']['smiles_variants'] == ['C1', 'C1 alternate']
    product_two = next(row for row in index['products'] if row['identity']['inchi'] == 'InChI=product-two')
    assert product_two['target_status'] == 'off_target'
    assert index['identity_unknown'][0]['orientation'] == 3
    assert index['non_products'] == [{'orientation': 4, 'reason': 'identity_unchanged'}]
    assert [row['product_role'] for row in result['results'][:2]] == ['duplicate_route', 'representative']
    assert result['results'][3]['product_id'] is None
    assert result['results'][3]['product_role'] == 'identity_unknown'
    assert not list(root.glob('orientation_*/reaction_path'))

    digest = sha256_file(root / 'products/index.json')
    request = read_json(root / 'request.json')
    from pyar.workflows.scan_layout import ScanLayout
    layout = ScanLayout(root)
    layout.product_registry(list(reversed(result['results'])), request)
    assert sha256_file(root / 'products/index.json') == digest

    resumed = workflow.run_scan_bond(fragment_a, fragment_b, (0, 0), 5, params, root,
                                      scan_end=1.0, scan_points=3, through=None)
    assert resumed['status'] == 'complete'
    assert len(scan_calls) == 5  # Resume reuses every committed scan.
    assert sha256_file(root / 'products/index.json') == digest

    path_calls = []
    def characterize(**kwargs):
        path_calls.append(kwargs)
        from pyar.workflows.path_layout import PathLayout
        PathLayout(kwargs['output'], kwargs['start'], kwargs['end'], None,
                   kwargs['qc_params'], kwargs['options'], kwargs['layout_dependency'],
                   path_seed=kwargs['path_seed'])
        if len(path_calls) == 1:
            return {'status': 'scientific_gate_failed', 'failed_stage': 'neb',
                    'completed_stages': ['relax']}
        return {'status': 'complete', 'completed_stages': ['relax', 'neb'],
                'reactant_product_connection_confirmed': True}
    monkeypatch.setattr(path_module, 'run_path_stages', characterize)
    characterized = workflow.run_scan_bond(fragment_a, fragment_b, (0, 0), 5, params, root,
        scan_end=1.0, scan_points=3, through='all', characterize_unique_products=True)
    assert characterized['status'] == 'failed'
    assert len(path_calls) == 2  # One representative for each unique product, not five routes.
    assert len(scan_calls) == 5
    assert all(call['path_seed'].is_file() for call in path_calls)
    product_statuses = {row['characterization_status'] for row in
                        read_json(root / 'products/index.json')['products']}
    assert product_statuses == {'path_scientific_gate_failed', 'reaction_connected'}

    record_path = root / product_one['product_record_path']
    record = read_json(record_path)
    record['orientation_ids'] = [1]
    write_json(record_path, record)
    with pytest.raises(ValueError, match='Immutable scan product record changed'):
        layout.product_registry(result['results'], request)


def test_unique_product_characterization_uses_one_scan_seeded_path_and_reuses_scan(scan_run, monkeypatch):
    from pyar.backends import write_xyz
    root, run, scan_calls = scan_run
    relax_calls = []
    def optimize(molecule, parameters):
        relax_calls.append(1)
        destination = Path('job_relaxed/result_relaxed.xyz')
        destination.parent.mkdir()
        write_xyz(molecule.atoms_list, molecule.coordinates, destination)
        return True
    monkeypatch.setattr(workflow, 'optimise', optimize)
    monkeypatch.setattr(workflow, 'separated_reactant_identity', lambda *args: {'inchi': 'reactants'})
    monkeypatch.setattr(workflow, 'molecule_identity_from_xyz', lambda *args: {'inchi': 'product'})
    monkeypatch.setattr(workflow, 'reaction_product_changed', lambda *args: True)
    calls = []
    def path_stages(**kwargs):
        calls.append(kwargs)
        assert kwargs['path_seed'].is_file()
        assert kwargs.get('ts_guess') is None
        from pyar.workflows.path_layout import PathLayout
        PathLayout(kwargs['output'], kwargs['start'], kwargs['end'], None,
                   kwargs['qc_params'], kwargs['options'], kwargs['layout_dependency'],
                   path_seed=kwargs['path_seed'])
        if len(calls) == 1:
            return {'status': 'complete', 'completed_stages': ['relax', 'neb'],
                    'reactant_product_connection_confirmed': True}
        return {'status': 'scientific_gate_failed', 'failed_stage': 'neb',
                'completed_stages': ['relax']}
    monkeypatch.setattr('pyar.workflows.scan_path.run_path_stages', path_stages)

    result = run('all', characterize_unique_products=True)
    assert result['status'] == 'complete'
    assert len(scan_calls) == len(relax_calls) == len(calls) == 1
    assert calls[0]['through'] == 'all'
    index = read_json(root / 'products/index.json')
    assert len(index['products']) == 1
    product = index['products'][0]
    assert product['characterization_status'] == 'reaction_connected'
    assert (root / product['characterization_path']).is_file()
    path_output = calls[0]['output']

    changed = dict(images=5)
    resumed = run('all', characterize_unique_products=True, reaction_options=changed)
    assert resumed['status'] == 'failed'
    assert len(scan_calls) == len(relax_calls) == 1
    assert len(calls) == 2
    assert calls[1]['output'] == path_output
    updated_index = read_json(root / 'products/index.json')
    assert updated_index['products'][0]['characterization_status'] == 'path_scientific_gate_failed'
    assert (root / product['product_record_path']).is_file()  # Product survives a path gate failure.
    assert (root / 'products' / product['product_id'] / 'characterization' /
            'attempts' / 'attempt_001' / 'request.json').is_file()


def test_unique_product_mode_with_no_product_skips_path_service(scan_run, monkeypatch):
    from pyar.backends import write_xyz
    root, run, scan_calls = scan_run
    def optimize(molecule, parameters):
        destination = Path('job_relaxed/result_relaxed.xyz')
        destination.parent.mkdir()
        write_xyz(molecule.atoms_list, molecule.coordinates, destination)
        return True
    monkeypatch.setattr(workflow, 'optimise', optimize)
    monkeypatch.setattr(workflow, 'separated_reactant_identity', lambda *args: {'inchi': 'same'})
    monkeypatch.setattr(workflow, 'molecule_identity_from_xyz', lambda *args: {'inchi': 'same'})
    monkeypatch.setattr(workflow, 'reaction_product_changed', lambda *args: False)
    monkeypatch.setattr('pyar.workflows.scan_path.run_path_stages',
                        lambda **kwargs: pytest.fail('no products means no path work'))
    result = run('all', characterize_unique_products=True)
    assert result['status'] == 'complete'
    assert len(scan_calls) == 1
    assert read_json(root / 'products/index.json')['products'] == []
