"""Single-molecule scan planning, product gates and immutable continuation."""
import importlib
import json
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

from pyar import modern_cli
from pyar.backends import write_xyz
from pyar.backends.orca_scan import OrcaBondScanResult
from pyar.core.molecule import Molecule
from pyar.persistence import read_json, sha256_file
from pyar.workflows.scan_bond import intramolecular_scan_plan, run_scan_bond


@pytest.mark.parametrize('distance,direction', [(.74, 'break'), (2., 'form')])
def test_automatic_plan_uses_current_distance_and_finite_grid(distance, direction):
    molecule = Molecule(['H', 'H'], [[0., 0., 0.], [distance, 0., 0.]])
    plan = intramolecular_scan_plan(molecule, (0, 1))
    assert plan['scan_direction'] == direction
    assert plan['start_distance_angstrom'] == distance
    assert plan['step_angstrom'] <= .1 + 1e-12
    assert plan['intervals'] == plan['n_points'] - 1
    if direction == 'break':
        assert plan['end_distance_angstrom'] > plan['safe_separation_threshold_angstrom']
    else:
        assert plan['end_distance_angstrom'] < plan['bond_threshold_angstrom']


@pytest.mark.parametrize('atoms,end', [((0, 0), None), ((0, 2), None), ((-1, 1), None),
                                     ((0, 1), .5), ((0, 1), float('nan'))])
def test_bad_intramolecular_plans_are_rejected(atoms, end):
    molecule = Molecule(['H', 'H'], [[0., 0., 0.], [.74, 0., 0.]])
    with pytest.raises(ValueError):
        intramolecular_scan_plan(molecule, atoms, scan_end=end)


def test_explicit_grid_and_end_are_used_without_changing_start():
    molecule = Molecule(['H', 'H'], [[0., 0., 0.], [.74, 0., 0.]])
    plan = intramolecular_scan_plan(molecule, (1, 0), scan_end=3., scan_points=7)
    assert plan['start_distance_angstrom'] == .74
    assert plan['end_distance_angstrom'] == 3.
    assert plan['n_points'] == 7 and plan['endpoint_source'] == 'user_explicit'


@pytest.fixture
def cli_case(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('molecule.xyz').write_text('2\nH2\nH 0 0 0\nH .74 0 0\n')
    monkeypatch.setattr('pyar.cli.shutil.which', lambda name: '/bin/backend')
    return ['scan-bond', 'molecule.xyz', '--atoms', '0', '1', '--backend', 'xtb']


def test_single_input_cli_infers_mode_and_does_not_duplicate_charge(cli_case):
    with mock.patch('pyar.scripts.scan_bond.run_scan_bond') as run:
        modern_cli.main([*cli_case, '--charge', '-2'])
    args = run.call_args.args
    assert args[:4] == ('molecule.xyz', None, [0, 1], 1)
    assert args[4]['charge_a'] == -2
    assert 'charge_b' not in args[4]
    assert args[4]['multiplicity_a'] == 1


def test_single_input_preflight_is_readonly_and_shows_plan(cli_case, capsys):
    before = set(Path('.').iterdir())
    with mock.patch('pyar.scripts.scan_bond.run_scan_bond') as run, mock.patch(
        'subprocess.Popen', side_effect=AssertionError('preflight cannot launch a program')
    ):
        modern_cli.main([*cli_case, '--check'])
    run.assert_not_called()
    assert before == set(Path('.').iterdir())
    output = capsys.readouterr().out
    assert 'Intramolecular break' in output and '0.7400 ->' in output
    assert 'points' in output and 'then free relaxation' in output


def test_single_input_dry_run_records_the_proposed_plan(cli_case, capsys):
    with mock.patch('pyar.scripts.scan_bond.run_scan_bond') as run:
        modern_cli.main([*cli_case, '--dry-run', '--json'])
    value = json.loads(capsys.readouterr().out)
    assert 'Plan intramolecular bond formation or breaking' in value['plan']['stages']
    assert 'Relax the scanned endpoint without constraints' in value['plan']['stages']
    run.assert_not_called()


def test_legacy_single_input_does_not_require_fragment_orientation_count(cli_case):
    from pyar.scripts.scan_bond import build_parser
    args = build_parser().parse_args(['molecule.xyz', '--atoms', '0', '1'])
    assert len(args.inputs) == 1


@pytest.mark.parametrize('options', [
    ['--orientations', '2'], ['--charge', '0', '0'], ['--atoms', '0', '0'],
    ['--atoms', '0', '2'], ['--total-multiplicity', '3'], ['--scan-end', '.5'],
])
def test_single_input_invalid_requests_do_not_start(cli_case, options):
    with mock.patch('pyar.scripts.scan_bond.run_scan_bond') as run:
        with pytest.raises(SystemExit):
            modern_cli.main([*cli_case, *options])
    run.assert_not_called()


@pytest.fixture
def intramolecular_run(tmp_path, monkeypatch):
    workflow = importlib.import_module('pyar.workflows.scan_bond')
    source = tmp_path / 'molecule.xyz'
    write_xyz(['H', 'H'], [[0., 0., 0.], [.74, 0., 0.]], source)
    calls = []
    reactant = {'inchi': 'reactant', 'smiles': '[H][H]'}
    product = {'inchi': 'product', 'smiles': '[H].[H]'}
    controls = {'changed': True, 'relax_failure': False, 'identity_failure': False}

    def identity(path):
        if Path(path).resolve() == source.resolve():
            return reactant
        if controls['identity_failure']:
            raise ValueError('identity unavailable')
        return product if controls['changed'] else reactant

    monkeypatch.setattr(workflow, 'molecule_identity_from_xyz', identity)
    monkeypatch.setattr(workflow, '_make_orientations', lambda *args: pytest.fail('No fragment placement allowed'))
    monkeypatch.setattr(workflow, 'separated_reactant_identity', lambda *args: pytest.fail('Do not split the input'))

    def scan(molecule, request, directory, parameters):
        calls.append(('scan', molecule.charge, molecule.multiplicity, list(molecule.atoms_list)))
        directory.mkdir(parents=True, exist_ok=True)
        grid = np.linspace(request.start_distance_angstrom, request.end_distance_angstrom, request.n_points)
        frames = [(['H', 'H'], np.array([[0., 0., 0.], [d, 0., 0.]]), '') for d in grid]
        profile = [dict(target_distance_angstrom=float(d), energy_hartree=-1. + .1 * (i == 1))
                   for i, d in enumerate(grid)]
        from pyar.neb import _write_xyz_trajectory
        trajectory, final = directory / 'scan_trajectory.xyz', directory / 'final_scan.xyz'
        _write_xyz_trajectory(trajectory, ['H', 'H'], [f[1] for f in frames], [p['energy_hartree'] for p in profile])
        write_xyz(['H', 'H'], frames[-1][1], final)
        result = OrcaBondScanResult(True, trajectory, final, frames[-1][1], directory / 'scan_request.json',
                                   directory / 'scan_summary.json', 'success', frames, profile)
        controls['scan_result'] = result
        return result

    def relax(molecule, settings):
        calls.append(('relax', molecule.charge, molecule.multiplicity, list(molecule.atoms_list)))
        if controls['relax_failure']:
            return False
        if not controls['changed']:
            molecule.coordinates[1] = [.74, 0., 0.]
        molecule.energy = -1.
        root = Path('job_relaxed')
        root.mkdir()
        write_xyz(molecule.atoms_list, molecule.coordinates, root / 'result_relaxed.xyz', energy=-1.)
        return True

    monkeypatch.setattr('pyar.backends.bond_scan.run_bond_scan', scan)
    monkeypatch.setattr('pyar.backends.bond_scan.load_bond_scan_result', lambda *args: controls['scan_result'])
    monkeypatch.setattr(workflow, 'optimise', relax)
    root = tmp_path / 'run'
    def run(through='scan', **options):
        return run_scan_bond(source, None, (0, 1), 1, {'software': 'xtb', 'opt_cycles': 10}, root,
                             scan_points=4, through=through, **options)
    return root, run, calls, controls


def test_intramolecular_scan_always_relaxes_and_commits_one_input(intramolecular_run):
    root, run, calls, _ = intramolecular_run
    row = run()['results'][0]
    assert [item[0] for item in calls] == ['scan', 'relax']
    assert row['discovery_status'] == 'new_product'
    assert row['target_change_achieved']
    assert row['scan_direction'] == 'break'
    assert not row['target_bond_present_after_relaxation']
    assert len(read_json(root / 'inputs/metadata.json')['inputs']) == 1
    assert read_json(root / 'products/index.json')['products'][0]['target_status'] == 'target'
    commit = sha256_file(root / 'orientation_000/scan_commit.json')
    assert run()['results'][0]['product_identity_changed']
    assert len(calls) == 2
    assert sha256_file(root / 'orientation_000/scan_commit.json') == commit


@pytest.mark.parametrize('changed', [False, True])
def test_followup_uses_changed_relaxed_product_and_full_scan_trace(intramolecular_run, monkeypatch, changed):
    root, run, calls, controls = intramolecular_run
    controls['changed'] = changed
    followed = []
    def path(**kwargs):
        followed.append(kwargs)
        assert Path(kwargs['end']).name == 'result_relaxed.xyz'
        assert Path(kwargs['path_seed']).name == 'scan_trajectory.xyz'
        assert kwargs.get('ts_guess') is None
        return dict(status='complete', completed_stages=[], reactant_product_connection_confirmed=True)
    monkeypatch.setattr('pyar.workflows.scan_path.run_path_stages', path)
    result = run('all')
    assert len(followed) == int(changed)
    assert result['results'][0]['product_identity_changed'] is changed
    assert result['results'][0]['reactant_product_connection_confirmed'] is changed
    assert len(calls) == 2
    if changed:
        index = read_json(root / 'products/index.json')
        assert index['products'][0]['characterization_status'] == 'reaction_connected'


@pytest.mark.parametrize('problem', ['relax_failure', 'identity_failure'])
def test_failed_relaxation_or_unknown_identity_never_launches_path(intramolecular_run, monkeypatch, problem):
    _, run, _, controls = intramolecular_run
    controls[problem] = True
    monkeypatch.setattr('pyar.workflows.scan_path.run_path_stages', lambda **kwargs: pytest.fail('No product certificate'))
    assert run('all')['status'] == 'failed'


def test_intramolecular_formation_uses_the_actual_unplaced_input(intramolecular_run):
    root, run, _, _ = intramolecular_run
    source = root.parent / 'molecule.xyz'
    write_xyz(['H', 'H'], [[0., 0., 0.], [2., 0., 0.]], source)
    result = run()
    row = result['results'][0]
    assert row['scan_direction'] == 'form'
    assert row['start_distance_angstrom'] == 2.
    assert row['target_bond_present_after_relaxation']
    assert row['target_change_achieved']
    np.testing.assert_allclose(Molecule.from_xyz(root / 'orientation_000/start.xyz').coordinates,
                               Molecule.from_xyz(source).coordinates)


def test_single_molecule_resume_rejects_new_geometry_without_new_calculation(intramolecular_run):
    root, run, calls, _ = intramolecular_run
    run()
    write_xyz(['H', 'H'], [[0., 0., 0.], [.8, 0., 0.]], root.parent / 'molecule.xyz')
    with pytest.raises(FileExistsError, match='different request'):
        run()
    assert len(calls) == 2


def test_unavailable_input_identity_stops_before_scan_allocation(intramolecular_run, monkeypatch):
    root, run, calls, _ = intramolecular_run
    def missing(path):
        raise FileNotFoundError('identity executable missing')
    monkeypatch.setattr('pyar.workflows.scan_bond.molecule_identity_from_xyz', missing)
    with pytest.raises(FileNotFoundError, match='identity executable missing'):
        run()
    assert calls == []
    assert not root.exists()
