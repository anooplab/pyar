import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from pyar.modern_cli import main
from pyar.geometry_utilities import orient_structures
from pyar.core.molecule import parse_xyz
from pyar.selection.deduplication import deduplicate_structures
from pyar.selection.energy_window import select_structures
from pyar.selection.reports import HARTREE_TO_KCAL_MOL
from pyar.structure_comparison.models import ComparisonResult
from pyar.structure_inspection import format_formula
from pyar.utility_io import load_structures


def xyz(tmp_path, name, atoms=('O', 'H', 'H'), coords=None, energy=None):
    coords = coords if coords is not None else [[0, 0, 0], [.95, 0, 0], [-.24, .92, 0]]
    path = tmp_path / name
    path.write_text(f'{len(atoms)}\n' + ('geometry only' if energy is None else str(energy)) + '\n' + ''.join(
        f'{s} {x:.17g} {y:.17g} {z:.17g}\n' for s, (x, y, z) in zip(atoms, coords)))
    return str(path)


@pytest.mark.parametrize('command', ['orient', 'deduplicate', 'select', 'identify', 'split', 'trace'])
def test_command_help(command, capsys):
    with pytest.raises(SystemExit) as exc:
        main([command, '--help'])
    assert exc.value.code == 0
    assert f'pyar {command}' in capsys.readouterr().out


def test_top_level_help(capsys):
    with pytest.raises(SystemExit):
        main(['--help'])
    output = capsys.readouterr().out
    assert all(command in output for command in ['orient', 'deduplicate', 'select', 'identify', 'split', 'trace'])


def test_orient_defaults_deterministic_and_inputs_untouched(tmp_path, monkeypatch):
    a = xyz(tmp_path, 'a.xyz')
    b = xyz(tmp_path, 'b.xyz')
    inputs = [Path(a).read_bytes(), Path(b).read_bytes()]
    monkeypatch.chdir(tmp_path)
    # None of the utility paths has permission to use computational execution.
    with patch('pyar.optimiser.bulk_optimize', side_effect=AssertionError('calculation')), patch(
            'pyar.sampling.trial_generator.merge_two_molecules', wraps=__import__(
                'pyar.sampling.trial_generator', fromlist=['merge_two_molecules']).merge_two_molecules) as merge:
        main(['orient', a, b])
        assert merge.call_count == 8
        assert all(call.kwargs['distance_scaling'] == 1.5 for call in merge.call_args_list)
    first = tmp_path / 'orientations'
    assert len(list(first.glob('*.xyz'))) == 8
    assert (first / 'trial_vectors.dat').exists()
    main(['orient', a, b, '--output', 'again'])
    assert {p.name: p.read_bytes() for p in first.iterdir()} == {
        p.name: p.read_bytes() for p in (tmp_path / 'again').iterdir()}
    assert inputs == [Path(a).read_bytes(), Path(b).read_bytes()]
    with pytest.raises(SystemExit):
        main(['orient', a, b])


def test_orient_atomic_and_offset(tmp_path):
    a = xyz(tmp_path, 'a.xyz')
    b = xyz(tmp_path, 'atom.xyz', ('H',), [[0, 0, 0]])
    first = orient_structures(a, b, orientations=20, output=tmp_path/'one')
    assert len(first['files']) == 20 and first['distance_scale'] == 1.2
    vectors = np.loadtxt(first['trial_vectors'])
    assert np.all(vectors[:, 3:] == 0)
    second = orient_structures(a, b, orientations=20, sequence_offset=2, output=tmp_path/'two')
    assert not np.array_equal(vectors, np.loadtxt(second['trial_vectors']))
    assert np.all(np.loadtxt(second['trial_vectors'])[:, 3:] == 0)


@pytest.mark.parametrize('option,value', [('--orientations', '0'), ('--distance-scale', 'nan'),
                                        ('--sequence-offset', '-1')])
def test_orient_invalid_no_output(tmp_path, option, value):
    a = xyz(tmp_path, 'a.xyz')
    directory = tmp_path / 'out'
    with pytest.raises(SystemExit):
        main(['orient', a, a, option, value, '--output', str(directory)])
    assert not directory.exists()


def test_dedup_no_energy_and_rigid_transform(tmp_path, monkeypatch, capsys):
    a = xyz(tmp_path, 'a.xyz')
    coords = np.array([[0, 0, 0], [.95, 0, 0], [-.24, .92, 0]])
    b = xyz(tmp_path, 'b.xyz', coords=coords @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + 3)
    monkeypatch.chdir(tmp_path)
    before = set(tmp_path.iterdir())
    main(['deduplicate', a, b])
    assert 'Unique structures: 1' in capsys.readouterr().out
    assert set(tmp_path.iterdir()) == before
    result = deduplicate_structures(load_structures([a, b]))
    assert [m.name for m in result['kept']] == [a]
    main(['deduplicate', a, b, '--output', 'unique'])
    assert [p.name for p in (tmp_path/'unique').iterdir()] == ['a.xyz']


def test_dedup_energy_representative_and_missing_partial(tmp_path):
    a = xyz(tmp_path, 'high.xyz', energy=-9)
    b = xyz(tmp_path, 'low.xyz', energy=-10)
    assert [m.name for m in deduplicate_structures(load_structures([a, b]))['kept']] == [b]
    c = xyz(tmp_path, 'missing.xyz')
    assert [m.name for m in deduplicate_structures(load_structures([a, c, b]))['kept']] == [a]


def test_dedup_topology_and_distinct_geometry(tmp_path):
    a = xyz(tmp_path, 'a.xyz', ('C', 'O'), [[0, 0, 0], [1.3, 0, 0]])
    b = xyz(tmp_path, 'b.xyz', ('C', 'O'), [[0, 0, 0], [4, 0, 0]])
    c = xyz(tmp_path, 'c.xyz', ('C', 'O'), [[0, 0, 0], [1.62, 0, 0]])
    result = deduplicate_structures(load_structures([a, b, c]))
    assert len(result['kept']) == 3


@pytest.mark.parametrize('result', [RuntimeError('unavailable'), ComparisonResult(
    True, None, None, 'graph', metadata={'comparison_complete': False})])
def test_dedup_failure_retains(tmp_path, result):
    files = [xyz(tmp_path, name) for name in ['a.xyz', 'b.xyz']]
    with patch('pyar.structure_comparison.GraphFirstDeduplicationComparator.compare') as compare:
        if isinstance(result, Exception):
            compare.side_effect = result
        else:
            compare.return_value = result
        report = deduplicate_structures(load_structures(files))
    assert len(report['kept']) == 2 and report['diagnostics']


def test_select_window_top_unique(tmp_path):
    files = [xyz(tmp_path, f'{index}.xyz', energy=relative/HARTREE_TO_KCAL_MOL)
             for index, relative in enumerate([0., 2., 5., 5.1])]
    molecules = load_structures(files[::-1], require_energy=True)
    assert [m.name for m in select_structures(molecules, within=5)['kept']] == files[:3]
    assert [m.name for m in select_structures(molecules, top=2)['kept']] == files[:2]
    assert [m.name for m in select_structures(molecules, within=5, top=2)['kept']] == files[:2]
    assert [m.name for m in select_structures(molecules, within=5, unique=True)['kept']] == files[:1]


def test_select_missing_energy_and_output(tmp_path, capsys, monkeypatch):
    a = xyz(tmp_path, 'a.xyz', energy=-10)
    b = xyz(tmp_path, 'missing.xyz')
    with pytest.raises(SystemExit):
        main(['select', a, b, '--within', '5'])
    captured = capsys.readouterr()
    assert b in captured.err and not captured.out
    monkeypatch.chdir(tmp_path)
    before = set(tmp_path.iterdir())
    main(['select', a, '--top', '1'])
    assert set(tmp_path.iterdir()) == before
    main(['select', a, '--within', '5', '--output', 'selected'])
    assert [p.name for p in (tmp_path/'selected').iterdir()] == ['a.xyz']


@pytest.mark.parametrize('options', [[], ['--unique'], ['--within', '-1'], ['--within', 'nan'], ['--top', '0']])
def test_select_invalid(tmp_path, options):
    a = xyz(tmp_path, 'a.xyz', energy=-10)
    with pytest.raises(SystemExit):
        main(['select', a, *options])


def test_identify_basic_optional_and_multiple(tmp_path, capsys):
    a = xyz(tmp_path, 'a.xyz', energy=-76)
    b = xyz(tmp_path, 'b.xyz')
    with patch('pyar.structure_comparison.chemical_identity._rdkit_api', side_effect=ImportError('no rdkit')):
        main(['identify', a, b, '--charge', '0', '-1', '--json'])
    results = json.loads(capsys.readouterr().out)['structures']
    assert len(results) == 2
    assert results[0]['formula'] == 'H2O' and results[0]['atom_count'] == 3
    assert results[0]['component_count'] == 1 and results[0]['energy_hartree'] == -76
    assert results[1]['energy_hartree'] is None
    assert results[1]['chemical_identity']['charge_used'] == -1
    assert results[0]['chemical_identity']['charge_source'] == 'explicit'
    assert results[0]['chemical_identity']['installation_hint'] == 'pip install "pyar-chem[identity]"'


def test_identify_rdkit_and_failure(tmp_path, capsys):
    pytest.importorskip('rdkit.Chem.rdDetermineBonds')
    a = xyz(tmp_path, 'a.xyz')
    main(['identify', a, '--json'])
    identity = json.loads(capsys.readouterr().out)['structures'][0]['chemical_identity']
    assert identity['canonical_smiles'] == 'O' and identity['charge_source'] == 'assumed'
    with patch('rdkit.Chem.rdDetermineBonds.DetermineBonds', side_effect=ValueError('bad valence')):
        main(['identify', a, '--json'])
    result = json.loads(capsys.readouterr().out)['structures'][0]
    assert result['component_count'] == 1
    assert result['chemical_identity']['status'] == 'failed'


@pytest.mark.parametrize('composition,formula', [({'C': 2, 'H': 6, 'O': 1}, 'C2H6O'),
                                               ({'C': 1, 'H': 4, 'N': 2, 'O': 1}, 'CH4N2O'),
                                               ({'H': 2, 'O': 1}, 'H2O'), ({'Na': 1, 'Cl': 1}, 'ClNa')])
def test_hill_formula(composition, formula):
    assert format_formula(composition) == formula


def test_split_order_coordinates_and_provenance(tmp_path, monkeypatch):
    atoms = ['O', 'C', 'H', 'H', 'H']
    coords = np.array([[0, 0, 0], [10, 0, 0], [.9, 0, 0], [10.9, 0, 0], [-.2, .9, 0]])
    a = xyz(tmp_path, 'complex.xyz', atoms, coords)
    monkeypatch.chdir(tmp_path)
    main(['split', a])
    files = sorted((tmp_path/'complex_fragments').glob('*.xyz'))
    assert len(files) == 2
    for path, indices in zip(files, [[0, 2, 4], [1, 3]]):
        symbols, coordinates, _, title, _ = parse_xyz(str(path))
        assert symbols == [atoms[index] for index in indices]
        assert np.array_equal(coordinates, coords[indices])
        assert 'original atom indices:' in title
        assert 'energy unavailable' in title
    with pytest.raises(SystemExit):
        main(['split', a])


def test_split_single_and_bond_scale(tmp_path, capsys, monkeypatch):
    a = xyz(tmp_path, 'a.xyz', ('H', 'H'), [[0, 0, 0], [.8, 0, 0]])
    monkeypatch.chdir(tmp_path)
    main(['split', a, '--bond-scale', '2'])
    assert 'nothing to split' in capsys.readouterr().out
    assert not (tmp_path/'a_fragments').exists()
    main(['split', a])
    assert len(list((tmp_path/'a_fragments').glob('*.xyz'))) == 2


@pytest.mark.parametrize('command', ['deduplicate', 'select'])
def test_copy_safety(tmp_path, command):
    a = xyz(tmp_path, 'a.xyz', energy=-10)
    out = tmp_path/'out'
    out.mkdir()
    existing = out/'unrelated.txt'
    existing.write_text('keep')
    options = ['--top', '1'] if command == 'select' else []
    with pytest.raises(SystemExit):
        main([command, a, *options, '--output', str(out)])
    assert existing.read_text() == 'keep' and list(out.iterdir()) == [existing]


@pytest.mark.parametrize('options', [[], ['--plot'], ['--plot-only', '--plot-directory', 'plots'],
                                     ['--max-force', '.2', '--exclude-energy-outliers', '3.5']])
def test_trace_shared_execution(options, capsys):
    from pyar.scripts.reaction_trace import main as legacy
    with patch('pyar.scripts.reaction_trace.analyse_reaction_trace', return_value={'count': 2}) as analyze, patch(
            'pyar.scripts.reaction_trace.plot_reaction_trace', return_value={'files': ['plot.png']}) as plot:
        main(['trace', 'RUN', *options])
        modern_output = capsys.readouterr().out
        modern_calls = (analyze.call_args_list[:], plot.call_args_list[:])
        analyze.reset_mock()
        plot.reset_mock()
        legacy(['RUN', *options])
        assert capsys.readouterr().out == modern_output
        assert (analyze.call_args_list, plot.call_args_list) == modern_calls


def test_generated_components_have_no_fabricated_energy(tmp_path):
    from pyar.geometry_utilities import split_structure
    a = xyz(tmp_path, 'a.xyz', ('H', 'H'), [[0, 0, 0], [3, 0, 0]], energy=-10)
    result = split_structure(a, output=tmp_path/'fragments')
    assert all(m.energy is None for m in load_structures(result['files']))


def test_unknown_merge_state_and_historical_state():
    from pyar.core.molecule import Molecule
    first = Molecule(['H'], [[0, 0, 0]], charge=None, multiplicity=None)
    second = Molecule(['H'], [[2, 0, 0]], charge=None, multiplicity=None)
    merged = first.merged_with(second)
    assert merged.charge is None and merged.multiplicity is None
    first.charge, first.multiplicity = -1, 2
    second.charge, second.multiplicity = 0, 2
    merged = first.merged_with(second)
    assert merged.charge == -1 and merged.multiplicity == 1


def test_collision_before_copying(tmp_path):
    from pyar.utility_io import copy_structures
    (tmp_path/'one').mkdir()
    (tmp_path/'two').mkdir()
    a = xyz(tmp_path/'one', 'a.xyz')
    b = xyz(tmp_path/'two', 'a.xyz')
    with pytest.raises(ValueError, match='colliding'):
        copy_structures(load_structures([a, b]), tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_inclusive_window_with_large_negative_energies(tmp_path):
    minimum = -1000.
    files = [xyz(tmp_path, f'{index}.xyz', energy=minimum + relative/HARTREE_TO_KCAL_MOL)
             for index, relative in enumerate([0., 5., 5.0001])]
    selected = select_structures(load_structures(files), within=5)['kept']
    assert [m.name for m in selected] == files[:2]


def test_trace_real_fixture(tmp_path, capsys):
    from pyar.reaction_trace import ReactionTraceRecorder
    recorder = ReactionTraceRecorder(str(tmp_path))
    for step in range(3):
        recorder.record(symbols=['H', 'H'], coordinates_angstrom=[[0, 0, 0], [1.5-step*.2, 0, 0]],
                        backend_energy_hartree=[-1., -.9, -1.01][step],
                        afir_energy_hartree=0., total_energy_hartree=[-1., -.9, -1.01][step],
                        backend_forces_hartree_per_bohr=np.zeros((2, 3)),
                        afir_forces_hartree_per_bohr=np.zeros((2, 3)),
                        total_forces_hartree_per_bohr=np.zeros((2, 3)),
                        backend_force_norm=0., afir_force_norm=0., total_force_norm=0., max_force=0.,
                        fragment_indices=[[0], [1]])
    main(['trace', str(tmp_path)])
    summary = json.loads(capsys.readouterr().out)
    assert 'analysis' in summary
    assert (tmp_path/'path_summary.csv').is_file()
    from pyar.scripts.reaction_trace import main as legacy
    first_artifact = (tmp_path/'path_summary.csv').read_bytes()
    legacy([str(tmp_path/'reaction_trace')])
    assert (tmp_path/'path_summary.csv').read_bytes() == first_artifact


@pytest.mark.parametrize('alias', ['--orientations', '-N'])
def test_orient_count_alias_and_molecular_rotations(tmp_path, alias):
    a = xyz(tmp_path, 'a.xyz')
    out = tmp_path/'out'
    main(['orient', a, a, alias, '20', '--output', str(out)])
    assert len(list(out.glob('*.xyz'))) == 20
    assert np.any(np.loadtxt(out/'trial_vectors.dat')[:, 3:] != 0)


def test_identity_extra_optional():
    import tomllib
    project = tomllib.loads((Path(__file__).parents[1]/'pyproject.toml').read_text())['project']
    assert project['optional-dependencies']['identity'] == ['rdkit']
    assert 'rdkit' in project['optional-dependencies']['conformer']
    assert not any(dependency.startswith('rdkit') for dependency in project['dependencies'])
