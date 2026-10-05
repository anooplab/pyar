import json
from unittest.mock import patch

import numpy as np
import pytest

from pyar.modern_cli import main
from pyar.selection.reports import HARTREE_TO_KCAL_MOL
from pyar.structure_inspection import compare_structures
from pyar.structure_comparison.models import ComparisonResult


def xyz(tmp_path, name, atoms=('C', 'O', 'H'), coords=None, comment='-10.0'):
    coords = coords if coords is not None else [(0, 0, 0), (1.3, 0, 0), (0, 1, 0)]
    path = tmp_path / name
    path.write_text(f'{len(atoms)}\n{comment}\n' + ''.join(
        f'{atom} {x} {y} {z}\n' for atom, (x, y, z) in zip(atoms, coords)))
    return str(path)


def test_energies_rank_and_legacy(tmp_path, capsys):
    from pyar.scripts.energy_table import main as legacy
    files = [xyz(tmp_path, f'{i}.xyz', comment=energy)
             for i, energy in enumerate(['energy = -9.8e0', '-1.0E1', '-9.9'])]
    main(['energies', *files, '--json'])
    result = json.loads(capsys.readouterr().out)
    assert result['minimum'] == files[1]
    assert [r['file'] for r in result['structures']] == [files[1], files[2], files[0]]
    assert result['structures'][0]['relative_energy_kcal_mol'] == 0
    assert result['structures'][1]['relative_energy_kcal_mol'] == pytest.approx(.1 * HARTREE_TO_KCAL_MOL)
    main(['energies', *files])
    modern = capsys.readouterr().out
    legacy(files)
    assert capsys.readouterr().out == modern


@pytest.mark.parametrize('bad', ['no energy', ''])
def test_missing_energy_no_partial_table(tmp_path, capsys, bad):
    a = xyz(tmp_path, 'a.xyz')
    b = xyz(tmp_path, 'bad.xyz', comment=bad)
    with pytest.raises(SystemExit) as exc:
        main(['energies', a, b])
    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert not captured.out
    assert b in captured.err


@pytest.mark.parametrize('command', ['energies', 'compare'])
def test_help(command, capsys):
    with pytest.raises(SystemExit) as exc:
        main([command, '--help'])
    assert exc.value.code == 0
    assert f'pyar {command}' in capsys.readouterr().out


def test_top_help(capsys):
    with pytest.raises(SystemExit):
        main(['--help'])
    help_text = capsys.readouterr().out
    for command in ('energies', 'compare', 'clustering', 'optimize', 'scan-bond'):
        assert command in help_text


def test_identical_and_energy(tmp_path):
    a = xyz(tmp_path, 'a.xyz')
    r = compare_structures(a, a)
    assert r['connectivity_match'] and r['same_composition']
    assert r['rmsd_angstrom'] == pytest.approx(0, abs=1e-12)
    assert r['delta_energy_kcal_mol'] == 0
    assert r['added_edges'] == r['removed_edges'] == []
    b = xyz(tmp_path, 'b.xyz', comment='-9.99')
    assert compare_structures(a, b)['delta_energy_kcal_mol'] == pytest.approx(.01 * HARTREE_TO_KCAL_MOL)


def test_rotation_translation_and_permutation(tmp_path):
    coords = np.array([(0, 0, 0), (1.3, 0, 0), (0, 1, 0)])
    a = xyz(tmp_path, 'a.xyz', coords=coords)
    rotated = coords @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + 7
    b = xyz(tmp_path, 'b.xyz', atoms=('O', 'C', 'H'), coords=rotated[[1, 0, 2]])
    r = compare_structures(a, b, atom_mode='all')
    assert r['rmsd_angstrom'] == pytest.approx(0, abs=1e-12)
    assert not r['same_atom_order']
    assert 'added_edges' not in r


def test_equivalent_atom_permutation(tmp_path):
    coords = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    a = xyz(tmp_path, 'a.xyz', atoms=('O', 'H', 'H'), coords=coords)
    b = xyz(tmp_path, 'b.xyz', atoms=('O', 'H', 'H'), coords=[coords[0], coords[2], coords[1]])
    assert compare_structures(a, b, atom_mode='all')['rmsd_angstrom'] == pytest.approx(0, abs=1e-12)


def test_geometry_change_without_connectivity_change(tmp_path):
    a = xyz(tmp_path, 'a.xyz')
    b = xyz(tmp_path, 'b.xyz', coords=[(0, 0, 0), (1.4, 0, 0), (0, 1, 0)])
    r = compare_structures(a, b)
    assert r['connectivity_match']
    assert r['rmsd_angstrom'] > 0
    assert r['geometry_equivalent_under_threshold'] is None


def test_edge_changes_and_hydrogen_fallback(tmp_path):
    a = xyz(tmp_path, 'a.xyz', atoms=('H', 'H'), coords=[(0, 0, 0), (2, 0, 0)])
    b = xyz(tmp_path, 'b.xyz', atoms=('H', 'H'), coords=[(0, 0, 0), (.6, 0, 0)])
    r = compare_structures(a, b)
    assert r['added_edges'] == [[0, 1]]
    assert (r['input_component_count'], r['output_component_count']) == (2, 1)
    assert r['connectivity_match'] is False and r['rmsd_angstrom'] is None
    assert r['rmsd_atom_mode'] == 'all'
    assert compare_structures(b, a)['removed_edges'] == [[0, 1]]


def test_composition_mismatch_and_missing_energy(tmp_path):
    a = xyz(tmp_path, 'a.xyz', comment='no energy')
    b = xyz(tmp_path, 'b.xyz', atoms=('H',), coords=[(0, 0, 0)])
    r = compare_structures(a, b)
    assert r['first_energy_hartree'] is None
    assert not r['same_composition']
    assert r['connectivity_match'] is None and r['rmsd_angstrom'] is None
    assert any('different compositions' in item for item in r['limitations'])
    r = compare_structures(a, a)
    assert r['connectivity_match'] and r['delta_energy_hartree'] is None


def test_incomplete_mapping(tmp_path, capsys):
    a = xyz(tmp_path, 'a.xyz')
    incomplete = ComparisonResult(True, None, None, 'graph', metadata={
        'connectivity_match': True, 'comparison_complete': False, 'isomorphisms_evaluated': 10})
    with patch('pyar.structure_comparison.pairwise.GraphRMSDComparator.compare', return_value=incomplete):
        main(['compare', a, a])
    text = capsys.readouterr().out
    assert 'incomplete' in text and 'limit reached' in text
    assert 'Connectivity match: yes' in text


def test_json_and_no_side_effects(tmp_path, capsys, monkeypatch):
    a = xyz(tmp_path, 'a.xyz')
    b = xyz(tmp_path, 'b.xyz', comment='geometry only')
    monkeypatch.chdir(tmp_path)
    before = set(tmp_path.iterdir())
    main(['compare', a, b, '--json'])
    r = json.loads(capsys.readouterr().out)
    assert r['second_energy_hartree'] is None
    assert set(tmp_path.iterdir()) == before


@pytest.mark.parametrize('options', [['--bond-scale', 'nan'], ['--maximum-mappings', '0'],
                                     ['--rmsd-threshold', '-1']])
def test_invalid_options(tmp_path, options):
    a = xyz(tmp_path, 'a.xyz')
    with pytest.raises(SystemExit) as exc:
        main(['compare', a, a, *options])
    assert exc.value.code == 2


@pytest.mark.parametrize('command', ['compare', 'energies'])
@pytest.mark.parametrize('bad', ['missing', 'malformed', 'nonfinite'])
def test_invalid_geometry(tmp_path, capsys, command, bad):
    a = xyz(tmp_path, 'a.xyz')
    b = str(tmp_path / 'bad.xyz')
    if bad == 'malformed':
        (tmp_path / 'bad.xyz').write_text('2\n-10\nH 0 0 0\n')
    if bad == 'nonfinite':
        (tmp_path / 'bad.xyz').write_text('1\n-10\nH nan 0 0\n')
    with pytest.raises(SystemExit) as exc:
        main([command, a, b])
    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert not captured.out and 'bad.xyz' in captured.err


def test_single_energy_and_required_input(tmp_path, capsys):
    a = xyz(tmp_path, 'a.xyz')
    main(['energies', a, '--json'])
    result = json.loads(capsys.readouterr().out)
    assert result['minimum'] == a
    assert result['structures'][0]['relative_energy_kcal_mol'] == 0
    with pytest.raises(SystemExit) as exc:
        main(['energies'])
    assert exc.value.code == 2
