"""Optional RDKit perception and independent coordinate analysis regressions."""

import json
from unittest.mock import patch

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.modern_cli import main
from pyar.structure_comparison.chemical_identity import (
    INSTALL_HINT, perceive_chemical_identity, resolve_charge,
)
from pyar.structure_comparison.pairwise import compare_structures


def water():
    return Molecule(['O', 'H', 'H'], [[0, 0, 0], [.9572, 0, 0], [-.24, .927, 0]],
                    charge=None, multiplicity=None)


def write_xyz(tmp_path, name, molecule):
    path = tmp_path / name
    path.write_text(f'{len(molecule.atoms_list)}\n-10\n' + ''.join(
        f'{s} {x} {y} {z}\n' for s, (x, y, z) in zip(molecule.atoms_list, molecule.coordinates)))
    return str(path)


@pytest.mark.parametrize('smiles', ['O', 'C', 'CCO', 'C=C', 'C#C', 'C=O'])
def test_rdkit_xyz_perception_and_perturbation(smiles):
    Chem = pytest.importorskip('rdkit.Chem')
    AllChem = pytest.importorskip('rdkit.Chem.AllChem')
    pytest.importorskip('rdkit.Chem.rdDetermineBonds')
    reference = Chem.MolFromSmiles(smiles)
    fixture = Chem.AddHs(reference)
    assert AllChem.EmbedMolecule(fixture, randomSeed=123) == 0
    symbols = [atom.GetSymbol() for atom in fixture.GetAtoms()]
    coordinates = np.array(fixture.GetConformer().GetPositions())
    molecule = Molecule(symbols, coordinates, charge=None, multiplicity=None)
    result = perceive_chemical_identity(molecule)
    assert result['perception_success'], result
    assert result['canonical_smiles'] == Chem.MolToSmiles(reference, canonical=True, isomericSmiles=True)
    assert result['charge_used'] == 0 and result['charge_source'] == 'assumed'
    perturbed = Molecule(symbols, coordinates + np.random.default_rng(5).normal(0, .002, coordinates.shape))
    assert perceive_chemical_identity(perturbed)['canonical_smiles'] == result['canonical_smiles']


def test_charge_sensitive_hydroxide():
    pytest.importorskip('rdkit.Chem.rdDetermineBonds')
    molecule = Molecule(['O', 'H'], [[0, 0, 0], [.97, 0, 0]])
    result = perceive_chemical_identity(molecule, charge=-1)
    assert result['perception_success'], result
    assert result['canonical_smiles'] == '[OH-]'
    assert result['charge_source'] == 'explicit'
    # Neutral cannot be made valid by silently searching for another charge.
    neutral = perceive_chemical_identity(molecule)
    assert not neutral['perception_success']
    assert neutral['charge_used'] == 0 and neutral['canonical_smiles'] is None


def test_charge_sources():
    assert resolve_charge() == (0, 'assumed')
    assert resolve_charge(known_charge=-1) == (-1, 'known')
    assert resolve_charge(0, known_charge=-1) == (0, 'explicit')
    with pytest.raises(ValueError):
        resolve_charge(.5)


@pytest.mark.parametrize('charges,expected', [(['0'], [0, 0]), (['0', '-1'], [0, -1]),
                                            (['-1'], [-1, -1])])
def test_cli_charge_resolution(tmp_path, capsys, charges, expected):
    first = write_xyz(tmp_path, 'a.xyz', water())
    second = write_xyz(tmp_path, 'b.xyz', water())
    with patch('pyar.structure_comparison.chemical_identity._rdkit_api', side_effect=ImportError('no rdkit')):
        main(['compare', first, second, '--charge', *charges, '--json'])
    chemical = json.loads(capsys.readouterr().out)['chemical_identity']
    assert [chemical[side]['charge_used'] for side in ('first', 'second')] == expected
    assert all(chemical[side]['charge_source'] == 'explicit' for side in ('first', 'second'))


def test_invalid_charge_count(tmp_path):
    first = write_xyz(tmp_path, 'a.xyz', water())
    with pytest.raises(SystemExit) as exc:
        main(['compare', first, first, '--charge', '0', '0', '0'])
    assert exc.value.code == 2


def test_known_charge_in_reusable_api():
    with patch('pyar.structure_comparison.chemical_identity._rdkit_api', side_effect=ImportError('no rdkit')):
        result = compare_structures(water(), water(), known_charges=(0, -1))
    assert result['chemical_identity']['first']['charge_source'] == 'known'
    assert result['chemical_identity']['second']['charge_used'] == -1


def test_optional_rdkit_unavailable(tmp_path, capsys):
    first = write_xyz(tmp_path, 'a.xyz', water())
    with patch('pyar.structure_comparison.chemical_identity._rdkit_api', side_effect=ModuleNotFoundError('RDKit is not installed')):
        main(['compare', first, first, '--json'])
        result = json.loads(capsys.readouterr().out)
        main(['compare', first, first])
    assert result['connectivity_match'] and result['rmsd_angstrom'] == pytest.approx(0, abs=1e-10)
    assert result['delta_energy_kcal_mol'] == 0
    identity = result['chemical_identity']['first']
    assert identity['status'] == 'unavailable' and identity['canonical_smiles'] is None
    assert identity['installation_hint'] == INSTALL_HINT
    assert result['chemical_identity']['canonical_smiles_match'] is None
    assert INSTALL_HINT in capsys.readouterr().out


def test_perception_failure_independent_of_topology(tmp_path, capsys):
    pytest.importorskip('rdkit.Chem.rdDetermineBonds')
    first = write_xyz(tmp_path, 'a.xyz', water())
    with patch('rdkit.Chem.rdDetermineBonds.DetermineBonds', side_effect=ValueError('unsupported valence')) as determine:
        main(['compare', first, first, '--json'])
    result = json.loads(capsys.readouterr().out)
    assert determine.call_count == 2  # exactly one attempt per structure; no charge search
    identity = result['chemical_identity']['first']
    assert identity['status'] == 'failed' and 'unsupported valence' in identity['reason']
    assert identity['canonical_smiles'] is None
    assert result['connectivity_match'] and result['comparison_complete']


def test_smiles_match_independent_of_raw_atom_order(tmp_path, capsys):
    pytest.importorskip('rdkit.Chem.rdDetermineBonds')
    molecule = water()
    first = write_xyz(tmp_path, 'a.xyz', molecule)
    reordered = Molecule(['H', 'O', 'H'], molecule.coordinates[[1, 0, 2]], charge=None, multiplicity=None)
    second = write_xyz(tmp_path, 'b.xyz', reordered)
    main(['compare', first, second, '--json'])
    result = json.loads(capsys.readouterr().out)
    assert result['chemical_identity']['canonical_smiles_match'] is True
    assert not result['same_atom_order'] and 'added_edges' not in result
    assert result['rmsd_angstrom'] == pytest.approx(0, abs=1e-10)
