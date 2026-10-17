import json
from pathlib import Path

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.workflows import association_reference as association


def test_separated_reactant_reference_optimizes_and_validates_fragments_once(tmp_path, monkeypatch):
    reactants = [Molecule(['H', 'H'], [[0, 0, 0], [0.74, 0, 0]]),
                 Molecule(['He'], [[0, 0, 0]])]
    reactants[0].charge, reactants[0].multiplicity = 0, 1
    reactants[1].charge, reactants[1].multiplicity = 0, 1
    calculators, optimizations, frequencies = [], [], []
    monkeypatch.setattr(association, 'write_xyz', association.write_xyz)
    monkeypatch.setattr('pyar.backends.geometric.PyarGeometricCalculator',
                        lambda qc_params: calculators.append(dict(qc_params)) or object())
    monkeypatch.setattr('pyar.reaction_identity.molecule_identity_from_xyz', lambda path: {'inchi': Path(path).name})
    monkeypatch.setattr('pyar.reaction_identity.same_molecular_identity', lambda first, second: True)
    monkeypatch.setattr('pyar.reaction_identity.separated_reactant_identity', lambda first, second: {'inchi': 'AB'})

    def optimize(symbols, coordinates, calculator, output, label, max_cycles, **kwargs):
        optimizations.append((label, max_cycles, kwargs['fmax']))
        return [np.asarray(coordinates)], [-1.0], True

    def frequency(symbols, coordinates, calculator, output, label, threshold):
        frequencies.append(label)
        return {'minimum_confirmed': True, 'energy_hartree': -0.5 if label == 'fragment_1' else -0.25,
                'evaluated_coordinates_angstrom': np.asarray(reactants[len(frequencies) - 1].coordinates).tolist(),
                'frequencies_cm-1': [100.0], 'imaginary_frequencies_cm-1': []}

    monkeypatch.setattr('pyar.neb._optimize_geometry', optimize)
    monkeypatch.setattr('pyar.neb._frequency', frequency)
    params = {'software': 'xtb', 'xtb_model': 'gfn2', 'nprocs': 1}
    options = {'product_relaxation_max_steps': 75, 'product_relaxation_fmax': 0.03}
    first = association.prepare_separated_reactant_reference(reactants, params, tmp_path, options)
    assert first['status'] == 'complete'
    assert first['reference_energy_hartree'] == -0.75
    assert [row['multiplicity'] for row in calculators] == [1, 1]
    assert [item[1] for item in optimizations] == [75, 75]
    assert frequencies == ['fragment_1', 'fragment_2']
    second = association.prepare_separated_reactant_reference(reactants, params, tmp_path, options)
    assert second == first
    assert len(optimizations) == len(frequencies) == 2
    assert json.loads((tmp_path / 'reference.json').read_text())['status'] == 'complete'


def test_separated_reactant_reference_fails_closed_for_nonminimum_fragment(tmp_path, monkeypatch):
    reactants = [Molecule(['H'], [[0, 0, 0]]), Molecule(['He'], [[0, 0, 0]])]
    for molecule in reactants:
        molecule.charge, molecule.multiplicity = 0, 1
    monkeypatch.setattr('pyar.backends.geometric.PyarGeometricCalculator', lambda **kwargs: object())
    monkeypatch.setattr('pyar.reaction_identity.molecule_identity_from_xyz', lambda path: {'inchi': 'same'})
    monkeypatch.setattr('pyar.reaction_identity.same_molecular_identity', lambda first, second: True)
    monkeypatch.setattr('pyar.reaction_identity.separated_reactant_identity', lambda first, second: {'inchi': 'AB'})
    monkeypatch.setattr('pyar.neb._optimize_geometry',
                        lambda symbols, coordinates, calculator, output, label, max_cycles, **kwargs:
                        ([np.asarray(coordinates)], [-1.0], True))
    monkeypatch.setattr('pyar.neb._frequency',
                        lambda *args: {'minimum_confirmed': False, 'energy_hartree': -1.0})
    result = association.prepare_separated_reactant_reference(
        reactants, {'software': 'xtb'}, tmp_path, {})
    assert result['status'] == 'failed'
    assert result['reference_energy_hartree'] is None
    assert 'not a confirmed minimum' in result['error']


def test_association_endpoint_matching_ignores_relative_fragment_pose(tmp_path, monkeypatch):
    from pyar.backends import write_xyz
    symbols = ['H', 'H', 'He']
    reactants = [Molecule(['H', 'H'], [[0, 0, 0], [1, 0, 0]]),
                 Molecule(['He'], [[0, 0, 0]])]
    reference_dir = tmp_path / 'reference'
    reference_dir.mkdir()
    write_xyz(['H', 'H'], [[0, 0, 0], [1, 0, 0]], reference_dir / 'a.xyz')
    write_xyz(['He'], [[0, 0, 0]], reference_dir / 'b.xyz')
    product_path = tmp_path / 'product.xyz'
    product = np.array([[0., 0., 0.], [1.2, 0., 0.], [1.5, 0., 0.]])
    write_xyz(symbols, product, product_path)
    reference = {'status': 'complete', 'request': {'combined_reactant_identity': {'inchi': 'AB'}},
                 'fragments': [
        {'identity': {'inchi': 'A'}, 'optimized_geometry': 'a.xyz'},
        {'identity': {'inchi': 'B'}, 'optimized_geometry': 'b.xyz'},
    ]}

    def identity(path):
        path = Path(path)
        if path.name == 'product.xyz':
            return {'inchi': 'P'}
        elements = path.read_text().splitlines()
        count = int(elements[0])
        rows = elements[2:2 + count]
        if count == 2:
            return {'inchi': 'A'}
        if count == 1:
            return {'inchi': 'B'}
        coords = np.asarray([[float(x) for x in row.split()[1:4]] for row in rows])
        return {'inchi': 'P' if np.linalg.norm(coords[2] - coords[1]) < 1.0 else 'AB'}

    monkeypatch.setattr('pyar.reaction_identity.molecule_identity_from_xyz', identity)
    monkeypatch.setattr('pyar.reaction_identity.same_molecular_identity',
                        lambda first, second: first.get('inchi') == second.get('inchi'))
    reactant_branch = np.array([[10., 4., 0.], [10., 5., 0.], [100., -40., 2.]])
    result = association.classify_association_endpoints(
        symbols, [product, reactant_branch], product_path, reference, reactants, 0.01,
        reference_directory=reference_dir)
    assert result['status'] == 'matched'
    assert result['reactant_product_assignment'] == [1, 0]
    assert result['branch_results'][1]['classification'] == 'separated_fragments_minima'
    assert result['reactant_product_connection_confirmed']

    def wrong_fragment_identity(path):
        if Path(path).name == 'fragment_2.xyz':
            return {'inchi': 'wrong-fragment'}
        return identity(path)

    monkeypatch.setattr('pyar.reaction_identity.molecule_identity_from_xyz', wrong_fragment_identity)
    wrong = association.classify_association_endpoints(
        symbols, [product, reactant_branch], product_path, reference, reactants, 0.01,
        reference_directory=reference_dir)
    assert wrong['status'] == 'unmatched'
    assert wrong['branch_results'][1]['fragment_identity_matches'] == [True, False]


def test_association_endpoint_accepts_stable_pre_reactive_complex(tmp_path, monkeypatch):
    from pyar.backends import write_xyz
    symbols = ['H', 'He']
    reactants = [Molecule(['H'], [[0, 0, 0]]), Molecule(['He'], [[0, 0, 0]])]
    reference_dir = tmp_path / 'reference'
    reference_dir.mkdir()
    write_xyz(['H'], [[0, 0, 0]], reference_dir / 'a.xyz')
    write_xyz(['He'], [[0, 0, 0]], reference_dir / 'b.xyz')
    product_path = tmp_path / 'product.xyz'
    write_xyz(symbols, [[0, 0, 0], [2.0, 0, 0]], product_path)
    reference = {'status': 'complete', 'request': {'combined_reactant_identity': {'inchi': 'R'}},
                 'fragments': [{'identity': {'inchi': 'A'}, 'optimized_geometry': 'a.xyz'},
                               {'identity': {'inchi': 'B'}, 'optimized_geometry': 'b.xyz'}]}

    def identity(path):
        name = Path(path).name
        if name == 'endpoint.xyz':
            rows = Path(path).read_text().splitlines()[2:4]
            distance = np.linalg.norm(np.array([float(x) for x in rows[0].split()[1:4]]) -
                                      np.array([float(x) for x in rows[1].split()[1:4]]))
            return {'inchi': 'P' if distance > 1.0 else 'R'}
        return {'inchi': {'product.xyz': 'P', 'fragment_1.xyz': 'A',
                          'fragment_2.xyz': 'B'}.get(name, 'R')}

    monkeypatch.setattr('pyar.reaction_identity.molecule_identity_from_xyz', identity)
    monkeypatch.setattr('pyar.reaction_identity.same_molecular_identity',
                        lambda first, second: first.get('inchi') == second.get('inchi'))
    precomplex = np.array([[0., 0, 0], [.6, 0, 0]])
    product_branch = np.array([[0., 0, 0], [2., 0, 0]])
    result = association.classify_association_endpoints(
        symbols, [product_branch, precomplex], product_path, reference, reactants, .01,
        reference_directory=reference_dir)
    assert result['status'] == 'matched'
    assert result['branch_results'][1]['classification'] == 'stable_pre_reactive_minimum'
    assert result['reactant_product_connection_confirmed']


@pytest.mark.parametrize("second_distance", [0.1, 1.2])
def test_two_irc_branches_to_one_product_do_not_confirm_association(tmp_path, monkeypatch, second_distance):
    from pyar.backends import write_xyz
    symbols = ['H', 'He']
    reactants = [Molecule(['H'], [[0, 0, 0]]), Molecule(['He'], [[0, 0, 0]])]
    reference_dir = tmp_path / 'reference'
    reference_dir.mkdir()
    write_xyz(['H'], [[0, 0, 0]], reference_dir / 'a.xyz')
    write_xyz(['He'], [[0, 0, 0]], reference_dir / 'b.xyz')
    product_path = tmp_path / 'product.xyz'
    write_xyz(symbols, [[0, 0, 0], [0.1, 0, 0]], product_path)
    reference = {'status': 'complete', 'request': {'combined_reactant_identity': {'inchi': 'R'}},
                 'fragments': [
        {'identity': {'inchi': 'A'}, 'optimized_geometry': 'a.xyz'},
        {'identity': {'inchi': 'B'}, 'optimized_geometry': 'b.xyz'}]}
    def identity(path):
        name = Path(path).name
        return {'inchi': {'product.xyz': 'P', 'fragment_1.xyz': 'A',
                          'fragment_2.xyz': 'B', 'endpoint.xyz': 'P'}.get(name, 'R')}
    monkeypatch.setattr('pyar.reaction_identity.molecule_identity_from_xyz', identity)
    monkeypatch.setattr('pyar.reaction_identity.same_molecular_identity',
                        lambda first, second: first.get('inchi') == second.get('inchi'))
    result = association.classify_association_endpoints(
        symbols, [np.array([[0., 0, 0], [.1, 0, 0]]),
                  np.array([[0., 0, 0], [second_distance, 0, 0]])],
        product_path, reference, reactants, .5, reference_directory=reference_dir)
    assert result['status'] == 'unmatched'
    assert result['branch_results'][1]['classification'] not in {
        'stable_pre_reactive_minimum', 'separated_fragments_minima'}
    assert not result['reactant_product_connection_confirmed']


def test_reactant_relaxation_failure_is_classified_without_claiming_a_connection(tmp_path, monkeypatch):
    from pyar.backends import write_xyz
    from pyar.workflows.association_reference import classify_reactant_relaxation_failure
    symbols = ['H', 'H']
    product = tmp_path / 'product.xyz'
    write_xyz(symbols, [[0, 0, 0], [0.8, 0, 0]], product)
    execution = tmp_path / 'execution'
    execution.mkdir()
    from pyar.neb import _relax_endpoint, _write_xyz_trajectory

    def optimize(symbols, coordinates, calculator, output, label, *args, **kwargs):
        frames = [np.asarray(coordinates), np.array([[0., 0, 0], [.8, 0, 0]])]
        energies = [-.5, -1.]
        _write_xyz_trajectory(output / f'{label}_path.xyz', symbols, frames, energies)
        return frames, energies, True

    monkeypatch.setattr('pyar.neb._optimize_geometry', optimize)
    with pytest.raises(RuntimeError, match='connectivity changed') as failure:
        _relax_endpoint(symbols, np.array([[0., 0, 0], [3., 0, 0]]), object(),
                        execution, 'reactant', .05, 20)
    failed = {'failed_stage': 'relax', 'error': str(failure.value)}
    monkeypatch.setattr('pyar.reaction_identity.molecule_identity_from_xyz', lambda _: {'inchi': 'product'})
    monkeypatch.setattr('pyar.reaction_identity.same_molecular_identity', lambda first, second: True)
    result = classify_reactant_relaxation_failure(failed, execution, product)
    assert result['status'] == 'product_formed_during_reactant_relaxation'
    assert result['classification'] == 'product'
    assert Path(result['trajectory']) == execution / 'reactant_path.xyz'
    assert not result['reactant_product_connection_confirmed']
    nonconverged = classify_reactant_relaxation_failure(
        {'error': 'Reactant endpoint did not converge during unbiased relaxation'}, execution, product)
    assert nonconverged['status'] == 'reactant_relaxation_nonconverged'
