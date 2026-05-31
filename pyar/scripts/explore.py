#!/usr/bin/env python3
"""Importable ``pyar-explore`` entrypoint.

This script generates candidate composite geometries by combining a seed
molecule, a monomer, and a formula-derived pathway. It is a convenience
driver around the trial-generation and molecule-merging helpers used by the
aggregation workflow.
"""

import argparse
from copy import deepcopy
import logging
import random
from collections import OrderedDict

import numpy as np

from pyar.core.molecule import Molecule
from pyar.data import new_atomic_data as atomic_data
from pyar.sampling.trial_generator import generate_points, merge_two_molecules

logging.basicConfig(filename='pyar_explore.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def parse_arguments():
    """Parse the `pyar-explore` command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate geometries based on formula.')
    parser.add_argument('seed_file', type=str, help='XYZ file for the seed monomer.')
    parser.add_argument('monomer_file', type=str, help='XYZ file for the monomer.')
    parser.add_argument('--formula', type=str, required=True, help='Target formula, e.g., C5H4.')
    parser.add_argument('--pop', type=int, required=True, help='Number of geometries to generate.')
    return parser.parse_args()


def parse_formula(formula):
    """Parse a chemical formula into an ordered atom-count mapping."""
    atom_counts = OrderedDict()
    current_atom = ''
    for char in formula:
        if char.isalpha():
            if current_atom:
                atom_counts[current_atom] = atom_counts.get(current_atom, 0) + 1
            current_atom = char
        else:
            atom_counts[current_atom] = atom_counts.get(current_atom, 0) + int(char)
    if current_atom:
        atom_counts[current_atom] = atom_counts.get(current_atom, 0) + 1
    return atom_counts


def generate_chemical_pathway(atom_counts, seed):
    """Return the list of atoms that still need to be added to the seed."""
    pathway = []
    seed_atom_counts = {}
    for atom in seed.atoms_list:
        seed_atom_counts[atom] = seed_atom_counts.get(atom, 0) + 1
    sorted_atoms = sorted(atom_counts.items(), key=lambda x: atomic_data.atomic_number[x[0]], reverse=True)
    for atom, count in sorted_atoms:
        remaining_count = count - seed_atom_counts.get(atom, 0)
        if remaining_count > 0:
            pathway.extend([atom] * remaining_count)
    return pathway


def _copy_molecule(molecule):
    """Return a copied molecule-like object."""
    if hasattr(molecule, "copy"):
        return molecule.copy()
    return deepcopy(molecule)


def create_composite_molecule_wrapper(seed, monomer, pathway, sequence_offset=0):
    """Create one candidate composite geometry for a pathway.

    The wrapper seeds the orientation generator with ``sequence_offset`` so
    repeated populations are reproducible while still producing distinct
    geometry sets.
    """
    composite = _copy_molecule(seed)
    points = generate_points(len(pathway), sequence_offset=sequence_offset)
    logging.info("Starting pathway with sequence offset %d: %s", sequence_offset, pathway)
    for atom, point in zip(pathway, points):
        if atom == monomer.atoms_list[0]:
            new_monomer = _copy_molecule(monomer)
        else:
            new_monomer = Molecule([atom], np.array([[0.0, 0.0, 0.0]], dtype=np.float64))
        composite = merge_two_molecules(np.array(point, dtype=np.float64), composite, new_monomer,
                                        freeze_fragments=False, distance_scaling=1.5)
    return composite


def main():
    """Generate one or more composite geometries from a formula."""
    args = parse_arguments()
    seed = Molecule.from_xyz(args.seed_file)
    monomer = Molecule.from_xyz(args.monomer_file)
    atom_counts = parse_formula(args.formula)
    base_pathway = generate_chemical_pathway(atom_counts, seed)
    geometries = []
    for population_index in range(args.pop):
        pathway = base_pathway.copy()
        random.shuffle(pathway)
        try:
            geometries.append(
                create_composite_molecule_wrapper(
                    seed,
                    monomer,
                    pathway,
                    sequence_offset=population_index,
                )
            )
        except Exception as e:
            logging.error(str(e))
    for i, geometry in enumerate(geometries):
        geometry.mol_to_xyz(f'geometry_{i:03d}.xyz')
        print(f'Generated geometry saved to geometry_{i:03d}.xyz. Atom count: {len(geometry.atoms_list)}')


if __name__ == '__main__':
    main()
