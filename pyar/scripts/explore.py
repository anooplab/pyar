#!/usr/bin/env python3
"""Importable `pyar-explore` entrypoint."""

import argparse
import copy
import logging
import random
from collections import OrderedDict

import numpy as np

from pyar.Molecule import Molecule
from pyar.data import new_atomic_data as atomic_data
from pyar.tabu import generate_points, merge_two_molecules

logging.basicConfig(filename='pyar_explore.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate geometries based on formula.')
    parser.add_argument('seed_file', type=str, help='XYZ file for the seed monomer.')
    parser.add_argument('monomer_file', type=str, help='XYZ file for the monomer.')
    parser.add_argument('--formula', type=str, required=True, help='Target formula, e.g., C5H4.')
    parser.add_argument('--pop', type=int, required=True, help='Number of geometries to generate.')
    return parser.parse_args()


def parse_formula(formula):
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


def create_composite_molecule_wrapper(seed, monomer, pathway):
    composite = copy.deepcopy(seed)
    points = generate_points(len(pathway), True, 0.95)
    logging.info(f"Starting pathway: {pathway}")
    for atom, point in zip(pathway, points):
        if atom == monomer.atoms_list[0]:
            new_monomer = copy.deepcopy(monomer)
        else:
            new_monomer = Molecule([atom], np.array([[0.0, 0.0, 0.0]], dtype=np.float64))
        composite = merge_two_molecules(np.array(point, dtype=np.float64), composite, new_monomer,
                                        freeze_fragments=False, distance_scaling=1.5)
    return composite


def main():
    args = parse_arguments()
    seed = Molecule.from_xyz(args.seed_file)
    monomer = Molecule.from_xyz(args.monomer_file)
    atom_counts = parse_formula(args.formula)
    base_pathway = generate_chemical_pathway(atom_counts, seed)
    geometries = []
    for _ in range(args.pop):
        pathway = copy.copy(base_pathway)
        random.shuffle(pathway)
        try:
            geometries.append(create_composite_molecule_wrapper(seed, monomer, pathway))
        except Exception as e:
            logging.error(str(e))
    for i, geometry in enumerate(geometries):
        geometry.mol_to_xyz(f'geometry_{i:03d}.xyz')
        print(f'Generated geometry saved to geometry_{i:03d}.xyz. Atom count: {len(geometry.atoms_list)}')


if __name__ == '__main__':
    main()
