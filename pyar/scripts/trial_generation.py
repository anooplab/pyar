#!/usr/bin/env python3
"""Command-line interface for trial geometry generation."""

import argparse
import os
import sys

import numpy as np

from pyar.Molecule import Molecule
import pyar.scan
import pyar.trial_generation as generation


def main():
    """Generate trial placements or composite configurations."""
    parser = argparse.ArgumentParser(
        description="Generate approach directions and trial molecular orientations."
    )
    parser.add_argument('-N', dest='number_of_trial_orientations', type=int,
                        required=True, metavar='n')
    parser.add_argument('-i', type=str, nargs=2, dest='file')
    parser.add_argument('--charge', type=int, default=0)
    parser.add_argument('--distance-scaling', metavar='x.x', default=1.5, type=float)
    parser.add_argument('--plot', action='store_true')
    task_group = parser.add_mutually_exclusive_group(required=False)
    task_group.add_argument('-best', type=int, nargs=2)
    task_group.add_argument('--make-configurations', action='store_true')
    task_group.add_argument('--make-composite', action='store_true')
    task_group.add_argument('-make-sphere', type=str)

    args = parser.parse_args()
    number_of_trial_orientations = args.number_of_trial_orientations
    d_scale = args.distance_scaling
    if args.file:
        seed = Molecule.from_xyz(args.file[0])
        monomer = Molecule.from_xyz(args.file[1])
    else:
        seed = None
        monomer = None

    points = generation.generate_points(number_of_trial_orientations)
    if args.plot:
        generation.plot_points(points, os.getcwd())

    if args.make_configurations:
        if not (seed and monomer):
            parser.print_help()
            sys.exit("Provide two molecules with '-i file1.xyz file2.xyz'.")
        for i, vector in enumerate(points):
            generation.merge_two_molecules(
                vector, seed, monomer, site=None, distance_scaling=d_scale
            ).mol_to_xyz(f'mol_{i:03d}.xyz')

    if args.make_composite:
        if not (seed and monomer):
            parser.print_help()
            sys.exit("Provide two molecules with '-i file1.xyz file2.xyz'.")
        for i in range(100):
            points = generation.generate_points(
                number_of_trial_orientations,
                sequence_offset=i,
            )
            result = generation.create_composite_molecule(
                seed, monomer, points, d_scale=d_scale
            )
            result.title = f"trial_{i:03d} sequence_offset={i}"
            result.mol_to_xyz(f'result_{i:03}.xyz')

    if args.make_sphere:
        import scipy.spatial.distance as sdist
        from pyar.data import new_atomic_data as atomic_data

        radii = [atomic_data.covalent_radius[n] for n in args.make_sphere.capitalize()]
        sphere_points = points[:, :3].copy()
        factor = sphere_points / 10
        while np.min(sdist.pdist(sphere_points)) < 2 * radii:
            sphere_points += factor
        atoms = [args.make_sphere for _ in range(number_of_trial_orientations)]
        result = Molecule(atoms, sphere_points)
        result.mol_to_xyz('mol.xyz')

    if args.best:
        if not (seed and monomer):
            parser.print_help()
            sys.exit("Provide two molecules with '-i file1.xyz file2.xyz'.")
        a = args.best[0]
        b = seed.number_of_atoms + args.best[1]
        pyar.scan.generate_guess_for_bonding(
            'xxx', seed, monomer, a, b, number_of_trial_orientations, d_scale=d_scale
        )


if __name__ == "__main__":
    main()
