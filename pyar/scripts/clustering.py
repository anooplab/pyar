#!/usr/bin/env python3
"""Importable ``pyar-clustering`` entrypoint.

This utility clusters or filters XYZ pools using the selection algorithms
implemented in :mod:`pyar.data_analysis.clustering`. It prints the energy
table for the input pool and then emits the selected geometries.
"""

import argparse
import sys

from pyar.data_analysis import clustering
from pyar.Molecule import Molecule


def main():
    """Cluster or filter the provided XYZ pool and print the selected files."""
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', type=str, nargs='+',
                        help="input xyz files for analysis")
    parser.add_argument('-m', '--mode', choices=['filter', 'cluster'],
                        default='cluster')
    parser.add_argument(
        '-a', '--algorithm',
        choices=['hybrid', 'hdbscan', 'optics', 'spectral', 'dbscan',
                 'kmeans', 'meanshift', 'affinity', 'agglomerative',
                 'gaussian_mixture', 'maxmin', 'max-min', 'max_min'],
        default='hybrid',
        help="selection algorithm used in cluster mode",
    )
    parser.add_argument(
        '-n', '--maximum-number-of-seeds',
        type=int,
        default=12,
        help="maximum number of geometries to keep",
    )
    args = parser.parse_args()
    input_files = args.input_files
    if len(input_files) < 2:
        print('Not enough files to cluster')
        sys.exit(0)

    mols = []
    for each_file in input_files:
        mol = Molecule.from_xyz(each_file)
        mol.energy = clustering.read_energy_from_xyz_file(each_file)
        mols.append(mol)

    clustering.print_energy_table(
        mols,
        stream=sys.stdout,
        title="Input pool energies:",
    )
    selected = []
    if args.mode == 'cluster':
        selected = clustering.choose_geometries(
            mols,
            maximum_number_of_seeds=args.maximum_number_of_seeds,
            algorithm=args.algorithm,
        )
    if args.mode == 'filter':
        selected = clustering.remove_similar(mols)
    clustering.print_energy_table(
        selected,
        stream=sys.stdout,
        title="Selected pool energies:",
    )
    print(' '.join(one.name + '.xyz' for one in selected))


if __name__ == '__main__':
    main()
