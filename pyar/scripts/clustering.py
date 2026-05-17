#!/usr/bin/env python3
"""Importable `pyar-clustering` entrypoint."""

import argparse
import sys

from pyar.data_analysis import clustering
from pyar.Molecule import Molecule


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', type=str, nargs='+',
                        help="input xyz files for analysis")
    parser.add_argument('-m', '--mode', choices=['filter', 'cluster'],
                        default='cluster')
    parser.add_argument('-cf', '--clustering_features',
                        choices=['fingerprint', 'scm', 'moi', 'fsmd', 'soap',
                                 'mbtr', 'ani', 'lmbtr', 'acsf', 'sinematrix',
                                 'vallornav'],
                        default='fingerprint')
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

    clustering.plot_energy_histogram(mols)
    selected = []
    if args.mode == 'cluster':
        selected = clustering.choose_geometries(mols, features=args.clustering_features)
    if args.mode == 'filter':
        selected = clustering.remove_similar(mols)
    print(' '.join(one.name + '.xyz' for one in selected))


if __name__ == '__main__':
    main()
