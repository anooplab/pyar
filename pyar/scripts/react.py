#!/usr/bin/env python3
"""Importable `pyar-react` entrypoint."""

import argparse
import logging
import sys
from collections import defaultdict

from pyar import reactor
from pyar.Molecule import Molecule
from pyar.reaction_state import ReactionStateError

logger = logging.getLogger('pyar-react')
handler = logging.FileHandler('pyar-react.log', 'a')


def argument_parse():
    parser = argparse.ArgumentParser(description="pyar-react - Command-line interface for PyAR reactor")
    parser.add_argument("input_files", metavar='file', type=str, nargs='+', help='input coordinate files in xyz format.')
    parser.add_argument('-N', dest='how_many_orientations', metavar='N', required=True, help='The number of orientations to be used')
    parser.add_argument('--gmin', type=float, required=True, help='minimum value of gamma')
    parser.add_argument('--gmax', type=float, required=True, help='maximum value of gamma')
    parser.add_argument('--software', type=str, required=True, help='Software for optimization e.g., orca-aiqm1')
    parser.add_argument('--index', type=int, help='Index for splitting the molecule (final index of the first reactant, i.e., number of atoms - 1). If not provided, it will be calculated from the first input file.')
    return parser.parse_args()


def calculate_index_from_xyz(filename):
    with open(filename, 'r') as f:
        num_atoms = int(f.readline().strip())
    return num_atoms - 1


def main():
    args = argument_parse()
    run_parameters = defaultdict(lambda: None, vars(args))

    input_molecules = []
    for file in run_parameters['input_files']:
        try:
            mol = Molecule.from_xyz(file)
            input_molecules.append(mol)
        except IOError:
            logger.critical(f"File {file} does not exist")
            sys.exit()

    if run_parameters['index'] is None:
        index = calculate_index_from_xyz(run_parameters['input_files'][0])
    else:
        index = run_parameters['index']

    logger.info(f"Using index: {index} (final index of the first reactant)")
    qc_params = {'software': run_parameters['software'], 'index': index}
    try:
        reactor.react(input_molecules[0], input_molecules[1],
                      run_parameters['gmin'], run_parameters['gmax'],
                      int(run_parameters['how_many_orientations']), qc_params,
                      None, 2.3)
    except (FileNotFoundError, ReactionStateError) as exc:
        logger.critical(str(exc))
        sys.exit(str(exc))


if __name__ == "__main__":
    main()
