#!/usr/bin/env python3
"""Importable `pyar-react` entrypoint."""

import argparse
import logging
import sys
from collections import defaultdict

from pyar.Molecule import Molecule
from pyar.backend_capabilities import (
    backend_supports_geometry_optimization,
    supported_geometry_backends,
)
from pyar.data import defualt_parameters
from pyar.state.reaction import ReactionStateError
from pyar.workflows import reaction as reaction_workflow

logger = logging.getLogger('pyar-react')
handler = logging.FileHandler('pyar-react.log', 'a')


def argument_parse():
    parser = argparse.ArgumentParser(description="pyar-react - Command-line interface for PyAR reactor")
    parser.add_argument("input_files", metavar='file', type=str, nargs='+', help='input coordinate files in xyz format.')
    parser.add_argument('-N', dest='how_many_orientations', metavar='N', required=True, help='The number of orientations to be used')
    parser.add_argument('--gmin', type=float, required=True, help='minimum value of gamma')
    parser.add_argument('--gmax', type=float, required=True, help='maximum value of gamma')
    parser.add_argument('--software', type=str, required=True, help='Backend used to evaluate energy and forces')
    parser.add_argument('--method', default=defualt_parameters.values['method'], help='Electronic-structure method')
    parser.add_argument('--basis', default=defualt_parameters.values['basis'], help='Basis set')
    parser.add_argument('--scf-cycles', type=int, default=defualt_parameters.values['scf_cycles'], help='Maximum SCF iterations')
    parser.add_argument('--nprocs', type=int, default=defualt_parameters.values['nprocs'], help='Number of backend processes or threads')
    parser.add_argument(
        '--geometry-optimizer',
        choices=['native', 'geometric'],
        help='Optimizer for the AFIR objective; defaults to geometric for backends with an energy-gradient provider',
    )
    parser.add_argument(
        '--opt-target',
        choices=['minimum', 'ts'],
        default='minimum',
        help='Optimization target; ts is reserved for future automatic TS searches',
    )
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
    geometry_optimizer = run_parameters['geometry_optimizer']
    if backend_supports_geometry_optimization(run_parameters['software']):
        if run_parameters['opt_target'] == 'ts':
            sys.exit(
                "Transition-state optimization is reserved for a future "
                "reaction-product workflow"
            )
        if geometry_optimizer is None:
            geometry_optimizer = 'geometric'
        elif geometry_optimizer != 'geometric':
            sys.exit(
                "AFIR reaction runs with "
                f"{', '.join(supported_geometry_backends())} require "
                "--geometry-optimizer geometric"
            )
    else:
        if geometry_optimizer == 'geometric':
            sys.exit(
                f"Backend '{run_parameters['software']}' cannot be used with geomeTRIC AFIR "
                "optimisation because it does not expose Cartesian energy and gradients."
            )
        geometry_optimizer = geometry_optimizer or 'native'
    qc_params = {
        'software': run_parameters['software'],
        'index': index,
        'geometry_optimizer': geometry_optimizer,
        'opt_target': run_parameters['opt_target'],
        'method': run_parameters['method'] or defualt_parameters.values['method'],
        'basis': run_parameters['basis'] or defualt_parameters.values['basis'],
        'scf_cycles': run_parameters['scf_cycles'] or defualt_parameters.values['scf_cycles'],
        'nprocs': run_parameters['nprocs'] or defualt_parameters.values['nprocs'],
    }
    try:
        reaction_workflow.react(
            input_molecules[0],
            input_molecules[1],
            run_parameters['gmin'],
            run_parameters['gmax'],
            int(run_parameters['how_many_orientations']),
            qc_params,
            None,
            2.3,
        )
    except (FileNotFoundError, ReactionStateError, ValueError) as exc:
        logger.critical(str(exc))
        sys.exit(str(exc))


if __name__ == "__main__":
    main()
