#!/usr/bin/env python3
"""Importable ``pyar-react`` entrypoint.

This script provides the historical reaction-search command-line interface
that wraps :func:`pyar.workflows.reaction.react`. It converts CLI arguments
into the backend parameter dictionary expected by the workflow and applies
the same geometry-optimizer guardrails as the main ``pyar-cli`` entrypoint.
"""

import argparse
import logging
import sys
from collections import defaultdict

from pyar.core.molecule import Molecule
from pyar.biases.softmin import resolve_softmin_beta
from pyar.biases.controller import resolve_controller_policy
from pyar.backend_capabilities import (
    backend_supports_geometry_optimization,
    normalize_backend_name,
    supported_geometry_backends,
)
from pyar.data import defualt_parameters
from pyar.state.reaction import ReactionStateError
from pyar.workflows import reaction as reaction_workflow

logger = logging.getLogger('pyar-react')
handler = logging.FileHandler('pyar-react.log', 'a')


def argument_parse():
    """Parse the reaction-search command-line arguments."""
    parser = argparse.ArgumentParser(description="pyar-react - Command-line interface for PyAR reactor")
    parser.add_argument("input_files", metavar='file', type=str, nargs='+', help='input coordinate files in xyz format.')
    parser.add_argument(
        '-N',
        dest='how_many_orientations',
        metavar='N',
        required=True,
        help='Number of trial orientations to generate.',
    )
    parser.add_argument(
        '--gmin', '--bias-min', dest='bias_min', type=float, required=True,
        help='minimum reaction-bias strength (legacy alias: --gmin)',
    )
    parser.add_argument(
        '--gmax', '--bias-max', dest='bias_max', type=float, required=True,
        help='maximum reaction-bias strength (legacy alias: --gmax)',
    )
    parser.add_argument('--bias-potential', choices=['afir', 'softmin'], default='afir')
    parser.add_argument(
        '--softmin-beta', type=resolve_softmin_beta, default=1.0,
        help='soft-min localization parameter in Bohr^-1 (default: 1.0)',
    )
    parser.add_argument(
        '--bias-controller', choices=['fixed', 'scheduled', 'adaptive'], default=None,
        help='bias force-scale controller (default: fixed)',
    )
    parser.add_argument('--bias-alpha-min', type=float, default=None,
                        help='lower bound for the applied bias scale (Ha/Bohr)')
    parser.add_argument('--bias-alpha-margin', type=float, default=None,
                        help='positive adaptive driving margin (Ha/Bohr; default: 0.001)')
    parser.add_argument('--bias-alpha-smoothing', type=float, default=None,
                        help='smoothing fraction for alpha decreases, 0 < value <= 1 (default: 1)')
    parser.add_argument('--bias-alpha-epsilon', type=float, default=None,
                        help='positive regularizer in the alpha_critical denominator')
    parser.add_argument('--bias-scheduled-alpha', type=float, default=None,
                        help='constant scale used by the scheduled controller (Ha/Bohr)')
    parser.add_argument('--software', type=str, required=True, help='Backend used to evaluate energy and forces')
    parser.add_argument('--method', default=defualt_parameters.values['method'], help='Electronic-structure method')
    parser.add_argument('--basis', default=defualt_parameters.values['basis'], help='Basis set')
    parser.add_argument('--scf-cycles', type=int, default=defualt_parameters.values['scf_cycles'], help='Maximum SCF iterations')
    parser.add_argument('--nprocs', type=int, default=defualt_parameters.values['nprocs'], help='Number of backend processes or threads')
    parser.add_argument(
        '--geometry-optimizer',
        choices=['native', 'geometric'],
        help='Optimizer for the reaction-bias objective; defaults to geometric for backends with an energy-gradient provider',
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
    """Return the zero-based index separating the two reactants in ``filename``."""
    with open(filename, 'r') as f:
        num_atoms = int(f.readline().strip())
    return num_atoms - 1


def main():
    """Run the reaction-search command-line workflow."""
    args = argument_parse()
    run_parameters = defaultdict(lambda: None, vars(args))
    try:
        controller_policy = resolve_controller_policy(
            run_parameters['bias_controller'],
            alpha_min=run_parameters['bias_alpha_min'],
            safety_margin=run_parameters['bias_alpha_margin'],
            smoothing=run_parameters['bias_alpha_smoothing'],
            epsilon=run_parameters['bias_alpha_epsilon'],
            scheduled_alpha=run_parameters['bias_scheduled_alpha'],
        )
    except ValueError as exc:
        sys.exit(str(exc))
    if run_parameters['bias_controller'] is not None or controller_policy != 'fixed':
        run_parameters['bias_controller'] = controller_policy
    try:
        softmin_beta = resolve_softmin_beta(run_parameters['softmin_beta'])
    except ValueError as exc:
        sys.exit(str(exc))

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
    software = normalize_backend_name(run_parameters['software'])
    if backend_supports_geometry_optimization(software):
        if run_parameters['opt_target'] == 'ts':
            sys.exit(
                "Transition-state optimization is reserved for a future "
                "reaction-product workflow"
            )
        if geometry_optimizer is None:
            geometry_optimizer = 'geometric'
        elif geometry_optimizer != 'geometric':
            sys.exit(
                "Reaction-bias runs with "
                f"{', '.join(supported_geometry_backends())} require "
                "--geometry-optimizer geometric"
            )
    else:
        if geometry_optimizer == 'geometric':
            sys.exit(
                f"Backend '{run_parameters['software']}' cannot be used with geomeTRIC reaction-bias "
                "optimisation because it does not expose Cartesian energy and gradients."
            )
        geometry_optimizer = geometry_optimizer or 'native'
    if controller_policy != 'fixed' and geometry_optimizer != 'geometric':
        sys.exit(
            "Non-fixed bias control requires a registered Cartesian energy-gradient backend "
            "and --geometry-optimizer geometric."
        )
    qc_params = {
        'software': run_parameters['software'],
        'index': index,
        'geometry_optimizer': geometry_optimizer,
        'opt_target': run_parameters['opt_target'],
        'bias_potential': run_parameters['bias_potential'],
        'softmin_beta': softmin_beta,
        'method': run_parameters['method'] or defualt_parameters.values['method'],
        'basis': run_parameters['basis'] or defualt_parameters.values['basis'],
        'scf_cycles': run_parameters['scf_cycles'] or defualt_parameters.values['scf_cycles'],
        'nprocs': run_parameters['nprocs'] or defualt_parameters.values['nprocs'],
    }
    for option in (
        'bias_controller', 'bias_alpha_min', 'bias_alpha_margin',
        'bias_alpha_smoothing', 'bias_alpha_epsilon', 'bias_scheduled_alpha',
    ):
        if run_parameters.get(option) is not None:
            qc_params[option] = run_parameters[option]
    try:
        reaction_workflow.react(
            input_molecules[0],
            input_molecules[1],
            run_parameters['bias_min'],
            run_parameters['bias_max'],
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
