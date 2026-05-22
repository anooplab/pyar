#!/usr/bin/env python3
"""Command-line entrypoint for PyAR."""

import argparse
import datetime
import logging
import os
import sys
import time
from collections import defaultdict

from pyar.data import defualt_parameters

logger = logging.getLogger('pyar')
handler = None


def _get_file_handler():
    """Create the CLI log file handler only when the command actually runs."""
    global handler
    if handler is None:
        handler = logging.FileHandler('pyar.log', 'a')
    return handler


def argument_parse():
    pyar_description = """pyar is a program to predict aggregation, reaction,
and clustering. Reactor explores several possible reactions between two given
molecules. Aggregator explores several possible geometries of weakly bound
molecular complexes or atomic clusters.
"""
    pyar_epilog = """Examples:
  pyar-cli -a C H -as 1 4 -N 8
  pyar-cli --aggregate --formula C5H4 -N 8
  pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000
  pyar-cli --scan-bond 1 2 A.xyz B.xyz -N 8

In aggregate mode, each input can be an XYZ file, a bare element symbol, or a
chemical formula.
"""

    parser = argparse.ArgumentParser(
        prog='pyar-cli',
        description=pyar_description,
        epilog=pyar_epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-v', '--verbosity',
                        choices=[0, 1, 2, 3, 4],
                        type=int,
                        default=1,
                        help="Choose output verbosity"
                             " (0=Debug; 1 = Info (default); "
                             "2 = Warning; 3 = Error; 4 = Critical)")
    parser.add_argument("input_files", metavar='input',
                        type=str, nargs='*',
                        help='input coordinate files in xyz format. In '
                             '--aggregate mode, bare element symbols and '
                             'chemical formulas are also accepted.')
    parser.add_argument('-N', dest='how_many_orientations', metavar='N',
                        required=True,
                        help='The number of orientations to be used')
    parser.add_argument('--tabu', choices=['y', 'n'], default='y',
                        help='Toggle Tabu search algorithm. Default is on (y)')
    parser.add_argument('--grid', choices=['y', 'n'], default='y',
                        help='Toggle the use of grid for search space.')
    parser.add_argument('--formula', type=str,
                        help='Chemical formula of the molecule to generate '
                             'in --aggregate mode')
    parser.add_argument('-nprocs', '--nprocs', metavar='n',
                        type=int, help='The number of processors/cores to be '
                                       'used by the quantum chemistry software.'
                                       'Not fully implemented with all '
                                       'interfaces. Write to '
                                       'anoop@chem.iitkgp.ac.in'
                                       ' if this does not work properly')
    parser.add_argument('-model', '--model', metavar='model',
                        type=str, help='The model to be used for the '
                                       'aggregation. Default is '
                                       'aimnet2_wb97m-d3_ens.jpt')
    parser.add_argument('-basis', '--basis', type=str,
                        help='Basis set (default=def2-SVP)')
    parser.add_argument('-method', '--method', type=str,
                        help='The method (default=BP86)')

    run_type_group = parser.add_mutually_exclusive_group(required=True)
    run_type_group.add_argument("-r", "--react", action='store_true',
                                help="Run a reactor calculation")
    run_type_group.add_argument("-s", "--solvate", action='store_true',
                                help="Add one solvent molecules to given solute molecules")
    run_type_group.add_argument("-a", "--aggregate", action='store_true',
                                help="Run a aggregator calculation")
    run_type_group.add_argument("--scan-bond", nargs=2, type=int,
                                metavar=('a', 'b'),
                                help="scan a bond between the given atoms of two fragments")

    parser.add_argument('-ss', '--solvation-size', type=int, metavar='n',
                        help='number of solvent molecules to be added')
    parser.add_argument('-mns', '--maximum-number-of-seeds', metavar='n',
                        type=int, help='maximum number of seeds')
    parser.add_argument('-f', '--features',
                        choices=['fingerprint', 'scm', 'moi', 'fsmd', 'soap', 'mbtr',
                                 'ani', 'lmbtr', 'acsf', 'sinematrix', 'vallornav'],
                        default='fingerprint',
                        help="Choose the features to be used for clustering")
    parser.add_argument('-as', '--aggregate-size', type=int, nargs='*',
                        metavar=('l', 'm',),
                        help='number of monomers in aggregate; defaults to '
                             '[1] for single-formula aggregate runs')
    parser.add_argument('--number-of-pathways', type=int, metavar='n',
                        help='How many pathways to be used in binary/ternary aggregation.')
    parser.add_argument('-gmin', type=float, help='minimum value of gamma')
    parser.add_argument('-gmax', type=float, help='maximum value of gamma')
    parser.add_argument('--site', type=int, nargs=2,
                        help='atom for site specific reaction')
    parser.add_argument("-c", "--charge", type=int, nargs='+', metavar='c',
                        help="Charge of the system")
    parser.add_argument("-m", "--multiplicity", type=int, nargs='+', metavar='m',
                        help="Multiplicity of the system")
    parser.add_argument("--scftype", type=str, nargs='+',
                        help="specify rhf or uhf (default=rhf)")
    parser.add_argument("--software", type=str,
                        choices=['gaussian', 'mopac', 'obabel', 'orca',
                                 'psi4', 'turbomole', 'xtb', 'xtb_turbo',
                                 'mlatom_aiqm1', 'aimnet_2', 'aiqm1_mlatom',
                                 'xtb-aimnet2', 'xtb-aiqm1'],
                        required=False, default=None, help="Software")
    parser.add_argument('--opt-threshold', type=str, default='normal',
                        choices=['loose', 'normal', 'tight'],
                        help='Optimization threshold')
    parser.add_argument('--opt-cycles', type=int, default=100,
                        help='Maximum optimization cycles')
    parser.add_argument('--scf-threshold', type=str, default='normal',
                        choices=['loose', 'normal', 'tight'],
                        help='SCF threshold')
    parser.add_argument('--scf-cycles', type=int, default=1000,
                        help='Maximum SCF cycles.')
    parser.add_argument('--custom-keywords', type=str,
                        help='Software related custom keywords.')
    return parser.parse_args()


def _resolve_aggregate_input(spec, aggregate_mode):
    """Resolve an input spec to a molecule, allowing formulas in aggregate mode."""
    from pyar import Molecule, aggregator

    if os.path.exists(spec):
        return Molecule.Molecule.from_xyz(spec)

    if aggregate_mode:
        if spec.lower().endswith(".xyz") or os.path.sep in spec or (
            os.path.altsep is not None and os.path.altsep in spec
        ):
            raise SystemExit(f"File {spec} does not exist")
        try:
            return aggregator.generate_molecule_from_formula(spec)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc

    raise SystemExit(f"File {spec} does not exist")


def _expand_formula_inputs(formula):
    """Expand a formula into aggregate fragment specs and multiplicities."""
    from pyar import aggregator

    try:
        return aggregator.expand_formula_to_aggregate_inputs(formula)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc


def _normalize_parameter_list(values, default_value, number_of_inputs, label):
    """Expand a single provided value or validate a full per-input list."""
    resolved = values or [default_value for _ in range(number_of_inputs)]
    if len(resolved) == 1 and number_of_inputs > 1:
        return resolved * number_of_inputs
    if len(resolved) != number_of_inputs:
        sys.exit(f'{label} are not specified for all input files')
    return resolved


def main():
    args = vars(argument_parse())
    from pyar import aggregator, reactor, scan

    run_parameters = defaultdict(lambda: None, defualt_parameters.values)

    for key, value in args.items():
        if value is not None and run_parameters[key] != value:
            run_parameters[key] = value

    if run_parameters['verbosity'] == 0:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(name)-12s %(filename)s %(funcName)s '
                                      '%(lineno)d %(levelname)-8s: %(message)s')
    elif run_parameters['verbosity'] == 1:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.INFO)
    elif run_parameters['verbosity'] == 2:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.WARNING)
    elif run_parameters['verbosity'] == 3:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.ERROR)
    else:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.CRITICAL)

    file_handler = _get_file_handler()
    file_handler.setFormatter(formatter)
    if file_handler not in logger.handlers:
        logger.addHandler(file_handler)

    time_now = datetime.datetime.now().strftime("%d %b %Y, %H:%M:%S")
    logger.info(
        r"""
+-+-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
|  _ \ _   _   / \  |  _ \
| |_) | | | | / _ \ | |_) |
|  __/| |_| |/ ___ \|  _ <
|_|    \__, /_/   \_\_| \_\
       |___/
                 Aggregation and Reaction
+-+-+-+-+ +-+-+-+-+-+-+-+-+-+ +-+-+-+ +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
"""
    )
    logger.info(
        f'==============================Starting at {time_now}==============================')
    logger.info(f'Job directory: {os.getcwd()}')
    logger.debug(f'Logging level is {{{logger.level}}}')

    input_files = run_parameters['input_files'] or []
    number_of_input_files = len(input_files)
    logger.debug(f"{number_of_input_files} input files")

    if run_parameters['formula']:
        if not run_parameters['aggregate']:
            message = '--formula can only be used with --aggregate in pyar-cli'
            logger.critical(message)
            sys.exit(message)
        if number_of_input_files != 0:
            message = 'Do not provide XYZ input files when using --formula'
            logger.critical(message)
            sys.exit(message)
    elif number_of_input_files == 0:
        message = 'Provide at least one XYZ input file, or use --formula with --aggregate'
        logger.critical(message)
        sys.exit(message)

    if run_parameters['react'] and number_of_input_files != 2:
        message = 'Reactor requires exactly two XYZ input files'
        logger.critical(message)
        sys.exit(message)
    if run_parameters['scan_bond'] and number_of_input_files != 2:
        message = 'Bond scanning requires exactly two XYZ input files'
        logger.critical(message)
        sys.exit(message)
    if run_parameters['solvate'] and number_of_input_files < 2:
        message = 'Solvation requires at least two XYZ input files'
        logger.critical(message)
        sys.exit(message)

    if run_parameters['formula']:
        input_specs, formula_aggregate_sizes = _expand_formula_inputs(run_parameters['formula'])
    else:
        input_specs = input_files
        formula_aggregate_sizes = None
    number_of_inputs = len(input_specs)

    charges = _normalize_parameter_list(run_parameters['charge'], 0, number_of_inputs, 'Charges')
    multiplicities = _normalize_parameter_list(run_parameters['multiplicity'], 1, number_of_inputs, 'Multiplicities')
    scftypes = _normalize_parameter_list(run_parameters['scftype'], 'rhf', number_of_inputs, 'SCF Types')

    logger.info("Parsing the following inputs: ")
    input_molecules = []
    for each_file, charge, multiplicity, scftype in zip(input_specs, charges, multiplicities, scftypes):
        try:
            mol = _resolve_aggregate_input(each_file, run_parameters['aggregate'])
            mol.charge = charge
            mol.multiplicity = multiplicity
            input_molecules.append(mol)
            logger.info(f" {each_file} {charge} {multiplicity} {scftype}")
        except SystemExit:
            raise

    if run_parameters['multiplicity'] is None:
        multiplicities = []
        for mol, charge in zip(input_molecules, charges):
            n_electrons = sum(mol.atomic_number) - charge
            multiplicity = 1 if n_electrons % 2 == 0 else 2
            multiplicities.append(multiplicity)
            mol.multiplicity = multiplicity

    if run_parameters['software'] is not None:
        for mol in input_molecules:
            n_electrons = sum(mol.atomic_number) - mol.charge
            if n_electrons % 2 == 0:
                if mol.multiplicity % 2 != 1:
                    sys.exit(f"{n_electrons} (even) electrons and multiplicty {mol.multiplicity} (odd) is not pssible for {mol.name}")
            else:
                if mol.multiplicity % 2 == 1:
                    sys.exit(f"{n_electrons} (odd) electrons and multiplicty {mol.multiplicity} (even) is not pssible for {mol.name}")

    custom_keywords = run_parameters['custom_keywords']
    quantum_chemistry_parameters = {
        'basis': run_parameters['basis'],
        'method': run_parameters['method'],
        'software': run_parameters['software'],
        'opt_cycles': run_parameters['opt_cycles'],
        'opt_threshold': run_parameters['opt_threshold'],
        'scf_cycles': run_parameters['scf_cycles'],
        'scf_threshold': run_parameters['scf_threshold'],
        'nprocs': run_parameters['nprocs'],
        'gamma': run_parameters['gamma'],
        'custom_keywords': custom_keywords,
        'custom_keyword': custom_keywords,
        'model': run_parameters['model']
    }

    logger.info(f'QM Software:   {quantum_chemistry_parameters["software"]}')

    number_of_orientations = run_parameters['how_many_orientations']
    logger.info(f'Number of orientations: {number_of_orientations}')
    maximum_number_of_seeds = run_parameters['maximum_number_of_seeds']
    logger.info(f'Maximum number of seeds: {maximum_number_of_seeds}')

    if run_parameters['site'] is None:
        site = None
    else:
        site = run_parameters['site']
        site = [site[0], input_molecules[0].number_of_atoms + site[1]]

    try:
        if run_parameters['aggregate']:
            size_of_aggregate = run_parameters['aggregate_size']
            if run_parameters['formula']:
                size_of_aggregate = formula_aggregate_sizes
            if size_of_aggregate is None or len(size_of_aggregate) != len(input_molecules):
                message = ('Error: For an Aggregation run, specify \nthe desired number of each monomers to be added \nusing the argument\n -as <int> <int> ...')
                logger.critical(message)
                sys.exit(message)
            if quantum_chemistry_parameters["software"] is None:
                logger.info(
                    "No --software specified: aggregate mode will generate trial "
                    "geometries only; no quantum-chemistry optimization will be run."
                )
            t1_0 = time.time()
            time_started = datetime.datetime.now()
            aggregator.aggregate(input_molecules, size_of_aggregate,
                                 number_of_orientations,
                                 quantum_chemistry_parameters,
                                 maximum_number_of_seeds,
                                 run_parameters['first_pathway'],
                                 run_parameters['number_of_pathways'],
                                 run_parameters['tabu'] == 'y',
                                 run_parameters['grid'] == 'y', site)
            logger.info('Total Time: {}'.format(time.time() - t1_0))
            logger.info("Started at {}\nEnded at {}".format(time_started, datetime.datetime.now()))

        if run_parameters['solvate']:
            number_of_solvent_molecules = run_parameters['solvation_size']
            if number_of_solvent_molecules is None:
                sys.exit('For this please provide the number of solvent\nmolecules to be added. Use the following option\n  -ss <int>')
            if len(input_molecules) == 1:
                sys.exit('Please provide more than two molecules.\nThe last input file will be considered as solvent\nand the other molecules as solutes to which solvent\nmolecules will be added.')
            monomer = input_molecules[-1]
            seeds = input_molecules[:-1]
            t1_0 = time.time()
            time_started = datetime.datetime.now()
            aggregator.solvate(seeds, monomer,
                               number_of_solvent_molecules,
                               number_of_orientations,
                               quantum_chemistry_parameters,
                               maximum_number_of_seeds,
                               run_parameters['tabu'] == 'y',
                               run_parameters['grid'] == 'y', site)
            logger.info('Total Time: {}'.format(time.time() - t1_0))
            logger.info("Started at {}\nEnded at {}".format(time_started, datetime.datetime.now()))

        if run_parameters['react']:
            minimum_gamma = run_parameters['gmin']
            maximum_gamma = run_parameters['gmax']
            if len(input_molecules) < 2:
                sys.exit('Missing arguments: provide at least two molecules')
            if minimum_gamma is None or maximum_gamma is None:
                sys.exit('missing arguments: -gmin <integer> -gmax <integer>')
            if number_of_orientations is None:
                sys.exit('Missing arguments: -N #')

            proximity_factor = 2.3
            zero_time = time.time()
            time_started = datetime.datetime.now()
            reactor.react(input_molecules[0], input_molecules[1],
                          minimum_gamma, maximum_gamma,
                          int(number_of_orientations),
                          quantum_chemistry_parameters,
                          site, proximity_factor,
                          run_parameters['tabu'] == 'y', run_parameters['grid'] == 'y')
            logger.info('Total run time: {}'.format(time.time() - zero_time))
            logger.info(f"Started at {time_started}\nEnded at {datetime.datetime.now()}")
            return

        if run_parameters['scan_bond']:
            if number_of_orientations is None:
                sys.exit('Missing arguments: -N #')
            scan.scan_distance(input_molecules, run_parameters['scan_bond'],
                               int(number_of_orientations),
                               quantum_chemistry_parameters)
    except FileNotFoundError as exc:
        logger.critical(str(exc))
        sys.exit(str(exc))


if __name__ == "__main__":
    main()
