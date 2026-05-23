#!/usr/bin/env python3
"""Command-line entrypoint for PyAR."""

import argparse
import datetime
import logging
import os
import sys
import time
from collections import Counter, defaultdict

from pyar.data import defualt_parameters

logger = logging.getLogger('pyar')
handler = None

BACKEND_PROFILES = {
    "gaussian": {
        "family": "dft_qc",
        "supports": {"method", "basis", "scf_cycles", "nprocs"},
    },
    "orca": {
        "family": "dft_qc",
        "supports": {"method", "basis", "scf_cycles", "nprocs"},
    },
    "psi4": {
        "family": "dft_qc",
        "supports": set(),
    },
    "turbomole": {
        "family": "dft_qc",
        "supports": {"method", "basis"},
    },
    "mopac": {
        "family": "semiempirical",
        "supports": set(),
    },
    "xtb": {
        "family": "semiempirical",
        "staged_optimization": True,
        "supports": {"opt_threshold", "nprocs"},
    },
    "xtb_turbo": {
        "family": "semiempirical",
        "staged_optimization": False,
        "supports": {"nprocs"},
    },
    "obabel": {
        "family": "semiempirical",
        "supports": set(),
    },
    "aimnet_2": {
        "family": "mlip",
        "supports": set(),
    },
    "aiqm1_mlatom": {
        "family": "mlip",
        "supports": set(),
    },
    "mlatom_aiqm1": {
        "family": "mlip",
        "supports": set(),
    },
    "xtb-aimnet2": {
        "family": "hybrid",
        "staged_optimization": True,
        "supports": {"opt_threshold", "nprocs"},
    },
    "xtb-aiqm1": {
        "family": "hybrid",
        "staged_optimization": True,
        "supports": {"opt_threshold", "nprocs"},
    },
}

QC_OPTION_ALIASES = {
    "basis": "--basis (basis set)",
    "method": "--method (functional/method)",
    "opt_threshold": "--opt-threshold",
    "opt_cycles": "--opt-cycles",
    "scf_threshold": "--scf-threshold",
    "scf_cycles": "--scf-cycles",
    "nprocs": "--nprocs",
    "custom_keywords": "--custom-keywords",
    "gamma": "--gmin/--gmax",
    "model": "--model",
}

QC_PARAMETER_KEYS = {
    "basis",
    "method",
    "opt_cycles",
    "opt_threshold",
    "scf_cycles",
    "scf_threshold",
    "nprocs",
    "gamma",
    "custom_keywords",
    "custom_keyword",
    "model",
}


def _get_file_handler():
    """Create the CLI log file handler only when the command actually runs."""
    global handler
    if handler is None:
        handler = logging.FileHandler('pyar.log', 'a')
    return handler


def _active_run_mode(run_parameters):
    """Return the active top-level workflow mode name."""
    if run_parameters['aggregate']:
        return 'aggregate'
    if run_parameters['solvate']:
        return 'solvate'
    if run_parameters['react']:
        return 'react'
    if run_parameters['scan_bond']:
        return 'scan-bond'
    return 'unknown'


def _verbosity_name(level):
    """Map CLI verbosity integer to logging level label."""
    return {
        0: 'DEBUG',
        1: 'INFO',
        2: 'WARNING',
        3: 'ERROR',
        4: 'CRITICAL',
    }.get(level, 'INFO')


def _provided_qc_options(args):
    """Return QC option keys explicitly provided by user on CLI."""
    provided = set()
    if args.get("basis") is not None:
        provided.add("basis")
    if args.get("method") is not None:
        provided.add("method")
    if args.get("model") is not None:
        provided.add("model")
    if args.get("custom_keywords") is not None:
        provided.add("custom_keywords")
    if args.get("nprocs") is not None:
        provided.add("nprocs")
    if "--opt-threshold" in sys.argv:
        provided.add("opt_threshold")
    if "--opt-cycles" in sys.argv:
        provided.add("opt_cycles")
    if "--scf-threshold" in sys.argv:
        provided.add("scf_threshold")
    if "--scf-cycles" in sys.argv:
        provided.add("scf_cycles")
    if args.get("gmin") is not None or args.get("gmax") is not None:
        provided.add("gamma")
    return provided


def _validate_backend_qc_options(software, provided_options):
    """Return backend family and unsupported provided options for a backend."""
    profile = BACKEND_PROFILES.get(software)
    if profile is None:
        return "unknown", []
    unsupported = sorted(provided_options - set(profile["supports"]))
    return profile["family"], unsupported


def _supports_staged_optimization(software):
    """Return True when the backend has wired loose/normal staged optimization."""
    profile = BACKEND_PROFILES.get(software)
    if profile is None:
        return False
    return bool(profile.get("staged_optimization", False))


def _mask_unsupported_qc_parameters(qc_params, software):
    """Return QC params with backend-unsupported options set to None."""
    masked = dict(qc_params)
    profile = BACKEND_PROFILES.get(software)
    if profile is None:
        return masked
    unsupported_options = QC_PARAMETER_KEYS - set(profile["supports"])
    for option in unsupported_options:
        if option in masked:
            masked[option] = None
        if option == "custom_keywords":
            masked["custom_keyword"] = None
    return masked


def _format_effective_qc_settings(qc_params):
    """Format active QC parameters without listing ignored/None values."""
    labels = []
    for key in sorted(QC_PARAMETER_KEYS):
        value = qc_params.get(key)
        if value is not None and key != "custom_keyword":
            labels.append(f"{key}={value}")
    return " ".join(labels) if labels else "none"


def _aggregate_stoichiometry_label(input_molecules, aggregate_sizes):
    """Return the aggregate stoichiometry label implied by the inputs."""
    parts = []
    for molecule, count in zip(input_molecules, aggregate_sizes):
        count = int(count)
        atoms_list = getattr(molecule, "atoms_list", None)
        if atoms_list:
            counts = Counter(atoms_list * count)
            fragment = []
            if 'C' in counts:
                carbon = counts.pop('C')
                fragment.append('C' if carbon == 1 else f'C{carbon}')
            if 'H' in counts:
                hydrogen = counts.pop('H')
                fragment.append('H' if hydrogen == 1 else f'H{hydrogen}')
            for element in sorted(counts):
                element_count = counts[element]
                fragment.append(element if element_count == 1 else f'{element}{element_count}')
            parts.append(''.join(fragment))
        else:
            label = getattr(molecule, "name", None) or getattr(molecule, "title", None) or "unknown"
            if count == 1:
                parts.append(label)
            elif len(label) == 1 and label.isalpha():
                parts.append(f"{label}{count}")
            else:
                parts.append(f"{label}x{count}")
    return ''.join(parts) if parts else 'unknown'


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
    from pyar.Molecule import Molecule
    from pyar import aggregator

    if os.path.exists(spec):
        return Molecule.from_xyz(spec)

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
    run_mode = _active_run_mode(run_parameters)
    logger.info(f'Run mode: {run_mode}')
    logger.info(f'Log level: {_verbosity_name(run_parameters["verbosity"])}')
    if run_parameters['formula']:
        logger.info(f'Formula input: {run_parameters["formula"]}')
    else:
        logger.info(f'Input specs: {input_files}')
    logger.info(
        "Trial sampling: Fibonacci approach directions with quaternion rotations "
        "for multi-atom monomers"
    )
    logger.debug(f'Parsed CLI options: {dict(run_parameters)}')

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
        logger.info('Multiplicity not provided: inferred defaults from electron parity.')

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
    provided_qc_options = _provided_qc_options(args)
    backend_family = "none"
    staged_optimization = False
    ignored_qc_options = []
    if run_parameters["software"] is not None:
        backend_family, ignored_qc_options = _validate_backend_qc_options(
            run_parameters["software"], provided_qc_options
        )
        staged_optimization = _supports_staged_optimization(run_parameters["software"])
        if ignored_qc_options:
            ignored_labels = [QC_OPTION_ALIASES.get(k, k) for k in ignored_qc_options]
            logger.warning(
                "Backend '%s' ignores unsupported options: %s",
                run_parameters["software"],
                ", ".join(ignored_labels),
            )
    effective_qc_options = sorted(provided_qc_options - set(ignored_qc_options))

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
    quantum_chemistry_parameters['_two_layer_optimization'] = staged_optimization
    quantum_chemistry_parameters = _mask_unsupported_qc_parameters(
        quantum_chemistry_parameters,
        run_parameters["software"],
    )

    logger.info(f'QM Software:   {quantum_chemistry_parameters["software"]}')
    logger.info(f'Backend family: {backend_family}')
    logger.info(
        "Optimization layers: %s",
        "two-stage" if staged_optimization else "single-stage",
    )
    if ignored_qc_options:
        logger.info(
            "Ignored QC options: %s",
            ", ".join(QC_OPTION_ALIASES.get(k, k) for k in ignored_qc_options),
        )
    else:
        logger.info("Ignored QC options: none")
    logger.info(
        "Effective QC options: %s",
        ", ".join(QC_OPTION_ALIASES.get(k, k) for k in effective_qc_options) if effective_qc_options else "defaults-only",
    )
    logger.info("QC settings: %s", _format_effective_qc_settings(quantum_chemistry_parameters))
    logger.debug(f'QC parameter object: {quantum_chemistry_parameters}')
    if run_parameters['software'] in {'xtb', 'xtb_turbo', 'xtb-aimnet2', 'xtb-aiqm1'}:
        if run_parameters['nprocs'] is not None:
            logger.info("xTB parallel threads: %s", run_parameters['nprocs'])
        else:
            logger.info("xTB parallel threads: not requested; xTB will run serially")

    number_of_orientations = run_parameters['how_many_orientations']
    logger.info(f'Number of orientations: {number_of_orientations}')
    maximum_number_of_seeds = run_parameters['maximum_number_of_seeds']
    logger.info(f'Maximum number of seeds: {maximum_number_of_seeds}')

    if run_parameters['site'] is None:
        site = None
    else:
        site = run_parameters['site']
        site = [site[0], input_molecules[0].number_of_atoms + site[1]]
    logger.info(f'Site constraint: {site}')
    if run_mode == 'aggregate':
        planned_aggregate_sizes = formula_aggregate_sizes if run_parameters['formula'] else run_parameters['aggregate_size']
        logger.info(
            f'Plan: aggregate fragments={len(input_molecules)} '
            f'sizes={planned_aggregate_sizes} pathways={run_parameters["number_of_pathways"]} '
            f'first_pathway={run_parameters["first_pathway"]}'
        )
    elif run_mode == 'solvate':
        logger.info(
            f'Plan: solvate seeds={max(len(input_molecules) - 1, 0)} '
            f'solvent_count={run_parameters["solvation_size"]}'
        )
    elif run_mode == 'react':
        logger.info(
            f'Plan: react gamma_range=({run_parameters["gmin"]}, {run_parameters["gmax"]}) '
            f'orientations={number_of_orientations}'
        )
    elif run_mode == 'scan-bond':
        logger.info(
            f'Plan: scan-bond pair={run_parameters["scan_bond"]} '
            f'orientations={number_of_orientations}'
        )

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
            selected_stoichiometry = _aggregate_stoichiometry_label(input_molecules, size_of_aggregate)
            logger.info(
                "Output hierarchy: aggregates/<aggregate-id>/selected/stoichiometry_%s/",
                selected_stoichiometry,
            )
            aggregator.aggregate(input_molecules, size_of_aggregate,
                                 number_of_orientations,
                                 quantum_chemistry_parameters,
                                 maximum_number_of_seeds,
                                 run_parameters['first_pathway'],
                                 run_parameters['number_of_pathways'],
                                 site)
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
                               site)
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
                          site, proximity_factor)
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
