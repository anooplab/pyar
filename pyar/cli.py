#!/usr/bin/env python3
"""Command-line entrypoint for PyAR."""

import argparse
import datetime
import importlib.util
import logging
import os
import shlex
import shutil
import sys
import time
from collections import Counter, defaultdict

from pyar.data import defualt_parameters
from pyar.backend_capabilities import (
    backend_family,
    backend_supports_geometry_optimization,
    backend_supports_staged_optimization,
    get_backend_capabilities,
    normalize_backend_name,
    unsupported_qc_options,
    supported_geometry_backends,
)

logger = logging.getLogger('pyar')
handler = None

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

_SUBCOMMAND_ALIASES = {
    "aggregate": "--aggregate",
    "react": "--react",
    "solvate": "--solvate",
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
    return backend_family(software), unsupported_qc_options(software, provided_options)


def _configure_reaction_optimizer(run_parameters, run_mode):
    """Select the external optimizer for supported AFIR reaction backends."""
    if run_mode != "react" or not backend_supports_geometry_optimization(run_parameters["software"]):
        return
    if run_parameters["opt_target"] == "ts":
        sys.exit(
            "Transition-state optimization is reserved for a future "
            "reaction-product workflow"
        )
    explicitly_selected = "--geometry-optimizer" in sys.argv
    if not explicitly_selected:
        run_parameters["geometry_optimizer"] = "geometric"
        logger.info(
            "Reaction optimizer: selected geomeTRIC for backend energy/forces plus AFIR bias."
        )
        return
    if run_parameters["geometry_optimizer"] != "geometric":
        sys.exit(
            "AFIR reaction runs with "
            f"{', '.join(supported_geometry_backends())} require "
            "--geometry-optimizer geometric"
        )


def _mask_unsupported_qc_parameters(qc_params, software):
    """Return QC params with backend-unsupported options set to None."""
    masked = dict(qc_params)
    for option in unsupported_qc_options(software, QC_PARAMETER_KEYS):
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


def _normalize_subcommand_argv(argv):
    """Translate subcommand style invocation into the existing flag-based CLI."""
    if not argv:
        return argv
    alias = _SUBCOMMAND_ALIASES.get(argv[0])
    if alias is None:
        return argv
    return [alias, *argv[1:]]


def _dispatch_trace_subcommand(argv):
    """Run the reaction-trace CLI from pyar-cli."""
    from pyar.scripts.reaction_trace import main as reaction_trace_main

    reaction_trace_main(argv)


def _dispatch_conformer_subcommand(argv):
    """Run the conformer-search CLI from pyar-cli."""
    from pyar.scripts.conformer import main as conformer_main

    conformer_main(argv)


def _log_workflow_result(result):
    """Log a structured workflow result when a workflow returns one."""
    if hasattr(result, "to_dict"):
        logger.info("Workflow result: %s", result.to_dict())


def argument_parse(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    else:
        argv = list(argv)
    argv = _normalize_subcommand_argv(argv)
    pyar_description = """pyar is a program to predict aggregation, reaction,
and clustering. Reactor explores several possible reactions between two given
molecules. Aggregator explores several possible geometries of weakly bound
molecular complexes or atomic clusters.
"""
    pyar_epilog = """Examples:
  pyar-cli aggregate C H -as 1 4 -N 8
  pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000
  pyar-cli solvate solute.xyz solvent.xyz -ss 10 -N 16
  pyar-cli -a C H -as 1 4 -N 8
  pyar-cli --aggregate --formula C5H4 -N 8
  pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000
  pyar-cli conformer "CCO" --num-conformers 50 --top-n 5
  pyar-cli trace .
  pyar-cli trace . --plot
  pyar-cli trace . --plot-only

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

    parser.add_argument('-ss', '--solvation-size', type=int, metavar='n',
                        help='number of solvent molecules to be added')
    parser.add_argument('-mns', '--maximum-number-of-seeds', metavar='n',
                        type=int, help='maximum number of seeds')
    parser.add_argument('-f', '--features',
                        choices=['fingerprint', 'scm', 'moi', 'fsmd', 'soap', 'mbtr',
                                 'ani', 'lmbtr', 'acsf', 'sinematrix', 'vallornav'],
                        default='fingerprint',
                        help="Choose the features to be used for clustering")
    parser.add_argument(
        "--connectivity-policy",
        choices=["auto", "off", "prefer", "strict"],
        default=defualt_parameters.values["connectivity_policy"],
        help=(
            "Connectivity handling for covalent-graph selection: "
            "auto lets the workflow decide, off never filters by connectivity, "
            "prefer keeps connected candidates when available, and strict "
            "discards disconnected candidates."
        ),
    )
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
    parser.add_argument('--geometry-optimizer', type=str, default='native',
                        choices=['native', 'geometric'],
                        help='Geometry optimizer used for backend gradients')
    parser.add_argument('--opt-target', type=str, default='minimum',
                        choices=['minimum', 'ts'],
                        help='Optimization target for the geometry optimizer')
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
    return parser.parse_args(argv)


def _resolve_aggregate_input(spec, aggregate_mode):
    """Resolve an input spec to a molecule, allowing formulas in aggregate mode."""
    from pyar.core.molecule import Molecule
    from pyar.workflows.aggregate import generate_molecule_from_formula

    if os.path.exists(spec):
        return Molecule.from_xyz(spec)

    if aggregate_mode:
        if spec.lower().endswith(".xyz") or os.path.sep in spec or (
            os.path.altsep is not None and os.path.altsep in spec
        ):
            raise SystemExit(f"File {spec} does not exist")
        try:
            return generate_molecule_from_formula(spec)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc

    raise SystemExit(f"File {spec} does not exist")


def _expand_formula_inputs(formula):
    """Expand a formula into aggregate fragment specs and multiplicities."""
    from pyar.workflows.aggregate import expand_formula_to_aggregate_inputs

    try:
        return expand_formula_to_aggregate_inputs(formula)
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


def _validate_cli_request_shape(run_parameters, number_of_input_files):
    """Validate CLI-only constraints before importing workflow modules."""
    if run_parameters['formula']:
        if not run_parameters['aggregate']:
            sys.exit('--formula can only be used with --aggregate in pyar-cli')
        if number_of_input_files != 0:
            sys.exit('Do not provide XYZ input files when using --formula')
    elif number_of_input_files == 0:
        sys.exit('Provide at least one XYZ input file, or use --formula with --aggregate')

    if run_parameters['react'] and number_of_input_files != 2:
        sys.exit('Reactor requires exactly two XYZ input files')
    if run_parameters['solvate'] and number_of_input_files < 2:
        sys.exit('Solvation requires at least two XYZ input files')


def _missing_module_message(module_name, install_hint):
    """Return a short install hint for a missing Python module."""
    if install_hint.startswith("http"):
        install_text = f"See {install_hint} for installation details."
    else:
        install_text = f"Install it with {install_hint}."
    return (
        f"Missing Python package '{module_name}'. "
        f"{install_text}"
    )


def _missing_executable_message(program, friendly_name, install_hint):
    """Return a short install hint for a missing executable."""
    if install_hint.startswith("http"):
        install_text = f"See {install_hint} for installation details."
    else:
        install_text = f"Install it with {install_hint}."
    return (
        f"Missing executable '{program}' for {friendly_name}. "
        f"{install_text} Make sure '{program}' is on PATH."
    )


_BACKEND_EXECUTABLE_REQUIREMENT_HINTS = {
    "g16": ("Gaussian", "https://www.gaussian.com/g16/g16src_install.pdf"),
    "orca": ("ORCA", "https://orca-manual.mpi-muelheim.mpg.de/contents/quickstartguide/installation.html"),
    "psi4": ("Psi4", "https://psicode.org/psi4manual/master/index.html"),
    "define": ("Turbomole", "https://www.turbomole.org/turbomole/turbomole-documentation/"),
    "xtb": ("xTB", "https://xtb-docs.readthedocs.io/en/latest/setup.html"),
    "mopac": ("MOPAC", "https://openmopac.net/download/installer/"),
    "obabel": ("OpenBabel", "https://openbabel.org/docs/Installation/install.html"),
    "obminimize": ("OpenBabel", "https://openbabel.org/docs/Installation/install.html"),
    "obenergy": ("OpenBabel", "https://openbabel.org/docs/Installation/install.html"),
}


_BACKEND_MODULE_REQUIREMENT_HINTS = {
    "mlatom": "https://mlatom.com/docs/installation.html",
    "torch": "https://pytorch.org/get-started/locally/",
    "torchani": "https://github.com/aiqm/torchani",
    "ase": "https://wiki.fysik.dtu.dk/ase/install.html",
    "geometric": "https://geometric.readthedocs.io/en/1.1/install.html",
}


def _workflow_requirement_messages(run_mode, software, geometry_optimizer):
    """Return user-facing messages for missing requirements on the active path."""
    missing = []
    capabilities = get_backend_capabilities(software)

    if geometry_optimizer == "geometric" and software is not None:
        if importlib.util.find_spec("ase") is None:
            missing.append(
                _missing_module_message(
                    "ase",
                    _BACKEND_MODULE_REQUIREMENT_HINTS["ase"],
                )
            )
        if importlib.util.find_spec("geometric") is None:
            missing.append(
                _missing_module_message(
                    "geometric",
                    _BACKEND_MODULE_REQUIREMENT_HINTS["geometric"],
                )
            )
        if shutil.which("geometric-optimize") is None:
            missing.append(
                _missing_executable_message(
                    "geometric-optimize",
                    "geomeTRIC",
                    "https://geometric.readthedocs.io/en/1.1/install.html",
                )
            )

    canonical_software = normalize_backend_name(software)
    for executable in capabilities.required_executables:
        if shutil.which(executable) is None:
            friendly_name, install_hint = _BACKEND_EXECUTABLE_REQUIREMENT_HINTS.get(
                executable,
                (canonical_software, 'python -m pip install "pyar-chem[all]"'),
            )
            missing.append(
                _missing_executable_message(
                    executable,
                    friendly_name,
                    install_hint,
                )
            )

    for module_name in capabilities.required_python_modules:
        if importlib.util.find_spec(module_name) is None:
                missing.append(
                    _missing_module_message(
                        module_name,
                        _BACKEND_MODULE_REQUIREMENT_HINTS.get(module_name, "https://pypi.org/"),
                    )
                )

    return missing


def _preflight_cli_requirements(run_mode, software, geometry_optimizer):
    """Fail fast when the selected workflow path is missing requirements."""
    missing = _workflow_requirement_messages(run_mode, software, geometry_optimizer)
    if missing:
        raise SystemExit(
            "Missing requirements for the selected workflow:\n- " + "\n- ".join(missing)
        )


def _merge_run_parameters(args):
    """Overlay parsed CLI arguments onto default run parameters."""
    run_parameters = defaultdict(lambda: None, defualt_parameters.values)
    for key, value in args.items():
        if value is not None and run_parameters[key] != value:
            run_parameters[key] = value
    return run_parameters


def _configure_cli_logging(verbosity):
    """Configure the CLI logger and file handler for this invocation."""
    if verbosity == 0:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(name)-12s %(filename)s %(funcName)s '
                                      '%(lineno)d %(levelname)-8s: %(message)s')
    elif verbosity == 1:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.INFO)
    elif verbosity == 2:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.WARNING)
    elif verbosity == 3:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.ERROR)
    else:
        formatter = logging.Formatter('%(message)s')
        logger.setLevel(logging.CRITICAL)

    file_handler = _get_file_handler()
    file_handler.setFormatter(formatter)
    if file_handler not in logger.handlers:
        logger.addHandler(file_handler)


def _log_startup_context(run_parameters, run_mode, input_files, number_of_input_files):
    """Write startup metadata that is independent of workflow dispatch."""
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
    logger.info("Command line: %s", " ".join(shlex.quote(arg) for arg in sys.argv))
    logger.debug(f'Logging level is {{{logger.level}}}')

    logger.debug(f"{number_of_input_files} input files")
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


def _resolve_cli_input_specs(run_parameters, input_files):
    """Return molecule input specs and formula-derived aggregate sizes."""
    if run_parameters['formula']:
        return _expand_formula_inputs(run_parameters['formula'])
    return input_files, None


def _infer_default_multiplicities(input_molecules, charges):
    """Infer spin multiplicity from electron parity when the CLI omits it."""
    multiplicities = []
    for mol, charge in zip(input_molecules, charges):
        n_electrons = sum(mol.atomic_number) - charge
        multiplicity = 1 if n_electrons % 2 == 0 else 2
        multiplicities.append(multiplicity)
        mol.multiplicity = multiplicity
    return multiplicities


def _load_input_molecules(run_parameters, input_specs):
    """Resolve input specs into molecules with validated charge/multiplicity lists."""
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
        multiplicities = _infer_default_multiplicities(input_molecules, charges)
        logger.info('Multiplicity not provided: inferred defaults from electron parity.')

    return input_molecules, charges, multiplicities, scftypes


def _validate_backend_spin_inputs(input_molecules):
    """Reject charge/multiplicity combinations with impossible electron parity."""
    for mol in input_molecules:
        n_electrons = sum(mol.atomic_number) - mol.charge
        if n_electrons % 2 == 0:
            if mol.multiplicity % 2 != 1:
                sys.exit(f"{n_electrons} (even) electrons and multiplicty {mol.multiplicity} (odd) is not pssible for {mol.name}")
        else:
            if mol.multiplicity % 2 == 1:
                sys.exit(f"{n_electrons} (odd) electrons and multiplicty {mol.multiplicity} (even) is not pssible for {mol.name}")


def _build_qc_parameters(run_parameters, args, run_mode):
    """Build backend-aware QC parameters and return reporting metadata."""
    custom_keywords = run_parameters['custom_keywords']
    provided_qc_options = _provided_qc_options(args)
    if run_mode == "react" and run_parameters["geometry_optimizer"] == "geometric":
        provided_qc_options.discard("gamma")
    backend_family_name = "none"
    staged_optimization = False
    ignored_qc_options = []
    if run_parameters["software"] is not None:
        backend_family_name, ignored_qc_options = _validate_backend_qc_options(
            run_parameters["software"], provided_qc_options
        )
        staged_optimization = backend_supports_staged_optimization(run_parameters["software"])
        if run_parameters["geometry_optimizer"] == "geometric" and not backend_supports_geometry_optimization(run_parameters["software"]):
            sys.exit(
                f"Backend '{run_parameters['software']}' cannot be used with geomeTRIC AFIR "
                "optimisation because it does not expose Cartesian energy and gradients."
            )
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
        'geometry_optimizer': run_parameters['geometry_optimizer'],
        'opt_target': run_parameters['opt_target'],
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
    return (
        quantum_chemistry_parameters,
        backend_family_name,
        ignored_qc_options,
        effective_qc_options,
        staged_optimization,
    )


def _log_qc_context(
    run_parameters,
    quantum_chemistry_parameters,
    backend_family_name,
    ignored_qc_options,
    effective_qc_options,
    staged_optimization,
):
    """Log backend and QC settings after unsupported options are masked."""
    logger.info(f'QM Software:   {quantum_chemistry_parameters["software"]}')
    logger.info(f'Geometry optimizer: {quantum_chemistry_parameters["geometry_optimizer"]}')
    logger.info(f'Optimization target: {quantum_chemistry_parameters["opt_target"]}')
    logger.info(f'Backend family: {backend_family_name}')
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


def _resolve_site_constraint(run_parameters, input_molecules):
    """Return the absolute site constraint used by workflow functions."""
    if run_parameters['site'] is None:
        return None
    site = run_parameters['site']
    return [site[0], input_molecules[0].number_of_atoms + site[1]]


def _log_workflow_plan(run_mode, run_parameters, input_molecules, formula_aggregate_sizes, number_of_orientations):
    """Log the user-facing workflow plan."""
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


def _run_aggregate_workflow(
    aggregate,
    run_parameters,
    input_molecules,
    formula_aggregate_sizes,
    number_of_orientations,
    maximum_number_of_seeds,
    quantum_chemistry_parameters,
    site,
):
    """Validate and dispatch aggregate workflow execution."""
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
    result = aggregate(
        input_molecules,
        size_of_aggregate,
        number_of_orientations,
        quantum_chemistry_parameters,
        maximum_number_of_seeds,
        run_parameters['first_pathway'],
        run_parameters['number_of_pathways'],
        site,
        connectivity_policy=run_parameters["connectivity_policy"],
    )
    _log_workflow_result(result)
    logger.info('Total Time: {}'.format(time.time() - t1_0))
    logger.info("Started at {}\nEnded at {}".format(time_started, datetime.datetime.now()))


def _run_solvation_workflow(
    solvate,
    run_parameters,
    input_molecules,
    number_of_orientations,
    maximum_number_of_seeds,
    quantum_chemistry_parameters,
    site,
):
    """Validate and dispatch solvation workflow execution."""
    number_of_solvent_molecules = run_parameters['solvation_size']
    if number_of_solvent_molecules is None:
        sys.exit('For this please provide the number of solvent\nmolecules to be added. Use the following option\n  -ss <int>')
    if len(input_molecules) == 1:
        sys.exit('Please provide more than two molecules.\nThe last input file will be considered as solvent\nand the other molecules as solutes to which solvent\nmolecules will be added.')
    monomer = input_molecules[-1]
    seeds = input_molecules[:-1]
    t1_0 = time.time()
    time_started = datetime.datetime.now()
    result = solvate(
        seeds,
        monomer,
        number_of_solvent_molecules,
        number_of_orientations,
        quantum_chemistry_parameters,
        maximum_number_of_seeds,
        site,
        connectivity_policy=run_parameters["connectivity_policy"],
    )
    _log_workflow_result(result)
    logger.info('Total Time: {}'.format(time.time() - t1_0))
    logger.info("Started at {}\nEnded at {}".format(time_started, datetime.datetime.now()))


def _run_reaction_workflow(
    reaction_workflow,
    run_parameters,
    input_molecules,
    number_of_orientations,
    quantum_chemistry_parameters,
    site,
):
    """Validate and dispatch reaction workflow execution."""
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
    result = reaction_workflow.react(
        input_molecules[0],
        input_molecules[1],
        minimum_gamma,
        maximum_gamma,
        int(number_of_orientations),
        quantum_chemistry_parameters,
        site,
        proximity_factor,
    )
    _log_workflow_result(result)
    logger.info('Total run time: {}'.format(time.time() - zero_time))
    logger.info(f"Started at {time_started}\nEnded at {datetime.datetime.now()}")


def main():
    if len(sys.argv) > 1 and sys.argv[1] in {"trace", "reaction-trace"}:
        _dispatch_trace_subcommand(sys.argv[2:])
        return
    if len(sys.argv) > 1 and sys.argv[1] in {"conformer", "conformers"}:
        _dispatch_conformer_subcommand(sys.argv[2:])
        return
    args = vars(argument_parse())

    run_parameters = _merge_run_parameters(args)
    run_mode = _active_run_mode(run_parameters)
    input_files = run_parameters['input_files'] or []
    number_of_input_files = len(input_files)
    _validate_cli_request_shape(run_parameters, number_of_input_files)
    _configure_reaction_optimizer(run_parameters, run_mode)
    _preflight_cli_requirements(
        run_mode,
        run_parameters["software"],
        run_parameters["geometry_optimizer"],
    )

    from pyar.state.aggregate import AggregateStateError
    from pyar.state.solvation import SolvationStateError
    from pyar.workflows.aggregate import aggregate
    from pyar.workflows import reaction as reaction_workflow
    from pyar.workflows.solvation import solvate
    from pyar.state.reaction import ReactionStateError

    _configure_cli_logging(run_parameters['verbosity'])
    _log_startup_context(run_parameters, run_mode, input_files, number_of_input_files)
    input_specs, formula_aggregate_sizes = _resolve_cli_input_specs(run_parameters, input_files)
    input_molecules, _, _, _ = _load_input_molecules(run_parameters, input_specs)
    if run_parameters['software'] is not None:
        _validate_backend_spin_inputs(input_molecules)

    (
        quantum_chemistry_parameters,
        backend_family_name,
        ignored_qc_options,
        effective_qc_options,
        staged_optimization,
    ) = _build_qc_parameters(run_parameters, args, run_mode)
    _log_qc_context(
        run_parameters,
        quantum_chemistry_parameters,
        backend_family_name,
        ignored_qc_options,
        effective_qc_options,
        staged_optimization,
    )

    number_of_orientations = run_parameters['how_many_orientations']
    logger.info(f'Number of orientations: {number_of_orientations}')
    maximum_number_of_seeds = run_parameters['maximum_number_of_seeds']
    logger.info(f'Maximum number of seeds: {maximum_number_of_seeds}')

    site = _resolve_site_constraint(run_parameters, input_molecules)
    logger.info(f'Site constraint: {site}')
    _log_workflow_plan(
        run_mode,
        run_parameters,
        input_molecules,
        formula_aggregate_sizes,
        number_of_orientations,
    )
    try:
        if run_parameters['aggregate']:
            _run_aggregate_workflow(
                aggregate,
                run_parameters,
                input_molecules,
                formula_aggregate_sizes,
                number_of_orientations,
                maximum_number_of_seeds,
                quantum_chemistry_parameters,
                site,
            )

        if run_parameters['solvate']:
            _run_solvation_workflow(
                solvate,
                run_parameters,
                input_molecules,
                number_of_orientations,
                maximum_number_of_seeds,
                quantum_chemistry_parameters,
                site,
            )

        if run_parameters['react']:
            _run_reaction_workflow(
                reaction_workflow,
                run_parameters,
                input_molecules,
                number_of_orientations,
                quantum_chemistry_parameters,
                site,
            )
            return

    except (FileNotFoundError, AggregateStateError, ReactionStateError, SolvationStateError, ValueError) as exc:
        logger.critical(str(exc))
        sys.exit(str(exc))


if __name__ == "__main__":
    main()
