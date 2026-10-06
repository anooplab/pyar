"""Modern conformer interface with read-only input and refinement preflight."""
from pathlib import Path
import shutil

from pyar.cli_runtime import prepared_run, started_run, finished_run
from pyar.backend_errors import BackendExecutionError
from pyar.core.molecule import Molecule
from pyar.conformer.request import ConformerRequest, ConformerRequestError
from pyar.modern_workflow import load_fragments, optional_backend, resolve_states, print_check
from pyar.optimization_request import preflight
from pyar.scripts.conformer import build_parser as legacy_parser, workflow_options
from pyar.workflows import conformer as workflow


def build_parser(prog=None):
    return legacy_parser(prog or 'pyar conformer', modern=True)


def resolve_request(args):
    """Validate chemistry without embedding conformers or converting XYZ files."""
    Chem, AllChem = workflow._rdkit_modules()
    fmt = workflow._resolve_input_format(args.input, args.input_format)
    if args.input_format == 'auto' and Path(args.input).suffix.lower() in {'.xyz', '.sdf', '.mol', '.sd'}:
        if not Path(args.input).is_file():
            raise ValueError(f'Input file does not exist: {args.input}')
    if fmt == 'xyz':
        if not Path(args.input).is_file():
            raise ValueError(f'XYZ input file {args.input!r} does not exist')
        if shutil.which('obabel') is None:
            raise ValueError('OpenBabel executable obabel is required for XYZ conformer input. '
                             'Install OpenBabel from https://openbabel.org/docs/Installation/install.html '
                             'and ensure obabel is on PATH.')
        mol = workflow._load_xyz_with_openbabel(args.input, None, Chem)
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol, addCoords=True)
        workflow._force_field_for_molecule(mol, AllChem, args.force_field)
        molecules = [Molecule([atom.GetSymbol() for atom in mol.GetAtoms()],
                              mol.GetConformer().GetPositions(), name='input')]
        default_charge = [0]  # XYZ does not reliably encode total charge.
    else:
        mol, fmt = workflow._load_rdkit_molecule(args.input, fmt, Chem, None)
        workflow._force_field_for_molecule(mol, AllChem, args.force_field)
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        # Parity validation needs atomic numbers, not generated coordinates.
        molecules = [Molecule(atoms, [[0., 0., 0.]] * len(atoms), name='input')]
        default_charge = [workflow._formal_charge(mol)]
    molecule = resolve_states(molecules, None if args.charge is None else [args.charge],
                              None if args.multiplicity is None else [args.multiplicity],
                              None if args.scftype is None else [args.scftype],
                              default_charges=default_charge)[0]
    args.charge, args.multiplicity, args.scftype = molecule.charge, molecule.multiplicity, molecule.scftype
    qc = optional_backend(args)
    if qc:
        # Conformer refinement calls optimise directly, rather than bulk_optimize.
        qc.update(charge=molecule.charge, multiplicity=molecule.multiplicity,
                  scftype=molecule.scftype, xtb_unpaired_electrons=molecule.multiplicity - 1)
    options = workflow_options(args, qc or None)
    options['input_format'] = fmt
    if args.max_iterations < 1 or args.num_threads < 0:
        raise ValueError('--max-iterations must be positive and --num-threads nonnegative')
    request = ConformerRequest.from_options(args.input, **options)
    # The existing conformer workflow does not resume a prior state.
    directory = Path('conformers')
    if directory.exists() and (not directory.is_dir() or (directory / 'state.json').exists()):
        raise ValueError('Existing conformers state cannot be resumed; start in a new directory.')
    return request, options, molecule


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        request, options, molecule = resolve_request(args)
        qc = dict(request.backend_parameters)
        requirements = ['RDKit'] + (['obabel'] if request.input_format == 'xyz' else [])
        if qc:
            requirements += preflight(qc, [molecule], check_example='pyar conformer INPUT --backend BACKEND --check')
        prepared_run(args, request=request.to_state_dict(), backend=qc,
                     molecules=[molecule], requirements=requirements, outputs=['conformers'], state_lists=False)
        if args.check:
            print_check('conformer', qc, requirements, [('Input', args.input), ('Format', request.input_format),
                        ('Conformers', request.num_conformers), ('Seeds', request.num_seeds),
                        ('Top N', request.top_n), ('Force field', request.force_field),
                        ('Torsion kicks', request.torsion_kicks), ('Charge', request.charge),
                        ('Multiplicity', request.multiplicity)])
            return
        started_run()
        result = workflow.conformer_search(args.input, **options)
    except (ValueError, OSError, ImportError, ConformerRequestError,
            workflow.ConformerWorkflowError, BackendExecutionError) as exc:
        parser.error(str(exc))
    finished_run(result)
    print(f'Conformer search {result.status}.\nSelected conformers: {len(result.selected_paths)}\nRun directory: {result.run_directory}')
    if result.status.startswith('failed'):
        raise SystemExit(1)
