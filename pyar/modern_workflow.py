"""Small shared request helpers for optional minimum-optimization workflows.

Scientific workflow requests, backend capabilities and parity rules remain
owned by their existing modules. These helpers perform no filesystem writes.
"""
from pathlib import Path

import numpy as np

from pyar import optimization_request
from pyar.core.molecule import Molecule, parse_xyz


def add_backend_arguments(parser):
    parser.add_argument('--backend', help='Optional backend refinement; omit for geometry-only generation')
    group = parser.add_argument_group('Backend optimization overrides')
    group.add_argument('--geometry-optimizer', choices=('native', 'geometric'))
    for name in ('method', 'basis', 'model', 'gxtb-wrapper', 'custom-keywords', 'scf-threshold'):
        group.add_argument('--' + name)
    group.add_argument('--xtb-model', choices=('gfn2', 'gxtb'))
    group.add_argument('--opt-threshold', choices=('loose', 'normal', 'tight'))
    for name in ('nprocs', 'opt-cycles', 'scf-cycles'):
        group.add_argument('--' + name, type=int)
    parser.set_defaults(opt_target='minimum')


def optional_backend(args):
    options = ('geometry_optimizer', 'method', 'basis', 'model', 'gxtb_wrapper',
               'custom_keywords', 'xtb_model', 'nprocs', 'opt_cycles', 'scf_cycles',
               'opt_threshold', 'scf_threshold')
    if not args.backend:
        if any(getattr(args, key, None) is not None for key in options):
            raise ValueError('Backend optimization options require --backend')
        return {}
    return optimization_request.resolve_settings(args)


def add_fragment_arguments(parser):
    parser.add_argument('--charge', nargs='+', type=int)
    parser.add_argument('--multiplicity', nargs='+', type=int)
    parser.add_argument('--scftype', nargs='+', choices=('rhf', 'uhf'))


def resolve_states(molecules, charges=None, multiplicities=None, scftypes=None, *, default_charges=None):
    """Apply one/per-fragment values using the established PyAR parity rules."""
    from pyar.cli import _normalize_parameter_list, _infer_default_multiplicities, _validate_backend_spin_inputs
    try:
        charges = (default_charges if charges is None and default_charges is not None else
                   _normalize_parameter_list(charges, 0, len(molecules), 'Charges'))
        spins = None if multiplicities is None else _normalize_parameter_list(multiplicities, 1, len(molecules), 'Multiplicities')
        types = None if scftypes is None else _normalize_parameter_list(scftypes, 'rhf', len(molecules), 'SCF types')
        for molecule, charge in zip(molecules, charges):
            molecule.charge = charge
            if not molecule.atoms_list or not np.isfinite(molecule.coordinates).all():
                raise ValueError(f'{molecule.name}: expected atoms with finite coordinates')
            if sum(molecule.atomic_number) - charge < 1:
                raise ValueError(f'{molecule.name}: charge leaves no electrons')
        if spins is None:
            spins = _infer_default_multiplicities(molecules, charges)
        for index, (molecule, spin) in enumerate(zip(molecules, spins)):
            if spin < 1 or spin - 1 > sum(molecule.atomic_number) - molecule.charge:
                raise ValueError('Invalid multiplicity')
            molecule.multiplicity = spin
            molecule.scftype = types[index] if types is not None else ('rhf' if spin == 1 else 'uhf')
        if types is not None and any(value not in {'rhf', 'uhf'} for value in types):
            raise ValueError('--scftype must be rhf or uhf')
        _validate_backend_spin_inputs(molecules)
    except SystemExit as exc:
        raise ValueError(f'Invalid electronic state: {exc}') from exc
    return molecules


def load_fragments(paths, *, allow_formula=False):
    molecules = []
    for spec in paths:
        if allow_formula and not Path(spec).exists():
            from pyar.workflows.aggregate import generate_molecule_from_formula
            try:
                # Formula packing is geometric initialization; make the modern
                # request stable so read-only restart validation is repeatable.
                molecule = generate_molecule_from_formula(spec, rng=np.random.default_rng(1))
            except ValueError as exc:
                raise ValueError(str(exc)) from exc
        else:
            atoms, coordinates, _, title, energy = parse_xyz(spec)
            try:
                molecule = Molecule(atoms, coordinates, name=Path(spec).stem, title=title, energy=energy)
            except KeyError as exc:
                raise ValueError(f'{spec}: unknown element {exc}') from exc
        molecules.append(molecule)
    return molecules


def add_selection_arguments(parser, *, maximum_seeds):
    parser.add_argument('--orientations', '-N', type=int, default=8)
    parser.add_argument('--maximum-number-of-seeds', type=int, default=maximum_seeds)
    parser.add_argument('--connectivity-policy', choices=('auto', 'off', 'prefer', 'strict'), default='auto')
    parser.add_argument('--selection-feature', default='auto')
    parser.add_argument('--selection-algorithm', default='auto')
    parser.add_argument('--selection-distance', default='euclidean')
    parser.add_argument('--selection-system-type', default='auto')
    parser.add_argument('--site', nargs=2, type=int, metavar=('I', 'J'),
                        help='0-based sites in the original seed and incoming fragment (grow)')
    parser.add_argument('--check', action='store_true', help='Read-only validation; generate no structures')


def print_check(task, settings, requirements, details):
    print(f'Preflight: {task}')
    for label, value in details:
        print(f'{label}: {value}')
    print('Backend: ' + (settings.get('software') or 'none (geometry-only)'))
    print('Requirements: ' + (', '.join(requirements) or 'none'))
    print('Ready to run.\n--check specified; no calculations were performed.')
