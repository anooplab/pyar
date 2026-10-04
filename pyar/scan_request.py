"""Read-only validation/preflight for the modern scan-bond interface."""

import math

import numpy as np

from pyar.core.molecule import Molecule, parse_xyz
from pyar.optimization_request import preflight


def _pair(values, label):
    if len(values) == 1:
        return values * 2
    if len(values) != 2:
        raise ValueError(f'{label} accepts one value or one value per fragment')
    return values


def validate_inputs(args):
    from pyar.cli import _infer_default_multiplicities, _validate_backend_spin_inputs
    charges = _pair(args.charge, 'charge')
    spins = None if args.multiplicity is None else _pair(args.multiplicity, 'multiplicity')
    scftypes = None if args.scftype is None else _pair(args.scftype, 'scftype')
    if spins is not None and any(value < 1 for value in spins):
        raise ValueError('multiplicities must be positive')
    if scftypes is not None and any(value not in {'rhf', 'uhf'} for value in scftypes):
        raise ValueError('scftype must be rhf or uhf')
    molecules = []
    for index, path in enumerate(args.inputs):
        try:
            atoms, xyz, name, title, energy = parse_xyz(path)
            mol = Molecule(atoms, xyz, name=name, title=title, energy=energy, charge=charges[index])
            if not atoms or not np.isfinite(xyz).all():
                raise ValueError('XYZ must contain atoms with finite coordinates')
            if sum(mol.atomic_number) - mol.charge < 1:
                raise ValueError('charge leaves no electrons')
            atom = args.atoms[index]
            if not 0 <= atom < mol.number_of_atoms:
                raise ValueError(f'Atom index {atom} is out of range for fragment {"AB"[index]} '
                                 f'(valid 0..{mol.number_of_atoms - 1})')
            molecules.append(mol)
        except (ValueError, KeyError) as exc:
            raise ValueError(f'{path}: {exc}') from exc
    if spins is None:
        spins = _infer_default_multiplicities(molecules, charges)
    for index, mol in enumerate(molecules):
        mol.multiplicity = spins[index]
        # Preserve the scan's established SCF default. Fragment spins are
        # combined by the workflow; two doublets may form a singlet, so a
        # fragment-derived UHF default must not leak into that merged system.
        mol.scftype = scftypes[index] if scftypes else 'rhf'
    try:
        _validate_backend_spin_inputs(molecules)
        merged = molecules[0].merged_with(molecules[1])
        _validate_backend_spin_inputs([merged])
    except SystemExit as exc:
        raise ValueError(f'Invalid charge/multiplicity: {exc}') from exc
    if args.software == 'xtb' and merged.multiplicity == 1 and merged.scftype != 'rhf':
        raise ValueError('Standalone xTB singlet scans require --scftype rhf; spin is determined by multiplicity')
    args.multiplicity = spins
    args.scftype = [mol.scftype for mol in molecules]
    for name in ('scan_end', 'scan_step'):
        value = getattr(args, name)
        if value is not None and (not math.isfinite(value) or value <= 0):
            raise ValueError(f'--{name.replace("_", "-")} must be finite and positive')
    if args.scan_points is not None and args.scan_points < 2:
        raise ValueError('--scan-points must be at least 2')
    return molecules


def preflight_scan(args, params, molecules):
    """Validate APIs for the exact route, without constructing calculators or stages."""
    native_scan_only = args.software == 'orca' and args.through == 'scan'
    settings = dict(params, geometry_optimizer='native' if native_scan_only else 'geometric')
    requirements = preflight(settings, molecules, require_optimizer_executable=False,
                             check_example=f'pyar scan-bond A.xyz B.xyz --atoms 0 0 --backend {args.software} --check')
    if not native_scan_only:
        # These are the in-process APIs used by bond_scan and scan_path/neb.
        from geometric.optimize import Optimizer  # noqa: F401
        from geometric.internal import DelocalizedInternalCoordinates  # noqa: F401
        from geometric.params import OptParams  # noqa: F401
        from pyar.backends.geometric import PyarGeometricCalculator  # noqa: F401
    if args.through != 'scan':
        from geometric.neb import ElasticBand, OptimizeChain  # noqa: F401
        from geometric.ase_engine import EngineASE  # noqa: F401
        if args.interpolation in {'idpp', 'geodesic'}:
            from ase.mep.neb import idpp_interpolate  # noqa: F401
        if args.interpolation == 'geodesic':
            from pyar.neb import _geodesic_api
            _geodesic_api()
            requirements.append('geodesic-interpolate')
        if args.through not in {'scan', 'neb'} and args.ts_optimizer == 'sella':
            from pyar.neb import _sella_api
            _sella_api()
            requirements.append('sella')
        if args.through not in {'neb', 'ts'}:
            from geometric.normal_modes import calc_cartesian_hessian, frequency_analysis  # noqa: F401
    return sorted(set(requirements))
