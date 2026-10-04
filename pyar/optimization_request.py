"""Focused resolution and read-only preflight for modern bulk optimization.

Backend capabilities and installation guidance remain owned by the existing
registry and requirement helpers. This is not a config/profile framework.
"""

from importlib import import_module
from pathlib import Path

import numpy as np

from pyar.backend_capabilities import (
    BACKEND_CAPABILITIES, backend_supports_geometry_optimization,
    get_backend_capabilities, normalize_backend_name, unsupported_qc_options,
)
from pyar.core.molecule import Molecule, parse_xyz


def resolve_settings(args):
    if not args.backend:
        raise ValueError('A backend is required.\nExample: pyar optimize *.xyz --backend xtb')
    backend = normalize_backend_name(args.backend.lower())
    if backend not in BACKEND_CAPABILITIES or backend == 'ani':
        raise ValueError(f'Unsupported optimization backend: {args.backend!r}')
    capabilities = get_backend_capabilities(backend)
    optimizer = args.geometry_optimizer or ('geometric' if backend == 'gaussian' else 'native')
    # The existing native Gaussian adapter writes a single-point route without
    # an Opt keyword. Use its established gradient route, without changing it.
    if backend == 'gaussian' and optimizer == 'native':
        raise ValueError('The Gaussian native adapter does not request geometry optimization; '
                         'use --geometry-optimizer geometric')
    if args.opt_target != 'minimum':
        raise ValueError('Transition-state optimization is not supported by the bulk optimize command')
    if optimizer == 'geometric':
        if not backend_supports_geometry_optimization(backend):
            raise ValueError(f'Backend {backend!r} does not support geomeTRIC energy/gradient optimization')
    elif not capabilities.native_optimization:
        raise ValueError(f'Backend {backend!r} does not support native optimization')
    for name in ('nprocs', 'opt_cycles', 'scf_cycles'):
        value = getattr(args, name)
        if value is not None and value < 1:
            raise ValueError(f'--{name.replace("_", "-")} must be positive')
    for name in ('method', 'basis'):
        value = getattr(args, name)
        if value is not None and not value.strip():
            raise ValueError(f'--{name} must not be empty')
    explicit = {name for name in (
        'method', 'basis', 'nprocs', 'opt_cycles', 'opt_threshold',
        'scf_cycles', 'scf_threshold', 'custom_keywords',
    ) if getattr(args, name) is not None}
    unsupported = set(unsupported_qc_options(backend, explicit))
    # Route-specific controls consumed by the existing adapters but not listed
    # in the historical registry's workflow-level option mask.
    if backend == 'orca' or optimizer == 'geometric':
        unsupported.discard('opt_cycles')
    if backend == 'orca':
        unsupported.discard('opt_threshold')
    if optimizer == 'geometric':
        unsupported.discard('opt_threshold')
    if unsupported:
        raise ValueError(f'Backend {backend!r} with {optimizer} does not support: '
                         + ', '.join('--' + name.replace('_', '-') for name in sorted(unsupported)))
    settings = {
        'software': backend, 'geometry_optimizer': optimizer,
        'opt_target': 'minimum', 'gamma': None, 'nprocs': 1,
        'opt_threshold': 'normal', 'opt_cycles': 100, 'scf_cycles': 1000,
        'scf_threshold': 'normal', 'method': None, 'basis': None,
        'custom_keywords': None, 'custom_keyword': None,
        '_electronic_state_from_molecule': True,
    }
    if backend in {'orca', 'gaussian', 'turbomole', 'xtb_turbo'}:
        settings.update(method='BP86', basis='def2-SVP')
    if backend == 'xtb':
        settings['xtb_model'] = 'gfn2'
    for name in explicit:
        settings[name] = getattr(args, name)
    if backend == 'orca':
        from pyar.backends.orca_methods import orca_method, orca_method_keywords, orca_external_method_block
        method, is_xtb = orca_method(settings['method'])
        builtin = method.lower() in {'r2scan-3c', 'b97-3c', 'pbeh-3c', 'hf-3c'}
        if is_xtb or builtin:
            if args.basis is not None:
                raise ValueError('The selected ORCA method defines its own basis; omit --basis')
            settings['basis'] = None
        if builtin:
            if optimizer == 'geometric':
                raise ValueError('ORCA composite methods currently require --geometry-optimizer native')
            settings['orca_builtin_method'] = True
        orca_method_keywords(settings)
        orca_external_method_block(settings)
    return settings


def load_inputs(args):
    # Reuse exactly the established PyAR parity inference/validation rules.
    from pyar.cli import _infer_default_multiplicities, _validate_backend_spin_inputs
    if args.multiplicity is not None and args.multiplicity < 1:
        raise ValueError('--multiplicity must be positive')
    molecules = []
    names = set()
    for filename in args.input_files:
        try:
            atoms, coordinates, _, title, energy = parse_xyz(filename)
            # Names are backend file stems, not input directory paths.
            name = Path(filename).stem
            if name in names:
                raise ValueError(f'Duplicate job name {name!r}; use distinct XYZ basenames')
            names.add(name)
            mol = Molecule(atoms, coordinates, name=name, title=title, energy=energy,
                           charge=args.charge, multiplicity=args.multiplicity or 1)
            if not atoms or not np.isfinite(coordinates).all():
                raise ValueError('XYZ must contain atoms with finite coordinates')
            if sum(mol.atomic_number) - args.charge < 1:
                raise ValueError('Charge leaves no electrons')
        except (ValueError, KeyError) as exc:
            raise ValueError(f'{filename}: {exc}') from exc
        molecules.append(mol)
    if args.multiplicity is None:
        _infer_default_multiplicities(molecules, [args.charge] * len(molecules))
    try:
        _validate_backend_spin_inputs(molecules)
    except SystemExit as exc:
        raise ValueError(f'Invalid charge/multiplicity: {exc}') from exc
    for mol in molecules:
        mol.scftype = args.scftype or ('rhf' if mol.multiplicity == 1 else 'uhf')
    return molecules


def preflight(settings, molecules, *, require_optimizer_executable=True, check_example=None):
    """Check only the selected route, without constructing a backend or writing files."""
    from pyar.cli import _workflow_requirement_messages
    backend = settings['software']
    optimizer = settings['geometry_optimizer']
    capabilities = get_backend_capabilities(backend)
    if not capabilities.supports_charge_multiplicity and any(
        mol.charge or mol.multiplicity != 1 for mol in molecules
    ):
        raise ValueError(f'Backend {backend!r} does not support charged/open-shell inputs')
    missing = _workflow_requirement_messages('optimize', backend, optimizer)
    modules = set(capabilities.required_python_modules)
    requirements = set(capabilities.required_executables) | modules
    if optimizer == 'geometric':
        modules.update({'ase', 'geometric'})
        requirements.update({'ase', 'geometric', 'geometric-optimize'})
        # Match the backend's existing active-environment executable fallback.
        import sys
        import os
        script = Path(sys.executable).with_name('geometric-optimize')
        if script.is_file() and os.access(script, os.X_OK):
            missing = [message for message in missing if "'geometric-optimize'" not in message]
        if not require_optimizer_executable:
            missing = [message for message in missing if "'geometric-optimize'" not in message]
            requirements.discard('geometric-optimize')
    if backend in {'aiqm1_mlatom', 'xtb-aiqm1'}:
        modules.add('mlatom')
        requirements.add('mlatom')
        missing.extend(_workflow_requirement_messages('optimize', 'gaussian', 'native'))
        requirements.add('g16')
    if backend in {'aimnet_2', 'xtb-aimnet2'}:
        from pyar.backends.aimnet2_assets import validate_aimnet2_runtime_assets
        validate_aimnet2_runtime_assets(include_script=optimizer == 'native')
        modules.update({'torch', 'ase'})
        if optimizer == 'native':
            modules.add('openbabel')
        requirements.add('AIMNet2 model assets')
    if missing:
        example = check_example or f'pyar optimize INPUT.xyz --backend {backend} --check'
        raise ValueError('\n'.join(missing) + '\nThen verify with ' + example)
    for module in sorted(modules):
        try:
            import_module(module)
        except ImportError as exc:
            extra = ('openbabel' if module == 'openbabel' else capabilities.optional_extra_hint
                     or ('xtb' if module in {'ase', 'geometric'} else 'ml' if module == 'mlatom' else 'all'))
            raise ImportError(f'Cannot import {module}: {exc}. Install with '
                              f'python -m pip install "pyar-chem[{extra}]"') from exc
    return sorted(requirements | modules)
