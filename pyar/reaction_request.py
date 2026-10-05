"""Read-only modern reaction resolution, two-phase preflight and restart checks.

The scientific workflow, controller, release gate and persisted request schema
remain owned by their existing modules. This resolver accepts plain mappings.
"""

from dataclasses import dataclass
from importlib import import_module
import math
from pathlib import Path

from pyar.backend_capabilities import (BACKEND_CAPABILITIES, normalize_backend_name,
                                      validate_backend_capability, get_backend_capabilities)
from pyar.biases.controller import resolve_controller_policy
from pyar.biases.softmin import resolve_softmin_beta
from pyar.data.defualt_parameters import values as defaults
from pyar.optimization_request import preflight
from pyar.state.reaction import ReactionRunState, ReactionStateError
from pyar.utility_io import load_structures, expand_charges
from pyar.workflows import reaction as reaction_workflow


@dataclass
class ReactionRequest:
    reactants: list
    qc_params: dict
    bias_min: float | None
    bias_max: float
    orientations: int
    site: list | None
    proximity_factor: float
    restart_request: dict


def resolve_reaction_request(input_files, options):
    """Resolve modern defaults independently of an argparse parser."""
    def option(name, default=None):
        value = options.get(name)
        return default if value is None else value

    if len(input_files) != 2:
        raise ValueError('Reaction search requires exactly two XYZ reactants')
    if not option('backend'):
        raise ValueError('A backend is required. Example: pyar react A.xyz B.xyz --backend xtb --bias-max 100')
    if option('bias_max') is None:
        raise ValueError('Reaction search requires --bias-max, a maximum bias strength.\n'
                         'Example: pyar react A.xyz B.xyz --backend xtb --bias-max 100')
    backend = normalize_backend_name(option('backend').lower())
    if backend not in BACKEND_CAPABILITIES:
        raise ValueError(f'Unknown reaction backend: {backend!r}')
    validate_backend_capability(backend, ('energy_gradient', 'biased_optimization', 'native_optimization'),
                                context='modern reaction search and native release')
    optimizer = option('geometry_optimizer', 'geometric')
    if optimizer != 'geometric':
        raise ValueError('Modern reaction bias requires --geometry-optimizer geometric (geomeTRIC/TRIC)')
    controller = resolve_controller_policy(
        option('bias_controller', 'adaptive'), alpha_min=option('bias_alpha_min'),
        safety_margin=option('bias_alpha_margin'), smoothing=option('bias_alpha_smoothing'),
        epsilon=option('bias_alpha_epsilon'), scheduled_alpha=option('bias_scheduled_alpha'))
    ceiling = option('bias_max')
    if not math.isfinite(ceiling) or ceiling < 0:
        raise ValueError('--bias-max must be finite and nonnegative')
    minimum = option('bias_min')
    if controller == 'adaptive':
        if minimum is not None and (not math.isfinite(minimum) or minimum < 0):
            raise ValueError('--bias-min must be finite and nonnegative')
        minimum = None  # Same canonical adaptive semantics: min is ignored.
        schedule = [float(ceiling)]
    else:
        if minimum is None:
            raise ValueError('Fixed or scheduled reaction bias requires --bias-min')
        schedule = reaction_workflow.build_gamma_schedule(minimum, ceiling)
    orientations = option('orientations', defaults['how_many_orientations'])
    if not isinstance(orientations, int) or orientations < 1:
        raise ValueError('--orientations must be a positive integer')
    proximity = option('proximity_factor', 2.3)
    if not math.isfinite(proximity) or proximity <= 0:
        raise ValueError('--proximity-factor must be finite and positive')
    retries = option('release_retry_limit', 2)
    factor = option('release_margin_factor', 2.0)
    if not isinstance(retries, int) or retries < 0 or not math.isfinite(factor) or factor <= 1:
        raise ValueError('--release-retry-limit must be nonnegative; --release-margin-factor must be finite and > 1')
    fraction = option('release_distance_fraction', .95)
    from pyar.release import ReleaseTracker
    ReleaseTracker(distance_fraction=fraction)
    molecules = load_structures(input_files)
    from pyar.cli import _infer_default_multiplicities, _validate_backend_spin_inputs
    charges = expand_charges(option('charge', [0]), 2)
    multiplicities = None if option('multiplicity') is None else expand_charges(option('multiplicity'), 2)
    scftypes = expand_charges(option('scftype', ['rhf']), 2)
    if any(value not in {'rhf', 'uhf'} for value in scftypes):
        raise ValueError('--scftype must be rhf or uhf')
    if multiplicities is not None and any(value < 1 for value in multiplicities):
        raise ValueError('--multiplicity must be positive')
    for index, molecule in enumerate(molecules):
        molecule.charge = charges[index]
        molecule.scftype = scftypes[index]
        if sum(molecule.atomic_number) - molecule.charge < 1:
            raise ValueError(f'{input_files[index]}: charge leaves no electrons')
    if multiplicities is None:
        multiplicities = _infer_default_multiplicities(molecules, charges)
    for molecule, spin in zip(molecules, multiplicities):
        molecule.multiplicity = spin
    try:
        _validate_backend_spin_inputs(molecules)
        merged = molecules[0].merged_with(molecules[1])
        _validate_backend_spin_inputs([merged])
    except SystemExit as exc:
        raise ValueError(f'Invalid charge/multiplicity: {exc}') from exc
    site = option('site')
    if site is not None:
        if len(site) != 2:
            raise ValueError('--site requires one fragment-local atom index per reactant')
        for index, atom in enumerate(site):
            if not isinstance(atom, int) or not 0 <= atom < molecules[index].number_of_atoms:
                raise ValueError(f'Site index {atom} is out of range for reactant {"AB"[index]} '
                                 f'(0..{molecules[index].number_of_atoms - 1})')
        site = [site[0], molecules[0].number_of_atoms + site[1]]
    settings = {'software': backend, 'geometry_optimizer': optimizer, 'opt_target': 'minimum',
                'bias_controller': controller, 'bias_potential': option('bias_potential', 'afir'),
                'softmin_beta': resolve_softmin_beta(option('softmin_beta')),
                'index': molecules[0].number_of_atoms - 1,
                'release_retry_limit': retries, 'release_margin_factor': factor,
                'release_distance_fraction': fraction, 'trace_enabled': True,
                'method': option('method'), 'basis': option('basis'),
                'nprocs': option('nprocs', 1), 'opt_cycles': option('opt_cycles', defaults['opt_cycles']),
                'opt_threshold': option('opt_threshold', defaults['opt_threshold']),
                'scf_cycles': option('scf_cycles', defaults['scf_cycles']),
                'charge': merged.charge, 'multiplicity': merged.multiplicity,
                'scftype': merged.scftype, 'xtb_unpaired_electrons': merged.multiplicity - 1}
    if settings['bias_potential'] not in {'afir', 'softmin'}:
        raise ValueError('--bias-potential must be afir or softmin')
    for name in ('nprocs', 'opt_cycles', 'scf_cycles'):
        if not isinstance(settings[name], int) or settings[name] < 1:
            raise ValueError(f'--{name.replace("_", "-")} must be a positive integer')
    for name in ('bias_alpha_min', 'bias_alpha_margin', 'bias_alpha_smoothing',
                 'bias_alpha_epsilon', 'bias_scheduled_alpha'):
        if option(name) is not None:
            settings[name] = option(name)
    if backend == 'xtb' and merged.multiplicity == 1 and merged.scftype != 'rhf':
        raise ValueError('Standalone xTB singlet reactions require --scftype rhf')
    if backend == 'xtb':
        from pyar.backends.xtb_utils import canonical_xtb_model
        settings['xtb_model'] = canonical_xtb_model()
    capabilities = get_backend_capabilities(backend)
    if not capabilities.supports_method_basis_options and (settings['method'] is not None or settings['basis'] is not None):
        raise ValueError(f'Backend {backend!r} does not support --method or --basis')
    for name in ('method', 'basis'):
        if settings[name] is not None and not settings[name].strip():
            raise ValueError(f'--{name} must not be empty')
    if capabilities.supports_method_basis_options:
        if not settings['method'] or not settings['method'].strip():
            raise ValueError(f'Backend {backend!r} requires an explicit --method')
        if backend == 'orca':
            from pyar.backends.orca_methods import orca_method, orca_method_keywords
            method, is_xtb = orca_method(settings['method'])
            if method == 'g-xTB':
                raise ValueError('ORCA external g-xTB has no supported Cartesian provider; use standalone xTB')
            if method.lower() in {'r2scan-3c', 'b97-3c', 'pbeh-3c', 'hf-3c'}:
                raise ValueError('ORCA composite methods are not supported by the current reaction gradient route')
            if is_xtb and settings['basis'] is not None:
                raise ValueError('ORCA xTB methods do not accept --basis')
            orca_method_keywords(settings)
        elif not settings['basis'] or not settings['basis'].strip():
            raise ValueError(f'Backend {backend!r} requires an explicit --basis')
    if backend == 'gaussian':
        raise ValueError('The current native Gaussian adapter does not request geometry optimization; '
                         'modern react requires a valid native unbiased release route. Use another backend.')
    if settings['opt_threshold'] not in {'loose', 'normal', 'tight'}:
        raise ValueError('--opt-threshold must be loose, normal, or tight for the shared biased/native route')
    request = reaction_workflow.build_reaction_request(*molecules, schedule, orientations, settings, site, proximity)
    return ReactionRequest(molecules, settings, minimum, ceiling, orientations, site, proximity, request)


def preflight_reaction(request):
    """Validate biased and unbiased requirements without constructing calculators."""
    from pyar.cli import _workflow_requirement_messages
    backend = request.qc_params['software']
    example = f'pyar react A.xyz B.xyz --backend {backend} --bias-max {request.bias_max} --check'
    try:
        requirements = preflight(request.qc_params, request.reactants, check_example=example)
        native = reaction_workflow.without_afir_bias(request.qc_params)
        requirements += preflight(native, [request.reactants[0].merged_with(request.reactants[1])],
                                  check_example=example)
    except (ValueError, ImportError) as exc:
        # Reuse central requirement messages; add the actual PyAR extra first
        # for Python requirements, without suggesting pip for external programs.
        if 'Python' in str(exc) or 'import' in str(exc):
            extra = get_backend_capabilities(backend).optional_extra_hint
            hints = 'pip install "pyar-chem[xtb]"'
            if extra:
                hints += f' and pip install "pyar-chem[{extra}]"'
            raise ValueError(f'Reaction Python requirements: install with {hints}.\n{exc}') from exc
        raise
    missing = _workflow_requirement_messages('react', 'obabel', 'native')
    # Only obabel converts product identities; force-field tools are not needed.
    missing = [message for message in missing if "'obabel'" in message]
    if missing:
        raise ValueError('\n'.join(missing) + '\nThen verify with ' + example)
    import_module('pyar.backends.geometric')
    import_module('pyar.backends.adaptive_geometric')
    import_module('pyar.backends.' + backend)
    return sorted(set(requirements) | {'obabel', 'Cartesian energy-gradient provider', 'native unbiased relaxation'})


def validate_restart(request, root=None):
    """Use the canonical persisted request and validation without writing state."""
    root = Path.cwd() if root is None else Path(root)
    state = ReactionRunState.load(root, request.restart_request)
    if state is not None:
        # Validate persisted geometry references before the workflow can add
        # identity/state metadata or resume its first remaining calculation.
        state.pending_molecules()
        state.current_survivor_molecules()
    if state is None and (root/'reaction').exists():
        raise ReactionStateError('Existing reaction directory has no resumable state; '
                                 'preserve it and start in a new directory.')
    if state is None and (root/'jobs.pkl').exists():
        # All modern routes use physical-v1, so canonical legacy migration rejects them.
        raise ReactionStateError('Legacy reaction checkpoints do not identify physical energies; '
                                 'start a new calculation in a new directory.')
    return state
