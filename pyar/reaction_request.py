"""Read-only modern reaction resolution, two-phase preflight and restart checks.

The scientific workflow, controller, release gate and persisted request schema
remain owned by their existing modules. This resolver accepts plain mappings.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING
from importlib import import_module
import math
from pathlib import Path

from pyar.backend_capabilities import (BACKEND_CAPABILITIES, normalize_backend_name,
                                      validate_backend_capability, get_backend_capabilities)
from pyar.biases.controller import resolve_controller_policy
from pyar.biases.softmin import resolve_softmin_beta
from pyar.data.defualt_parameters import values as defaults
from pyar.optimization_request import preflight
from pyar.state.reaction_v2 import ReactionRunStateV2 as ReactionRunState
from pyar.state.reaction import ReactionStateError
from pyar.utility_io import load_structures, expand_charges
from pyar.workflows import reaction as reaction_workflow


if TYPE_CHECKING:
    from pyar.core.molecule import Molecule


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


def _resolve_reactant_site(site, molecules):
    """Validate fragment-local site indices and map them into the combined system."""
    if site is None:
        return None
    if len(site) != 2:
        raise ValueError('--site requires one fragment-local atom index per reactant')
    for index, atom in enumerate(site):
        if not isinstance(atom, int) or not 0 <= atom < molecules[index].number_of_atoms:
            raise ValueError(f'Site index {atom} is out of range for reactant {"AB"[index]} '
                             f'(0..{molecules[index].number_of_atoms - 1})')
    return [site[0], molecules[0].number_of_atoms + site[1]]


@dataclass(frozen=True, slots=True)
class _BiasSchedule:
    controller: str
    minimum: float | None
    maximum: float
    values: list[float]


@dataclass(frozen=True, slots=True)
class _ReactionControls:
    backend: str
    optimizer: str
    bias: _BiasSchedule
    orientations: int
    proximity: float
    retry_limit: int
    margin_factor: float
    distance_fraction: float


@dataclass(frozen=True, slots=True)
class _ResolvedReactants:
    molecules: list["Molecule"]
    merged: "Molecule"
    site: list[int] | None


def _option(options, name, default=None):
    """Keep the established distinction between None and other falsy values."""
    value = options.get(name)
    return default if value is None else value


def _resolve_reaction_backend(options):
    if not _option(options, "backend"):
        raise ValueError(
            "A backend is required. Example: pyar react A.xyz B.xyz --backend xtb --bias-max 100"
        )
    if _option(options, "bias_max") is None:
        raise ValueError(
            "Reaction search requires --bias-max, a maximum bias strength.\n"
            "Example: pyar react A.xyz B.xyz --backend xtb --bias-max 100"
        )
    backend = normalize_backend_name(_option(options, "backend").lower())
    if backend not in BACKEND_CAPABILITIES:
        raise ValueError(f"Unknown reaction backend: {backend!r}")
    validate_backend_capability(
        backend,
        ("energy_gradient", "biased_optimization", "native_optimization"),
        context="modern reaction search and native release",
    )
    return backend


def _resolve_bias_schedule(options):
    controller = resolve_controller_policy(
        _option(options, "bias_controller", "adaptive"),
        alpha_min=_option(options, "bias_alpha_min"),
        safety_margin=_option(options, "bias_alpha_margin"),
        smoothing=_option(options, "bias_alpha_smoothing"),
        epsilon=_option(options, "bias_alpha_epsilon"),
        scheduled_alpha=_option(options, "bias_scheduled_alpha"),
    )
    ceiling = _option(options, "bias_max")
    if not math.isfinite(ceiling) or ceiling < 0:
        raise ValueError("--bias-max must be finite and nonnegative")
    minimum = _option(options, "bias_min")
    if controller == "adaptive":
        if minimum is not None and (not math.isfinite(minimum) or minimum < 0):
            raise ValueError("--bias-min must be finite and nonnegative")
        minimum = None  # Same canonical adaptive semantics: min is ignored.
        schedule = [float(ceiling)]
    else:
        if minimum is None:
            raise ValueError("Fixed or scheduled reaction bias requires --bias-min")
        schedule = reaction_workflow.build_gamma_schedule(minimum, ceiling)
    return _BiasSchedule(controller, minimum, ceiling, schedule)


def _resolve_release_settings(options):
    retries = _option(options, "release_retry_limit", 2)
    factor = _option(options, "release_margin_factor", 2.0)
    if (
        not isinstance(retries, int)
        or retries < 0
        or not math.isfinite(factor)
        or factor <= 1
    ):
        raise ValueError(
            "--release-retry-limit must be nonnegative; --release-margin-factor must be finite and > 1"
        )
    fraction = _option(options, "release_distance_fraction", 0.95)
    from pyar.release import ReleaseTracker

    ReleaseTracker(distance_fraction=fraction)
    return retries, factor, fraction


def _resolve_reaction_controls(options):
    backend = _resolve_reaction_backend(options)
    optimizer = _option(options, "geometry_optimizer", "geometric")
    if optimizer != "geometric":
        raise ValueError(
            "Modern reaction bias requires --geometry-optimizer geometric (geomeTRIC/TRIC)"
        )
    bias = _resolve_bias_schedule(options)
    orientations = _option(options, "orientations", defaults["how_many_orientations"])
    if not isinstance(orientations, int) or orientations < 1:
        raise ValueError("--orientations must be a positive integer")
    proximity = _option(options, "proximity_factor", 2.3)
    if not math.isfinite(proximity) or proximity <= 0:
        raise ValueError("--proximity-factor must be finite and positive")
    retries, factor, fraction = _resolve_release_settings(options)
    return _ReactionControls(
        backend=backend,
        optimizer=optimizer,
        bias=bias,
        orientations=orientations,
        proximity=proximity,
        retry_limit=retries,
        margin_factor=factor,
        distance_fraction=fraction,
    )


def _assign_reactant_states(molecules, input_files, options, infer_spins):
    """Assign fragment states before the existing fragment/combined validators."""
    charges = expand_charges(_option(options, "charge", [0]), 2)
    multiplicities = (
        None
        if _option(options, "multiplicity") is None
        else expand_charges(_option(options, "multiplicity"), 2)
    )
    scftypes = expand_charges(_option(options, "scftype", ["rhf"]), 2)
    if any(value not in {"rhf", "uhf"} for value in scftypes):
        raise ValueError("--scftype must be rhf or uhf")
    if multiplicities is not None and any(value < 1 for value in multiplicities):
        raise ValueError("--multiplicity must be positive")
    for index, molecule in enumerate(molecules):
        molecule.charge = charges[index]
        molecule.scftype = scftypes[index]
        if sum(molecule.atomic_number) - molecule.charge < 1:
            raise ValueError(f"{input_files[index]}: charge leaves no electrons")
    if multiplicities is None:
        multiplicities = infer_spins(molecules, charges)
    for molecule, spin in zip(molecules, multiplicities):
        molecule.multiplicity = spin


def _resolve_reactants(input_files, options):
    molecules = load_structures(input_files)
    from pyar.cli_validation import (
        _infer_default_multiplicities,
        _validate_backend_spin_inputs,
    )

    _assign_reactant_states(
        molecules, input_files, options, _infer_default_multiplicities
    )
    try:
        _validate_backend_spin_inputs(molecules)
        merged = molecules[0].merged_with(molecules[1])
        _validate_backend_spin_inputs([merged])
    except SystemExit as exc:
        raise ValueError(f"Invalid charge/multiplicity: {exc}") from exc
    site = _resolve_reactant_site(_option(options, "site"), molecules)
    return _ResolvedReactants(molecules, merged, site)


def _build_reaction_settings(options, controls, inputs):
    settings = {
        "software": controls.backend,
        "geometry_optimizer": controls.optimizer,
        "opt_target": "minimum",
        "bias_controller": controls.bias.controller,
        "bias_potential": _option(options, "bias_potential", "afir"),
        "softmin_beta": resolve_softmin_beta(_option(options, "softmin_beta")),
        "index": inputs.molecules[0].number_of_atoms - 1,
        "release_retry_limit": controls.retry_limit,
        "release_margin_factor": controls.margin_factor,
        "release_distance_fraction": controls.distance_fraction,
        "trace_enabled": True,
        "method": _option(options, "method"),
        "basis": _option(options, "basis"),
        "nprocs": _option(options, "nprocs", 1),
        "opt_cycles": _option(options, "opt_cycles", defaults["opt_cycles"]),
        "opt_threshold": _option(options, "opt_threshold", defaults["opt_threshold"]),
        "scf_cycles": _option(options, "scf_cycles", defaults["scf_cycles"]),
        "charge": inputs.merged.charge,
        "multiplicity": inputs.merged.multiplicity,
        "scftype": inputs.merged.scftype,
        "xtb_unpaired_electrons": inputs.merged.multiplicity - 1,
    }
    return settings


def _validate_execution_controls(settings):
    if settings["bias_potential"] not in {"afir", "softmin"}:
        raise ValueError("--bias-potential must be afir or softmin")
    for name in ("nprocs", "opt_cycles", "scf_cycles"):
        if not isinstance(settings[name], int) or settings[name] < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be a positive integer")


def _copy_controller_options(settings, options):
    for name in (
        "bias_alpha_min",
        "bias_alpha_margin",
        "bias_alpha_smoothing",
        "bias_alpha_epsilon",
        "bias_scheduled_alpha",
    ):
        if _option(options, name) is not None:
            settings[name] = _option(options, name)


def _resolve_backend_model(settings, options, merged):
    backend = settings["software"]
    if backend == "xtb" and merged.multiplicity == 1 and merged.scftype != "rhf":
        raise ValueError("Standalone xTB singlet reactions require --scftype rhf")
    if backend == "xtb":
        from pyar.backends.xtb_utils import canonical_xtb_model

        settings["xtb_model"] = canonical_xtb_model(_option(options, "xtb_model"))
    if backend != "xtb" and _option(options, "xtb_model") is not None:
        raise ValueError("--xtb-model requires --backend xtb")


def _validate_orca_reaction_method(settings):
    from pyar.backends.orca_methods import orca_method, orca_method_keywords

    method, is_xtb = orca_method(settings["method"])
    if method == "g-xTB":
        raise ValueError(
            "ORCA external g-xTB has no supported Cartesian provider; use standalone xTB"
        )
    if method.lower() in {"r2scan-3c", "b97-3c", "pbeh-3c", "hf-3c"}:
        raise ValueError(
            "ORCA composite methods are not supported by the current reaction gradient route"
        )
    if is_xtb and settings["basis"] is not None:
        raise ValueError("ORCA xTB methods do not accept --basis")
    orca_method_keywords(settings)


def _validate_method_basis(settings):
    backend = settings["software"]
    capabilities = get_backend_capabilities(backend)
    if not capabilities.supports_method_basis_options and (
        settings["method"] is not None or settings["basis"] is not None
    ):
        raise ValueError(f"Backend {backend!r} does not support --method or --basis")
    for name in ("method", "basis"):
        if settings[name] is not None and not settings[name].strip():
            raise ValueError(f"--{name} must not be empty")
    if capabilities.supports_method_basis_options:
        if not settings["method"] or not settings["method"].strip():
            raise ValueError(f"Backend {backend!r} requires an explicit --method")
        if backend == "orca":
            _validate_orca_reaction_method(settings)
        elif not settings["basis"] or not settings["basis"].strip():
            raise ValueError(f"Backend {backend!r} requires an explicit --basis")


def _validate_native_release_settings(settings):
    backend = settings["software"]
    if backend == "gaussian":
        raise ValueError(
            "The native Gaussian optimization route is not yet qualified; "
            "modern react requires a valid native unbiased release route. Use another backend."
        )
    if settings["opt_threshold"] not in {"loose", "normal", "tight"}:
        raise ValueError(
            "--opt-threshold must be loose, normal, or tight for the shared biased/native route"
        )


def _resolve_calculation_settings(options, controls, inputs):
    settings = _build_reaction_settings(options, controls, inputs)
    _validate_execution_controls(settings)
    _copy_controller_options(settings, options)
    _resolve_backend_model(settings, options, inputs.merged)
    _validate_method_basis(settings)
    _validate_native_release_settings(settings)
    return settings


def resolve_reaction_request(input_files, options):
    """Resolve controls, electronic states and backend settings before signing."""
    if len(input_files) != 2:
        raise ValueError("Reaction search requires exactly two XYZ reactants")
    controls = _resolve_reaction_controls(options)
    inputs = _resolve_reactants(input_files, options)
    settings = _resolve_calculation_settings(options, controls, inputs)
    request = reaction_workflow.build_reaction_request(
        *inputs.molecules,
        controls.bias.values,
        controls.orientations,
        settings,
        inputs.site,
        controls.proximity,
    )
    return ReactionRequest(
        inputs.molecules,
        settings,
        controls.bias.minimum,
        controls.bias.maximum,
        controls.orientations,
        inputs.site,
        controls.proximity,
        request,
    )

def preflight_reaction(request):
    """Validate biased and unbiased requirements without constructing calculators."""
    from pyar.cli_validation import _workflow_requirement_messages
    backend = request.qc_params['software']
    from pyar.cli_validation import _preflight_xtb_model
    _preflight_xtb_model(request.qc_params)
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
            hints = 'pip install "pyar-chem[reaction]"'
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
