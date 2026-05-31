"""Backend capability registry for PyAR workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from dataclasses import replace
from typing import FrozenSet


@dataclass(frozen=True)
class BackendCapabilities:
    """Describe what a backend can do for PyAR workflows."""

    family: str = "unknown"
    aliases: FrozenSet[str] = field(default_factory=frozenset)
    energy_gradient: bool = False
    native_optimization: bool = True
    supports_biased_optimization: bool = False
    supports_charge: bool = True
    supports_multiplicity: bool = True
    supports_charge_multiplicity: bool = True
    supports_method_basis_options: bool = True
    staged_optimization: bool = False
    required_executables: FrozenSet[str] = field(default_factory=frozenset)
    required_python_modules: FrozenSet[str] = field(default_factory=frozenset)
    optional_extra_hint: str = ""
    notes: str = ""
    supported_options: FrozenSet[str] = field(default_factory=frozenset)


BACKEND_CAPABILITIES = {
    "gaussian": BackendCapabilities(
        family="dft_qc",
        aliases=frozenset({"g16"}),
        energy_gradient=True,
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=True,
        supports_biased_optimization=True,
        required_executables=frozenset({"g16"}),
        supported_options=frozenset({"method", "basis", "scf_cycles", "nprocs"}),
        notes="Provides Cartesian energy/gradient evaluation for geomeTRIC AFIR."
    ),
    "orca": BackendCapabilities(
        family="dft_qc",
        aliases=frozenset({"orca16"}),
        energy_gradient=True,
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=True,
        supports_biased_optimization=True,
        required_executables=frozenset({"orca"}),
        supported_options=frozenset({"method", "basis", "scf_cycles", "nprocs"}),
        notes="Provides Cartesian energy/gradient evaluation for geomeTRIC AFIR."
    ),
    "psi4": BackendCapabilities(
        family="dft_qc",
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=True,
        required_executables=frozenset({"psi4"}),
        supported_options=frozenset(),
        notes="Native geometry optimisation only; no geomeTRIC energy/gradient provider."
    ),
    "turbomole": BackendCapabilities(
        family="dft_qc",
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=True,
        supports_biased_optimization=True,
        required_executables=frozenset({"define"}),
        supported_options=frozenset({"method", "basis"}),
        notes="Uses the internal Turbomole loop with optional AFIR coupling."
    ),
    "mopac": BackendCapabilities(
        family="semiempirical",
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        required_executables=frozenset({"mopac", "obabel"}),
        supported_options=frozenset(),
        notes="MOPAC input is generated through OpenBabel."
    ),
    "xtb": BackendCapabilities(
        family="semiempirical",
        energy_gradient=True,
        supports_biased_optimization=True,
        staged_optimization=True,
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        required_executables=frozenset({"xtb"}),
        supported_options=frozenset({"opt_threshold", "nprocs"}),
        notes="xTB can feed geomeTRIC through the Cartesian energy/gradient provider."
    ),
    "xtb_turbo": BackendCapabilities(
        family="semiempirical",
        staged_optimization=False,
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        supports_biased_optimization=True,
        aliases=frozenset({"xtbturbo"}),
        required_executables=frozenset({"define", "xtb"}),
        supported_options=frozenset({"nprocs"}),
        notes="Internal Turbomole loop using xTB gradients."
    ),
    "obabel": BackendCapabilities(
        family="semiempirical",
        supports_charge=False,
        supports_multiplicity=False,
        supports_charge_multiplicity=False,
        supports_method_basis_options=False,
        required_executables=frozenset({"obabel", "obminimize", "obenergy"}),
        supported_options=frozenset(),
        notes="OpenBabel helper used for force-field minimisation and format conversion."
    ),
    "ani": BackendCapabilities(
        family="mlip",
        aliases=frozenset({"torchani"}),
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        required_python_modules=frozenset({"torch", "torchani", "ase"}),
        optional_extra_hint="ml",
        supported_options=frozenset(),
        notes="TorchANI calculator/helper used outside the main workflow router."
    ),
    "aimnet_2": BackendCapabilities(
        family="mlip",
        energy_gradient=True,
        supports_biased_optimization=True,
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        required_python_modules=frozenset({"torch"}),
        optional_extra_hint="aimnet2",
        supported_options=frozenset(),
        notes="Uses external AIMNet2 model assets and the Cartesian energy/gradient provider."
    ),
    "aiqm1_mlatom": BackendCapabilities(
        family="mlip",
        aliases=frozenset({"mlatom_aiqm1"}),
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        required_python_modules=frozenset({"mlatom"}),
        optional_extra_hint="ml",
        supported_options=frozenset(),
        notes="MLatom-backed AIQM1 geometry optimisation."
    ),
    "xtb-aimnet2": BackendCapabilities(
        family="hybrid",
        aliases=frozenset({"xtb_aimnet2"}),
        energy_gradient=False,
        supports_biased_optimization=False,
        staged_optimization=True,
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        required_executables=frozenset({"xtb"}),
        supported_options=frozenset({"opt_threshold", "nprocs"}),
        notes="Two-stage xTB then AIMNet2 refinement."
    ),
    "xtb-aiqm1": BackendCapabilities(
        family="hybrid",
        aliases=frozenset({"xtb_aiqm1"}),
        energy_gradient=False,
        supports_biased_optimization=False,
        staged_optimization=True,
        supports_charge=True,
        supports_multiplicity=True,
        supports_charge_multiplicity=True,
        supports_method_basis_options=False,
        required_executables=frozenset({"xtb"}),
        supported_options=frozenset({"opt_threshold", "nprocs"}),
        notes="Two-stage xTB then AIQM1 refinement."
    ),
}

BACKEND_ALIASES = {}
for backend_name, capabilities in BACKEND_CAPABILITIES.items():
    for alias in capabilities.aliases:
        BACKEND_ALIASES[alias] = backend_name


def normalize_backend_name(software):
    """Return the canonical backend name for a possibly aliased identifier."""
    return BACKEND_ALIASES.get(software, software)


def _build_missing_feature_message(software, feature, context=None):
    """Return a consistent feature validation error message."""
    prefix = f"Backend {software!r}"
    if context:
        prefix += f" cannot be used with {context}"
    feature_messages = {
        "native_optimization": "does not support native optimisation",
        "energy_gradient": "does not expose Cartesian energy and gradients",
        "biased_optimization": "does not support biased optimisation",
        "staged_optimization": "does not support staged optimisation",
        "charge_multiplicity": "does not support charge/multiplicity inputs",
        "method_basis_options": "does not support method/basis options",
    }
    return f"{prefix} because it {feature_messages[feature]}."


def register_backend_capabilities(software, capabilities):
    """Register or replace capability metadata for a backend."""
    if not isinstance(capabilities, BackendCapabilities):
        raise TypeError("capabilities must be a BackendCapabilities instance")
    canonical_name = normalize_backend_name(software)
    if canonical_name != software and canonical_name in BACKEND_CAPABILITIES:
        BACKEND_CAPABILITIES[canonical_name] = replace(
            capabilities,
            aliases=BACKEND_CAPABILITIES[canonical_name].aliases | capabilities.aliases | frozenset({software}),
        )
        BACKEND_ALIASES[software] = canonical_name
    else:
        BACKEND_CAPABILITIES[software] = capabilities
        canonical_name = software
    for alias in capabilities.aliases:
        BACKEND_ALIASES[alias] = canonical_name


def backend_family(software):
    """Return the backend family label."""
    return get_backend_capabilities(software).family


def backend_supports_native_optimization(software):
    """Return True when the backend performs its own optimisation."""
    return get_backend_capabilities(software).native_optimization


def backend_supports_energy_gradient(software):
    """Return True when the backend can provide Cartesian energies and gradients."""
    return get_backend_capabilities(software).energy_gradient


def backend_supports_biased_optimization(software):
    """Return True when the backend participates in biased optimisation."""
    return get_backend_capabilities(software).supports_biased_optimization


def backend_supports_charge_multiplicity(software):
    """Return True when the backend accepts charge and multiplicity inputs."""
    return get_backend_capabilities(software).supports_charge_multiplicity


def backend_supports_method_basis_options(software):
    """Return True when the backend accepts method and basis options."""
    return get_backend_capabilities(software).supports_method_basis_options


def backend_supports_staged_optimization(software):
    """Return True when the backend has wired loose/normal staged optimization."""
    return get_backend_capabilities(software).staged_optimization


def get_backend_capabilities(software):
    """Return the capability profile for a backend, if known."""
    return BACKEND_CAPABILITIES.get(normalize_backend_name(software), BackendCapabilities())


def validate_backend_capability(software, required_features=(), context=None):
    """Validate that a backend advertises the requested capability flags."""
    canonical_name = normalize_backend_name(software)
    capabilities = get_backend_capabilities(canonical_name)
    if not required_features:
        return canonical_name, capabilities

    missing_features = []
    for feature in required_features:
        if feature == "native_optimization" and not capabilities.native_optimization:
            missing_features.append(feature)
        elif feature == "energy_gradient" and not backend_supports_geometry_optimization(canonical_name):
            missing_features.append(feature)
        elif feature == "biased_optimization" and not capabilities.supports_biased_optimization:
            missing_features.append(feature)
        elif feature == "staged_optimization" and not capabilities.staged_optimization:
            missing_features.append(feature)
        elif feature == "charge_multiplicity" and not capabilities.supports_charge_multiplicity:
            missing_features.append(feature)
        elif feature == "method_basis_options" and not capabilities.supports_method_basis_options:
            missing_features.append(feature)
    if missing_features:
        raise ValueError(_build_missing_feature_message(canonical_name, missing_features[0], context=context))
    return canonical_name, capabilities


def backend_supports_geometry_optimization(software):
    """Return True when geomeTRIC can resolve an energy-gradient provider."""
    from pyar.energy_gradient_providers import ENERGY_GRADIENT_PROVIDERS

    return (
        get_backend_capabilities(software).energy_gradient
        and normalize_backend_name(software) in ENERGY_GRADIENT_PROVIDERS
    )


def supported_geometry_backends():
    """Return the list of backends currently allowed on the AFIR/geomeTRIC path."""
    from pyar.energy_gradient_providers import ENERGY_GRADIENT_PROVIDERS

    preferred_order = {"xtb": 0, "aimnet_2": 1, "orca": 2, "gaussian": 3}
    return tuple(
        name
        for name in sorted(
            (
                backend_name
                for backend_name, capabilities in BACKEND_CAPABILITIES.items()
                if capabilities.energy_gradient and backend_name in ENERGY_GRADIENT_PROVIDERS
            ),
            key=lambda backend_name: (
                preferred_order.get(backend_name, 100),
                backend_name,
            ),
        )
    )


def unsupported_qc_options(software, provided_options):
    """Return explicitly provided QC options that are not supported by a backend."""
    profile = get_backend_capabilities(software)
    return sorted(provided_options - set(profile.supported_options))
