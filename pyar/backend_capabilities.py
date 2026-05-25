"""Backend capability registry for PyAR workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import FrozenSet


@dataclass(frozen=True)
class BackendCapabilities:
    """Describe what a backend can do for PyAR workflows."""

    family: str = "unknown"
    energy_gradient: bool = False
    native_optimization: bool = True
    supports_charge: bool = True
    supports_multiplicity: bool = True
    staged_optimization: bool = False
    supported_options: FrozenSet[str] = field(default_factory=frozenset)


BACKEND_CAPABILITIES = {
    "gaussian": BackendCapabilities(
        family="dft_qc",
        energy_gradient=True,
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset({"method", "basis", "scf_cycles", "nprocs"}),
    ),
    "orca": BackendCapabilities(
        family="dft_qc",
        energy_gradient=True,
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset({"method", "basis", "scf_cycles", "nprocs"}),
    ),
    "psi4": BackendCapabilities(
        family="dft_qc",
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset(),
    ),
    "turbomole": BackendCapabilities(
        family="dft_qc",
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset({"method", "basis"}),
    ),
    "mopac": BackendCapabilities(
        family="semiempirical",
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset(),
    ),
    "xtb": BackendCapabilities(
        family="semiempirical",
        energy_gradient=True,
        staged_optimization=True,
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset({"opt_threshold", "nprocs"}),
    ),
    "xtb_turbo": BackendCapabilities(
        family="semiempirical",
        staged_optimization=False,
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset({"nprocs"}),
    ),
    "obabel": BackendCapabilities(
        family="semiempirical",
        supports_charge=False,
        supports_multiplicity=False,
        supported_options=frozenset(),
    ),
    "aimnet_2": BackendCapabilities(
        family="mlip",
        energy_gradient=True,
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset(),
    ),
    "aiqm1_mlatom": BackendCapabilities(
        family="mlip",
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset(),
    ),
    "mlatom_aiqm1": BackendCapabilities(
        family="mlip",
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset(),
    ),
    "xtb-aimnet2": BackendCapabilities(
        family="hybrid",
        staged_optimization=True,
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset({"opt_threshold", "nprocs"}),
    ),
    "xtb-aiqm1": BackendCapabilities(
        family="hybrid",
        staged_optimization=True,
        supports_charge=True,
        supports_multiplicity=True,
        supported_options=frozenset({"opt_threshold", "nprocs"}),
    ),
}


def register_backend_capabilities(software, capabilities):
    """Register or replace capability metadata for a backend."""
    if not isinstance(capabilities, BackendCapabilities):
        raise TypeError("capabilities must be a BackendCapabilities instance")
    BACKEND_CAPABILITIES[software] = capabilities


def backend_family(software):
    """Return the backend family label."""
    return get_backend_capabilities(software).family


def backend_supports_staged_optimization(software):
    """Return True when the backend has wired loose/normal staged optimization."""
    return get_backend_capabilities(software).staged_optimization


def get_backend_capabilities(software):
    """Return the capability profile for a backend, if known."""
    return BACKEND_CAPABILITIES.get(software, BackendCapabilities())


def backend_supports_geometry_optimization(software):
    """Return True when geomeTRIC can resolve an energy-gradient provider."""
    from pyar.energy_gradient_providers import ENERGY_GRADIENT_PROVIDERS

    return (
        get_backend_capabilities(software).energy_gradient
        and software in ENERGY_GRADIENT_PROVIDERS
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
