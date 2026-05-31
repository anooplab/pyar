"""Compatibility shim for the legacy ``pyar.reactor`` path.

Use :mod:`pyar.workflows.reaction` instead.
"""

from pyar.workflows.reaction import (
    build_gamma_schedule,
    build_reaction_request,
    format_gamma_id,
    initialize_reaction_run,
    optimize_all,
    print_header,
    react,
    relax_without_afir_bias,
    saved_product_identities,
    with_gamma,
    without_afir_bias,
)

__all__ = [
    "build_gamma_schedule",
    "build_reaction_request",
    "format_gamma_id",
    "initialize_reaction_run",
    "optimize_all",
    "print_header",
    "react",
    "relax_without_afir_bias",
    "saved_product_identities",
    "with_gamma",
    "without_afir_bias",
]
