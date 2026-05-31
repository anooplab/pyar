"""Compatibility shim for the legacy ``pyar.aggregator`` path.

Use :mod:`pyar.workflows.aggregate` and :mod:`pyar.workflows.solvation`
instead.
"""

from pyar.workflows.aggregate import (
    aggregate,
    aggregate_from_formulas,
    expand_formula_to_aggregate_inputs,
    generate_molecule_from_formula,
)
from pyar.workflows.solvation import solvate

__all__ = [
    "aggregate",
    "aggregate_from_formulas",
    "expand_formula_to_aggregate_inputs",
    "generate_molecule_from_formula",
    "solvate",
]
