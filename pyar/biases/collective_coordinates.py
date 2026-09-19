"""Interfragment collective coordinates and contact diagnostics.

Coordinates accepted by this module are Cartesian coordinates in bohr.  The
returned coordinate is in bohr and its Cartesian derivative has shape
``(natoms, 3)`` in the original atom ordering.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math

from autograd import grad
import autograd.numpy as anp
import numpy as np

from pyar.biases.afir import get_covalent_radius


@dataclass(frozen=True)
class ContactPair:
    """One interfragment atom pair used by a collective coordinate."""

    atom_i: int
    atom_j: int
    symbol_i: str
    symbol_j: str
    distance_bohr: float
    covalent_radii_sum_bohr: float
    gap_bohr: float
    scaled_gap_bohr: float
    is_contact: bool


@dataclass(frozen=True)
class ContactDiagnostics:
    """Reporting-only interfragment contact information."""

    pair_count: int
    contacts: tuple[ContactPair, ...]
    contact_count: int
    minimum_distance_bohr: float
    minimum_gap_bohr: float
    closest_pair: tuple[int, int] | None
    effective_pair_count: float | None = None

    def as_dict(self):
        """Return a JSON-serializable representation."""
        return {
            "pair_count": self.pair_count,
            "contact_count": self.contact_count,
            "minimum_distance_bohr": self.minimum_distance_bohr,
            "minimum_gap_bohr": self.minimum_gap_bohr,
            "closest_pair": list(self.closest_pair) if self.closest_pair else None,
            "effective_pair_count": self.effective_pair_count,
            "contacts": [pair.__dict__ for pair in self.contacts],
        }


def _validate_inputs(fragment_indices, atom_symbols, coordinates, kind, beta, distance_power, contact_factor):
    coordinates = np.asarray(coordinates, dtype=float)
    symbols = list(atom_symbols)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or coordinates.shape[0] != len(symbols):
        raise ValueError("coordinates must have shape (len(atom_symbols), 3)")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("coordinates must contain only finite values")
    if kind not in {"softmin", "afir"}:
        raise ValueError("kind must be 'softmin' or 'afir'")
    if not fragment_indices or len(fragment_indices) != 2:
        raise ValueError("collective coordinate requires exactly two fragments")
    fragments = [tuple(int(index) for index in fragment) for fragment in fragment_indices]
    if not all(fragments):
        raise ValueError("collective coordinate fragments must be non-empty")
    flattened = [index for fragment in fragments for index in fragment]
    if len(set(flattened)) != len(flattened):
        raise ValueError("collective coordinate fragments must not overlap or repeat atoms")
    if any(index < 0 or index >= len(symbols) for index in flattened):
        raise ValueError("collective coordinate fragment index is out of range")
    for name, value, positive in (("beta", beta, True), ("distance_power", distance_power, True), ("contact_factor", contact_factor, False)):
        try:
            numeric = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be a finite number") from exc
        if not math.isfinite(numeric) or (positive and numeric <= 0.0) or (not positive and numeric < 0.0):
            comparator = "positive" if positive else "non-negative"
            raise ValueError(f"{name} must be finite and {comparator}")
    return coordinates, symbols, fragments


def evaluate_contact_coordinate(
    fragment_indices,
    atom_symbols,
    coordinates_bohr,
    *,
    kind="softmin",
    beta=1.0,
    distance_power=6.0,
    contact_factor=1.0,
):
    """Return ``(q, dq_dR, diagnostics)`` for two molecular fragments.

    ``kind='softmin'`` is the radius-adjusted soft minimum of interfragment
    distances. ``kind='afir'`` is the inverse-distance weighted separation
    used by PyAR's isotropic AFIR restraint.
    """
    coordinates, symbols, fragments = _validate_inputs(
        fragment_indices, atom_symbols, coordinates_bohr, kind, beta, distance_power, contact_factor
    )
    left, right = fragments
    radii = np.asarray([get_covalent_radius(symbol) for symbol in symbols], dtype=float)
    pairs = [(i, j) for i, j in product(left, right)]
    distances = np.asarray([np.linalg.norm(coordinates[i] - coordinates[j]) for i, j in pairs])
    if np.any(distances <= 0.0):
        raise ValueError("collective coordinate is undefined for coincident interfragment atoms")
    radii_sums = np.asarray([radii[i] + radii[j] for i, j in pairs])

    left_coordinates = anp.asarray(coordinates[list(left)])
    right_coordinates = anp.asarray(coordinates[list(right)])
    left_radii = anp.asarray(radii[list(left)])
    right_radii = anp.asarray(radii[list(right)])

    def coordinate(left_fragment, right_fragment):
        differences = left_fragment[:, None, :] - right_fragment[None, :, :]
        pair_distances = anp.sqrt(anp.sum(differences ** 2, axis=2)).reshape(-1)
        pair_radii = (left_radii[:, None] + right_radii[None, :]).reshape(-1)
        if kind == "softmin":
            scaled_gaps = -float(beta) * (pair_distances - pair_radii)
            maximum = anp.max(scaled_gaps)
            return -(maximum + anp.log(anp.mean(anp.exp(scaled_gaps - maximum)))) / float(beta)
        weights = (pair_radii / pair_distances) ** float(distance_power)
        return anp.sum(weights * pair_distances) / anp.sum(weights)

    q = float(coordinate(left_coordinates, right_coordinates))
    gradient_left = np.asarray(grad(coordinate, 0)(left_coordinates, right_coordinates), dtype=float)
    gradient_right = np.asarray(grad(coordinate, 1)(left_coordinates, right_coordinates), dtype=float)
    dq_dR = np.zeros_like(coordinates)
    dq_dR[list(left)] = gradient_left
    dq_dR[list(right)] = gradient_right

    gaps = distances - radii_sums
    scaled_gaps = distances - float(contact_factor) * radii_sums
    pair_records = tuple(
        ContactPair(i, j, symbols[i], symbols[j], float(distance), float(radius_sum), float(gap), float(scaled_gap), bool(scaled_gap <= 0.0))
        for (i, j), distance, radius_sum, gap, scaled_gap in zip(pairs, distances, radii_sums, gaps, scaled_gaps)
    )
    closest = int(np.argmin(distances))
    effective_pair_count = None
    if kind == "softmin":
        scaled = -float(beta) * gaps
        weights = np.exp(scaled - np.max(scaled))
        weights /= np.sum(weights)
        effective_pair_count = float(1.0 / np.sum(weights ** 2))
    diagnostics = ContactDiagnostics(
        pair_count=len(pair_records), contacts=pair_records,
        contact_count=sum(pair.is_contact for pair in pair_records),
        minimum_distance_bohr=float(distances[closest]), minimum_gap_bohr=float(np.min(gaps)),
        closest_pair=pairs[closest], effective_pair_count=effective_pair_count,
    )
    return q, dq_dR, diagnostics
