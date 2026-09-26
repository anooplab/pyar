"""Optional backend bond-order analysis.

Bond orders are diagnostics only. They are never used to build optimizer
topologies or to decide whether a backend calculation is valid.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import re


@dataclass(frozen=True)
class BondingAnalysis:
    """Backend-neutral result for an optional bonding analysis."""

    scheme: str | None
    bond_orders: dict[tuple[int, int], float] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    available: bool = True
    reason: str | None = None

    def bond_order(self, left, right):
        """Return a symmetric pair value, or ``None`` when it is absent."""
        pair = tuple(sorted((int(left), int(right))))
        return self.bond_orders.get(pair)


def unavailable_bonding_analysis(backend, reason="not_available"):
    """Return an explicit non-fatal result when no parser is available."""
    return BondingAnalysis(
        scheme=None,
        metadata={"backend": str(backend)},
        available=False,
        reason=str(reason),
    )


def _validate_pair(left, right, expected_atoms):
    left, right = int(left), int(right)
    if left == right or left < 0 or right < 0:
        raise ValueError(f"Invalid bond-order atom pair ({left}, {right})")
    if expected_atoms is not None and (left >= expected_atoms or right >= expected_atoms):
        raise ValueError(
            f"Bond-order atom pair ({left}, {right}) exceeds {expected_atoms} atoms"
        )
    return tuple(sorted((left, right)))


_ORCA_PAIR = re.compile(
    r"B\(\s*(\d+)\s*-[^,]+,\s*(\d+)\s*-[^)]+\)\s*:\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
)


def parse_orca_mayer_bond_orders(text, expected_atoms=None):
    """Parse ORCA's printed Mayer bond-order section.

    ORCA prints zero-based atom indices in entries such as
    ``B( 0-C , 1-C ) : 1.3850``. Values below ORCA's print threshold are not
    present in the output and therefore remain absent from the result.
    """
    text = str(text)
    markers = list(re.finditer(r"MAYER\s+POPULATION\s+ANALYSIS", text, re.IGNORECASE))
    if not markers:
        raise ValueError("ORCA Mayer population-analysis section is missing")
    # ORCA may print a population analysis after each optimization step. Only
    # the final complete block describes the final coordinates in stdout.
    section = text[markers[-1].end():]
    matches = list(_ORCA_PAIR.finditer(section))
    if not matches:
        if re.search(r"Mayer\s+bond\s+orders\s+larger\s+than", section, re.IGNORECASE):
            return BondingAnalysis(
                scheme="mayer",
                bond_orders={},
                metadata={"backend": "orca", "source": "Mayer population analysis"},
            )
        raise ValueError("ORCA Mayer bond-order entries are missing or malformed")
    bond_orders = {}
    for match in matches:
        pair = _validate_pair(match.group(1), match.group(2), expected_atoms)
        bond_orders[pair] = float(match.group(3))
    return BondingAnalysis(
        scheme="mayer",
        bond_orders=bond_orders,
        metadata={"backend": "orca", "source": "Mayer population analysis"},
    )


_XTB_ROW = re.compile(
    r"^\s*\d+\s+\S+\s+\S+\s+(.*?)\s*$"
)
_XTB_PAIR = re.compile(
    r"(?P<symbol>[A-Za-z][A-Za-z0-9]*)\s+(?P<index>\d+)\s+"
    r"(?P<value>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
)


def parse_xtb_wiberg_bond_orders(text, expected_atoms=None):
    """Parse xTB's ``Wiberg/Mayer (AO) data`` atom-wise printout.

    xTB atom numbers are one-based; the public result normalizes them to the
    repository's zero-based pair convention. Repeated symmetric entries are
    required to agree within numerical output precision.
    """
    text = str(text)
    marker = re.search(r"Wiberg/Mayer\s+\(AO\)\s+data", text, re.IGNORECASE)
    if marker is None:
        raise ValueError("xTB Wiberg/Mayer bond-order section is missing")
    section = text[marker.end():]
    bond_orders = {}
    row_count = 0
    for line in section.splitlines():
        if line.strip().startswith("molecular dipole"):
            break
        row = _XTB_ROW.match(line)
        if row is None:
            continue
        row_count += 1
        source_match = re.match(r"^\s*(\d+)\s+", line)
        if source_match is None:
            continue
        source = int(source_match.group(1)) - 1
        if expected_atoms is not None and not 0 <= source < expected_atoms:
            raise ValueError(f"xTB atom index {source + 1} exceeds {expected_atoms} atoms")
        for pair_match in _XTB_PAIR.finditer(row.group(1)):
            target = int(pair_match.group("index")) - 1
            pair = _validate_pair(source, target, expected_atoms)
            value = float(pair_match.group("value"))
            previous = bond_orders.get(pair)
            if previous is not None and abs(previous - value) > 1.0e-6:
                raise ValueError(f"Conflicting xTB bond-order values for pair {pair}")
            bond_orders[pair] = value
    if row_count == 0 or not bond_orders:
        raise ValueError("xTB Wiberg bond-order rows are missing or malformed")
    return BondingAnalysis(
        scheme="wiberg_ao",
        bond_orders=bond_orders,
        metadata={"backend": "xtb", "source": "Wiberg/Mayer (AO) data"},
    )


def parse_bonding_analysis(backend, text, expected_atoms=None):
    """Parse a supported backend output or return explicit unavailability."""
    backend = str(backend).lower()
    if backend in {"orca", "orca16"}:
        return parse_orca_mayer_bond_orders(text, expected_atoms)
    if backend in {"xtb", "xtb_turbo", "xtbturbo"}:
        return parse_xtb_wiberg_bond_orders(text, expected_atoms)
    return unavailable_bonding_analysis(backend, "no_stable_parser_in_current_integration")


# P2 audit record. These are deliberately descriptive capabilities, not
# optimizer feature flags. Values reflect the current PyAR integration.
BOND_ORDER_AUDIT = {
    "orca": {
        "available": True,
        "scheme": "Mayer",
        "extra_calculation": False,
        "per_step_feasible": True,
        "parser": "parse_orca_mayer_bond_orders",
        "caveats": "Output threshold omits small values; open-shell values include spin density.",
    },
    "xtb": {
        "available": True,
        "scheme": "Wiberg/Mayer (AO)",
        "extra_calculation": False,
        "per_step_feasible": True,
        "parser": "parse_xtb_wiberg_bond_orders",
        "caveats": "Requires --wbo/property output; values are not interchangeable with Mayer values.",
    },
    "gaussian": {
        "available": False,
        "scheme": "Mulliken/Wiberg output not wired",
        "extra_calculation": True,
        "per_step_feasible": False,
        "parser": None,
        "caveats": "Current PyAR Gaussian provider requests energy/forces only.",
    },
    "aimnet_2": {
        "available": False,
        "scheme": None,
        "extra_calculation": None,
        "per_step_feasible": False,
        "parser": None,
        "caveats": "Current calculator exposes energies/forces, not a bond-order observable.",
    },
}
