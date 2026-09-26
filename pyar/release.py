"""Conservative, geometry-only release evidence for adaptive reactions."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from pyar.reaction_trace import infer_bonds


RELEASE_STATES = {
    "DRIVING",
    "CANDIDATE",
    "FREE_RELAX_PROBE",
    "PRODUCT_CONFIRMED",
    "RELEASE_SURVIVED",
    "RELEASE_FAILED",
    "TERMINATED",
}


class BondOrderTrajectory:
    """Track optional per-pair bond-order trends at accepted geometries."""

    def __init__(self, stabilization_window=2, stabilization_tolerance=0.05,
                 growth_threshold=0.01):
        self.stabilization_window = int(stabilization_window)
        self.stabilization_tolerance = float(stabilization_tolerance)
        self.growth_threshold = float(growth_threshold)
        self.scheme = None
        self.metadata = {}
        self.available = False
        self.history = {}
        self.last_diagnostics = {}

    @staticmethod
    def _analysis_values(analysis):
        if analysis is None:
            return False, None, {}, {}
        if isinstance(analysis, dict):
            available = bool(analysis.get("available", True))
            return available, analysis.get("scheme"), dict(analysis.get("bond_orders", {})), dict(analysis.get("metadata", {}))
        return (
            bool(getattr(analysis, "available", False)),
            getattr(analysis, "scheme", None),
            dict(getattr(analysis, "bond_orders", {}) or {}),
            dict(getattr(analysis, "metadata", {}) or {}),
        )

    def observe(self, symbols, current_bonds, forming_pairs, analysis=None, fragment_indices=None):
        """Record bond orders for all relevant pairs at one accepted geometry."""
        available, scheme, raw_values, metadata = self._analysis_values(analysis)
        if available and self.scheme is not None and scheme != self.scheme:
            # Bond-order schemes have different scales and cannot share a
            # trajectory baseline.
            self.history.clear()
        self.available = available
        self.metadata = metadata
        if not available:
            self.last_diagnostics = {"available": False, "reason": "unavailable"}
            return self.last_diagnostics
        self.scheme = scheme
        values = {
            _parse_pair_key(pair): float(value)
            for pair, value in raw_values.items()
        }
        current_bonds = {tuple(sorted(pair)) for pair in current_bonds}
        forming_pairs = {tuple(sorted(pair)) for pair in forming_pairs}
        # Capture a pre-contact baseline for every cross-fragment pair.
        # Omitted printed values are unknown, not measured zeros.
        fragment_membership = {}
        for fragment, group in enumerate(fragment_indices or ()):
            for atom in group:
                fragment_membership[atom] = fragment
        relevant_pairs = {
            tuple(sorted((left, right)))
            for left in range(len(symbols))
            for right in range(left + 1, len(symbols))
            if fragment_membership
            and fragment_membership.get(left) != fragment_membership.get(right)
        } | current_bonds | forming_pairs
        for pair in relevant_pairs:
            value = values.get(pair)
            self.history.setdefault(pair, []).append(None if value is None else float(value))
        deltas = {}
        stabilized = {}
        for pair in forming_pairs:
            trajectory = self.history[pair]
            complete_trajectory = all(value is not None for value in trajectory)
            deltas[pair] = (
                trajectory[-1] - trajectory[0]
                if complete_trajectory and len(trajectory) >= 2 else None
            )
            recent = trajectory[-self.stabilization_window:]
            stabilized[pair] = (
                len(recent) >= self.stabilization_window
                and all(value is not None for value in recent)
                and max(value for value in recent if value is not None)
                - min(value for value in recent if value is not None)
                <= self.stabilization_tolerance
            )
        comparable = {}
        for pair in forming_pairs:
            element_pair = tuple(sorted((str(symbols[pair[0]]), str(symbols[pair[1]]))))
            comparable[pair] = []
            for other in current_bonds - forming_pairs:
                other_element_pair = tuple(sorted((str(symbols[other[0]]), str(symbols[other[1]]))))
                if other_element_pair == element_pair and other in values:
                    comparable[pair].append(float(values[other]))
        self.last_diagnostics = {
            "available": True,
            "scheme": scheme,
            "metadata": metadata,
            "bond_orders": {str(pair): values[pair] for pair in values if pair in current_bonds | forming_pairs},
            "delta_from_initial": {str(pair): value for pair, value in deltas.items()},
            "stabilized": {str(pair): value for pair, value in stabilized.items()},
            "comparable_bond_orders": {str(pair): values for pair, values in comparable.items()},
            "geometry_fallback_pairs": [str(pair) for pair in sorted(forming_pairs)
                                        if deltas[pair] is None],
            # Optional, incomplete electronic evidence must not permanently
            # veto a geometry-based free-relaxation probe. Complete evidence
            # still has to show growth and stabilization for its own pair.
            "permits_geometry_probe": bool(forming_pairs) and all(
                deltas[pair] is None or (
                    deltas[pair] > self.growth_threshold and stabilized[pair]
                ) for pair in forming_pairs
            ),
            "qualifies": bool(forming_pairs) and all(
                deltas[pair] is not None
                and deltas[pair] > self.growth_threshold
                and stabilized[pair]
                for pair in forming_pairs
            ),
        }
        return self.last_diagnostics

    def state_dict(self):
        return {
            "schema_version": 1,
            "stabilization_window": self.stabilization_window,
            "stabilization_tolerance": self.stabilization_tolerance,
            "growth_threshold": self.growth_threshold,
            "scheme": self.scheme,
            "metadata": self.metadata,
            "available": self.available,
            "history": {str(pair): values for pair, values in self.history.items()},
            "last_diagnostics": self.last_diagnostics,
        }

    def load_state_dict(self, state):
        if (
            state.get("schema_version") != 1
            or int(state.get("stabilization_window")) != self.stabilization_window
            or float(state.get("stabilization_tolerance")) != self.stabilization_tolerance
            or float(state.get("growth_threshold", 0.01)) != self.growth_threshold
        ):
            raise ValueError("Incompatible bond-order trajectory checkpoint configuration")
        self.scheme = state.get("scheme")
        self.metadata = dict(state.get("metadata") or {})
        self.available = bool(state.get("available"))
        self.history = {
            tuple(int(value) for value in key.strip("()").split(",")): list(values)
            for key, values in state.get("history", {}).items()
        }
        self.last_diagnostics = dict(state.get("last_diagnostics") or {})


def _pair_distance(pair, coordinates):
    left, right = pair
    return float(np.linalg.norm(np.asarray(coordinates[left]) - np.asarray(coordinates[right])))


def _inside_covalent_radius_sum(pair, symbols, coordinates, fraction=1.0):
    left, right = pair
    radius_sum = _covalent_radius(symbols[left]) + _covalent_radius(symbols[right])
    return _pair_distance(pair, coordinates) < float(fraction) * radius_sum


def _cross_fragment_pairs(bonds, fragment_indices):
    if not fragment_indices or len(fragment_indices) < 2:
        return set()
    fragments = {
        atom: fragment
        for fragment, group in enumerate(fragment_indices)
        for atom in group
    }
    return {
        tuple(pair) for pair in bonds
        if fragments.get(pair[0]) != fragments.get(pair[1])
    }


@dataclass
class ReleaseEvidence:
    """Serializable evidence and state for one accepted geometry."""

    accepted_step: int = 0
    segment_index: int = -1
    state: str = "DRIVING"
    reason: str = ""
    forming_pairs: tuple = ()
    forming_pair_distances: tuple = ()
    normalized_distances: tuple = ()
    persistence_counter: int = 0
    distance_delta: float | None = None
    alpha: float | None = None
    alpha_critical: float | None = None
    alpha_target: float | None = None
    bond_order_available: bool = False
    bond_order_scheme: str | None = None
    bond_order_delta: dict = None
    bond_order_stabilized: dict = None
    comparable_bond_orders: dict = None

    def as_dict(self):
        result = asdict(self)
        result["forming_pairs"] = [list(pair) for pair in self.forming_pairs]
        result["forming_pair_distances"] = list(self.forming_pair_distances)
        result["normalized_distances"] = list(self.normalized_distances)
        result["bond_order_delta"] = self.bond_order_delta or {}
        result["bond_order_stabilized"] = self.bond_order_stabilized or {}
        result["comparable_bond_orders"] = self.comparable_bond_orders or {}
        return result


class ReleaseTracker:
    """Require persistent, stabilized chemical contact before a free-relax probe."""

    def __init__(
        self,
        *,
        persistence_required=3,
        stabilization_window=2,
        distance_tolerance=0.05,
        alpha_critical_max=0.05,
        distance_fraction=0.95,
    ):
        self.persistence_required = int(persistence_required)
        self.stabilization_window = int(stabilization_window)
        self.distance_tolerance = float(distance_tolerance)
        self.alpha_critical_max = float(alpha_critical_max)
        self.distance_fraction = float(distance_fraction)
        if self.persistence_required < 2:
            raise ValueError("release persistence must be at least 2 accepted geometries")
        if self.stabilization_window < 2:
            raise ValueError("release stabilization window must be at least 2 accepted geometries")
        if (not np.isfinite(self.distance_tolerance) or self.distance_tolerance <= 0.0
                or not np.isfinite(self.alpha_critical_max) or self.alpha_critical_max < 0.0
                or not np.isfinite(self.distance_fraction) or not 0.0 < self.distance_fraction <= 1.0):
            raise ValueError("release tolerances must be positive and finite")
        self.accepted_step = 0
        self.segment_index = -1
        self.state = "DRIVING"
        self.reason = ""
        self.previous_bonds = set()
        self.pair_distances = {}
        self.pair_persistence = {}
        self.low_resistance_count = 0
        self.bond_order_tracker = BondOrderTrajectory(
            stabilization_window=self.stabilization_window,
        )
        self.evidence = ReleaseEvidence()

    def configuration(self):
        return {
            "persistence_required": self.persistence_required,
            "stabilization_window": self.stabilization_window,
            "distance_tolerance": self.distance_tolerance,
            "alpha_critical_max": self.alpha_critical_max,
            "distance_fraction": self.distance_fraction,
        }

    def state_dict(self):
        return {
            "schema_version": 1,
            "configuration": self.configuration(),
            "accepted_step": self.accepted_step,
            "segment_index": self.segment_index,
            "state": self.state,
            "reason": self.reason,
            "previous_bonds": [list(pair) for pair in sorted(self.previous_bonds)],
            "pair_distances": {str(pair): values for pair, values in self.pair_distances.items()},
            "pair_persistence": {str(pair): value for pair, value in self.pair_persistence.items()},
            "low_resistance_count": self.low_resistance_count,
            "bond_order_tracker": self.bond_order_tracker.state_dict(),
            "evidence": self.evidence.as_dict(),
        }

    def load_state_dict(self, state):
        if state.get("schema_version") != 1 or state.get("configuration") != self.configuration():
            raise ValueError("Incompatible release-state checkpoint configuration")
        if state.get("state") not in RELEASE_STATES:
            raise ValueError("Invalid release state in checkpoint")
        self.accepted_step = int(state["accepted_step"])
        self.segment_index = int(state["segment_index"])
        self.state = state["state"]
        self.reason = str(state.get("reason", ""))
        self.previous_bonds = {tuple(pair) for pair in state.get("previous_bonds", [])}
        self.pair_distances = {
            tuple(int(value) for value in key.strip("()").split(",")): list(values)
            for key, values in state.get("pair_distances", {}).items()
        }
        self.pair_persistence = {
            tuple(int(value) for value in key.strip("()").split(",")): int(value)
            for key, value in state.get("pair_persistence", {}).items()
        }
        self.low_resistance_count = int(state.get("low_resistance_count", 0))
        if state.get("bond_order_tracker"):
            self.bond_order_tracker.load_state_dict(state["bond_order_tracker"])
        evidence = dict(state.get("evidence") or {})
        if evidence:
            for key in ("forming_pairs", "forming_pair_distances", "normalized_distances"):
                evidence[key] = tuple(tuple(pair) if key == "forming_pairs" else pair
                                       for pair in evidence.get(key, ()))
            self.evidence = ReleaseEvidence(**evidence)

    def observe_accepted(
        self,
        symbols,
        coordinates_angstrom,
        fragment_indices,
        *,
        alpha=None,
        alpha_critical=None,
        alpha_target=None,
        segment_index=-1,
        bonding_analysis=None,
    ):
        """Advance evidence exactly once for an accepted geometry."""
        if self.state in {"PRODUCT_CONFIRMED", "RELEASE_SURVIVED", "RELEASE_FAILED", "TERMINATED"}:
            return self.evidence
        coordinates = np.asarray(coordinates_angstrom, dtype=float)
        bonds = infer_bonds(symbols, coordinates, self.previous_bonds)
        forming_pairs = _cross_fragment_pairs(bonds, fragment_indices)
        # The topology heuristic intentionally includes a tolerance above the
        # covalent-radius sum.  For release evidence, require a closer contact:
        # require penetration below the radius sum to avoid probing marginal
        # contacts that are likely to dissociate immediately.
        forming_pairs = {
            pair for pair in forming_pairs
            if _inside_covalent_radius_sum(
                pair, symbols, coordinates, fraction=self.distance_fraction
            )
        }
        distances = {
            pair: _pair_distance(pair, coordinates) for pair in forming_pairs
        }
        for pair in list(self.pair_persistence):
            if pair not in forming_pairs:
                self.pair_persistence[pair] = 0
        for pair in forming_pairs:
            self.pair_persistence[pair] = self.pair_persistence.get(pair, 0) + 1
            history = self.pair_distances.setdefault(pair, [])
            history.append(distances[pair])
            del history[:-self.stabilization_window]

        persistent_pairs = {
            pair for pair in forming_pairs
            if self.pair_persistence.get(pair, 0) >= self.persistence_required
        }
        stabilized_pairs = {
            pair for pair in persistent_pairs
            if len(self.pair_distances.get(pair, [])) >= self.stabilization_window
            and max(self.pair_distances[pair]) - min(self.pair_distances[pair]) <= self.distance_tolerance
        }
        resistance_low = (
            alpha_critical is not None
            and float(alpha_critical) <= self.alpha_critical_max
        )
        self.low_resistance_count = self.low_resistance_count + 1 if resistance_low else 0
        self.accepted_step += 1
        self.segment_index = int(segment_index)
        self.previous_bonds = bonds
        bond_order_diagnostics = self.bond_order_tracker.observe(
            symbols,
            bonds,
            forming_pairs,
            bonding_analysis,
            fragment_indices,
        )
        pair = sorted(forming_pairs)
        selected_distances = tuple(distances[item] for item in pair)
        normalized = tuple(
            distances[item] / (
                _covalent_radius(symbols[item[0]]) + _covalent_radius(symbols[item[1]])
            ) for item in pair
        )
        distance_delta = None
        if pair:
            history = self.pair_distances[pair[0]]
            if len(history) >= 2:
                distance_delta = history[-1] - history[-2]
        bond_order_ready = (
            not bond_order_diagnostics.get("available", False)
            or bool(bond_order_diagnostics.get("permits_geometry_probe"))
        )
        all_pairs_stable = bool(forming_pairs) and forming_pairs <= stabilized_pairs
        if all_pairs_stable and bond_order_ready:
            self.state = "CANDIDATE"
            self.reason = (
                "persistent_stabilized_contact_low_resistance"
                if self.low_resistance_count >= self.persistence_required
                else "persistent_stabilized_contact_probe_high_resistance"
            )
            if bond_order_diagnostics.get("geometry_fallback_pairs"):
                self.reason += "_incomplete_bond_order_geometry_fallback"
            elif bond_order_diagnostics.get("available"):
                self.reason += "_bond_order_supported"
        elif forming_pairs:
            self.state = "DRIVING"
            self.reason = "contact_not_yet_persistent_stabilized_or_bond_order_supported"
        else:
            self.state = "DRIVING"
            self.reason = "no_persistent_interfragment_contact"
        self.evidence = ReleaseEvidence(
            accepted_step=self.accepted_step,
            segment_index=self.segment_index,
            state=self.state,
            reason=self.reason,
            forming_pairs=tuple(pair),
            forming_pair_distances=selected_distances,
            normalized_distances=normalized,
            persistence_counter=max((self.pair_persistence.get(item, 0) for item in pair), default=0),
            distance_delta=distance_delta,
            alpha=None if alpha is None else float(alpha),
            alpha_critical=None if alpha_critical is None else float(alpha_critical),
            alpha_target=None if alpha_target is None else float(alpha_target),
            bond_order_available=bool(bond_order_diagnostics.get("available", False)),
            bond_order_scheme=bond_order_diagnostics.get("scheme"),
            bond_order_delta=bond_order_diagnostics.get("delta_from_initial", {}),
            bond_order_stabilized=bond_order_diagnostics.get("stabilized", {}),
            comparable_bond_orders=bond_order_diagnostics.get("comparable_bond_orders", {}),
        )
        return self.evidence


def _covalent_radius(symbol):
    from pyar.data import new_atomic_data

    return float(new_atomic_data.covalent_radius[str(symbol).capitalize()])


def _parse_pair_key(pair):
    if isinstance(pair, str):
        pair = pair.strip().strip("()").split(",")
    return tuple(sorted((int(pair[0]), int(pair[1]))))
