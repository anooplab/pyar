"""Policies for selecting an instantaneous reaction-bias force scale."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

DEFAULT_ADAPTIVE_MARGIN = 0.001  # Hartree/Bohr


def resolve_controller_policy(
    policy=None,
    *,
    alpha_min=None,
    safety_margin=None,
    smoothing=None,
    epsilon=None,
    scheduled_alpha=None,
):
    """Resolve CLI controller selection and reject options for another policy."""
    adaptive_options = (alpha_min, safety_margin, smoothing, epsilon)
    has_adaptive_options = any(value is not None for value in adaptive_options)
    has_scheduled_option = scheduled_alpha is not None
    if has_adaptive_options and has_scheduled_option:
        raise ValueError("Adaptive controller settings and scheduled alpha cannot be combined")

    requested = None if policy is None else str(policy).lower()
    inferred = "adaptive" if has_adaptive_options else (
        "scheduled" if has_scheduled_option else "fixed"
    )
    resolved = requested or inferred
    if resolved not in BiasController._POLICIES:
        choices = ", ".join(sorted(BiasController._POLICIES))
        raise ValueError(f"Unsupported bias-controller policy: {resolved!r}; choose one of {choices}")
    if has_adaptive_options and resolved != "adaptive":
        raise ValueError("--bias-alpha-* settings require --bias-controller adaptive")
    if has_scheduled_option and resolved != "scheduled":
        raise ValueError("--bias-scheduled-alpha requires --bias-controller scheduled")
    if resolved == "adaptive":
        BiasController(
            "adaptive",
            alpha_min=0.0 if alpha_min is None else alpha_min,
            safety_margin=safety_margin,
            smoothing=1.0 if smoothing is None else smoothing,
            epsilon=1.0e-12 if epsilon is None else epsilon,
        )
    elif resolved == "scheduled":
        BiasController("scheduled", scheduled_alpha=scheduled_alpha)
    return resolved


@dataclass(frozen=True)
class BiasControllerDecision:
    """The bounded force scale and diagnostics chosen for one evaluation."""

    policy: str
    alpha: float
    alpha_critical: float | None
    alpha_target: float


class BiasController:
    """Choose fixed, scheduled, or resistance-based bias strengths.

    The adaptive policy uses ``max(0, -(grad_E . grad_q)/(||grad_q||^2 +
    epsilon))`` plus a positive driving margin. Smoothing damps decreases;
    increases are applied immediately so filtering cannot cancel the drive.
    """

    _POLICIES = {"fixed", "scheduled", "adaptive"}

    def __init__(
        self,
        policy="fixed",
        *,
        alpha_min=0.0,
        safety_margin=None,
        smoothing=1.0,
        epsilon=1.0e-12,
        scheduled_alpha=None,
    ):
        self.policy = str(policy).lower()
        if self.policy not in self._POLICIES:
            raise ValueError(f"Unsupported bias-controller policy: {policy!r}")
        if safety_margin is None:
            safety_margin = DEFAULT_ADAPTIVE_MARGIN if self.policy == "adaptive" else 0.0
        self.alpha_min = self._finite_nonnegative(alpha_min, "alpha_min")
        self.safety_margin = self._finite_nonnegative(safety_margin, "safety_margin")
        if self.policy == "adaptive" and self.safety_margin == 0.0:
            raise ValueError("adaptive safety_margin (--bias-alpha-margin) must be positive")
        self.smoothing = self._finite_nonnegative(smoothing, "smoothing")
        if not 0.0 < self.smoothing <= 1.0:
            raise ValueError("smoothing must be greater than 0 and at most 1")
        self.epsilon = self._finite_positive(epsilon, "epsilon")
        self.scheduled_alpha = None if scheduled_alpha is None else self._finite_nonnegative(
            scheduled_alpha, "scheduled_alpha"
        )
        self._previous_alpha = None
        self.decision = None
        self.energy_offset = 0.0
        self.segment_index = -1
        self.alpha_max = None

    def configuration(self):
        """JSON-safe parameters required to reproduce controller decisions."""
        parameters = {name: getattr(self, name) for name in (
            "policy", "alpha_min", "safety_margin", "smoothing", "epsilon", "scheduled_alpha"
        )}
        if self.policy == "adaptive":
            parameters["smoothing_mode"] = "decrease_only"
        return parameters

    def state_dict(self):
        """Serialize the active, accepted segment (never a trial proposal)."""
        return {
            "schema_version": 2,
            "configuration": self.configuration(),
            "decision": None if self.decision is None else asdict(self.decision),
            "energy_offset_hartree": self.energy_offset,
            "segment_index": self.segment_index,
            "alpha_max": self.alpha_max,
        }

    def load_state_dict(self, state, alpha_max):
        """Restore history only when the checkpoint matches this configuration."""
        if state.get("schema_version") not in {1, 2} or state.get("configuration") != self.configuration():
            raise ValueError("Incompatible bias-controller checkpoint configuration")
        restored_alpha_max = self._finite_nonnegative(
            state.get("alpha_max", alpha_max), "alpha_max"
        )
        configured_alpha_max = self._finite_nonnegative(alpha_max, "alpha_max")
        if not math.isclose(restored_alpha_max, configured_alpha_max, rel_tol=0.0, abs_tol=1.0e-14):
            raise ValueError("Checkpoint alpha_max does not match configured bias strength")
        decision = BiasControllerDecision(**state["decision"])
        for name in ("alpha", "alpha_target"):
            value = self._finite_nonnegative(getattr(decision, name), name)
            if not self.alpha_min <= value <= alpha_max:
                raise ValueError("Checkpoint bias strength is outside configured bounds")
        if decision.policy != self.policy:
            raise ValueError("Checkpoint bias policy does not match")
        if decision.alpha_critical is not None:
            self._finite_nonnegative(decision.alpha_critical, "alpha_critical")
        offset = float(state["energy_offset_hartree"])
        index = state["segment_index"]
        if not math.isfinite(offset) or type(index) is not int or index < 0:
            raise ValueError("Invalid bias-controller checkpoint offset or segment index")
        self.decision = decision
        self._previous_alpha = decision.alpha
        self.energy_offset = offset
        self.segment_index = index
        self.alpha_max = restored_alpha_max

    def start_segment(self, physical_gradient, coordinate_gradient, q, alpha_max):
        """Commit an update at an accepted geometry, preserving its energy.

        Trial evaluations must use ``decision`` without calling this method.
        The offset makes old_alpha*q + old_offset equal new_alpha*q + new_offset.
        """
        q = float(q)
        if not math.isfinite(q):
            raise ValueError("The accepted contact coordinate must be finite")
        decision = self.select(physical_gradient, coordinate_gradient, alpha_max)
        if self.decision is not None:
            self.energy_offset += (self.decision.alpha - decision.alpha) * q
        self.decision = decision
        self._previous_alpha = decision.alpha
        self.segment_index += 1
        self.alpha_max = self._finite_nonnegative(alpha_max, "alpha_max")
        return decision

    @staticmethod
    def _finite_nonnegative(value, name):
        value = float(value)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be a finite non-negative number")
        return value

    @staticmethod
    def _finite_positive(value, name):
        value = float(value)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a finite positive number")
        return value

    def select(self, physical_gradient, coordinate_gradient, alpha_max):
        """Propose a bounded decision without advancing accepted-step history."""
        alpha_max = self._finite_nonnegative(alpha_max, "alpha_max")
        if self.alpha_min > alpha_max:
            raise ValueError("alpha_min must not exceed alpha_max")
        if self.policy == "fixed":
            return BiasControllerDecision("fixed", alpha_max, None, alpha_max)
        if self.policy == "scheduled":
            target = alpha_max if self.scheduled_alpha is None else self.scheduled_alpha
            target = float(np.clip(target, self.alpha_min, alpha_max))
            return BiasControllerDecision("scheduled", target, None, target)

        gradient = np.asarray(physical_gradient, dtype=float)
        coordinate = np.asarray(coordinate_gradient, dtype=float)
        if gradient.shape != coordinate.shape or not np.all(np.isfinite(gradient)) or not np.all(np.isfinite(coordinate)):
            raise ValueError("physical and coordinate gradients must be finite arrays with identical shapes")
        denominator = float(np.sum(coordinate * coordinate) + self.epsilon)
        alpha_critical = max(0.0, -float(np.sum(gradient * coordinate)) / denominator)
        target = float(np.clip(alpha_critical + self.safety_margin, self.alpha_min, alpha_max))
        alpha = target if self._previous_alpha is None else (
            self.smoothing * target + (1.0 - self.smoothing) * self._previous_alpha
        )
        # A lagging increase can leave alpha below the physical resistance.
        # Preserve at least the current driving target, subject to the ceiling.
        alpha = max(alpha, target)
        alpha = float(np.clip(alpha, self.alpha_min, alpha_max))
        return BiasControllerDecision("adaptive", alpha, alpha_critical, target)


class FixedBiasController(BiasController):
    """Explicit fixed-force controller; retained ``BiasController`` API is compatible."""

    def __init__(self):
        super().__init__("fixed")
