"""Describe only explicitly accepted observations; never select recovery actions."""
from dataclasses import dataclass
import math

from .models import StepObservation


@dataclass(frozen=True)
class TrajectorySummary:
    accepted_steps: int
    rejected_steps: int
    unknown_steps: int
    effective_energy_change_hartree: float | None
    gradient_rms_change_hartree_per_bohr: float | None


def analyze_accepted_steps(steps: tuple[StepObservation, ...]):
    accepted = [step for step in steps if step.accepted is True]
    if any(not math.isfinite(step.effective_energy_hartree) for step in accepted):
        raise ValueError('Accepted energy observations must be finite')
    energy_change = accepted[-1].effective_energy_hartree - accepted[0].effective_energy_hartree if len(accepted) > 1 else None
    gradients = [s.gradient_rms_hartree_per_bohr for s in accepted if s.gradient_rms_hartree_per_bohr is not None]
    if any(not math.isfinite(g) for g in gradients):
        raise ValueError('Accepted gradient observations must be finite')
    return TrajectorySummary(len(accepted), sum(s.accepted is False for s in steps),
                             sum(s.accepted is None for s in steps), energy_change,
                             gradients[-1] - gradients[0] if len(gradients) > 1 else None)
