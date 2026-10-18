"""Small fail-closed SAFE policy for standalone unbiased closed-shell GFN2."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class RecoveryPolicy:
    mode: str = 'off'
    max_attempts: int = 3
    max_cycles: int = 700
    max_evaluations: int = 1000
    wall_seconds: float = 3600.

    def __post_init__(self):
        if self.mode not in {'off', 'safe', 'guided'}:
            raise ValueError('Recovery mode must be off, safe or guided')
        for key, maximum in (('max_attempts', 8), ('max_cycles', 10000), ('max_evaluations', 10000)):
            value = getattr(self, key)
            if type(value) is not int or not 1 <= value <= maximum:
                raise ValueError(f'{key} must be an integer between 1 and {maximum}')
        if type(self.wall_seconds) not in (int, float) or not math.isfinite(self.wall_seconds) or not 0 < self.wall_seconds <= 86400:
            raise ValueError('Recovery wall_seconds must be finite, positive and at most 86400')


def validate_settings(settings):
    from pyar.backend_capabilities import optimization_recovery_capabilities
    if (not optimization_recovery_capabilities(settings.get('software'), settings.get('geometry_optimizer'))
            or settings.get('xtb_model', 'gfn2') != 'gfn2'
            or settings.get('opt_target', 'minimum') != 'minimum'
            or settings.get('multiplicity', 1) != 1 or settings.get('scftype', 'rhf') != 'rhf'):
        raise ValueError('recovery_unsupported: SAFE pilot requires closed-shell geomeTRIC/xTB GFN2 minima')
    forbidden = ('constraints', 'gamma', 'microsolvation_wall', 'custom_keywords', 'custom_keyword',
                 'bias_controller_restart', 'release_retry_attempt', 'microsolvation_evaluation_budget_path')
    if any(settings.get(key) for key in forbidden):
        raise ValueError('recovery_unsupported: bias, constraints, custom controls or workflow restart')
    # Unknown controls may carry hidden physical or restart semantics.
    allowed = {'software', 'geometry_optimizer', 'xtb_model', 'opt_target', 'charge', 'multiplicity',
        'scftype', 'xtb_unpaired_electrons', 'fragments', 'nprocs', 'opt_threshold', 'opt_cycles',
        'scf_cycles', 'scf_threshold', 'method', 'basis', '_electronic_state_from_molecule',
        '_ore_safe', '_ore_evaluation_limit', '_ore_initial_trust', '_ore_recovery_action',
        'workflow_role', *forbidden}
    if set(settings) - allowed:
        raise ValueError('recovery_unsupported: unverified controls: ' + ', '.join(sorted(set(settings)-allowed)))
    if settings.get('method') is not None or settings.get('basis') is not None:
        raise ValueError('recovery_unsupported: xTB does not consume method/basis controls')
    trust = settings.get('_ore_initial_trust')
    if trust is not None and (type(trust) not in (int, float) or not math.isfinite(trust) or not .005 <= trust <= .1):
        raise ValueError('recovery_unsupported: invalid bounded geomeTRIC trust radius')
    if settings.get('_ore_recovery_action') not in (None, 'trust_hessian_restart'):
        raise ValueError('recovery_unsupported: unrecognized internal recovery action')
    for key in ('opt_cycles', 'nprocs'):
        value = settings.get(key, 100 if key == 'opt_cycles' else 1)
        if type(value) is not int or value < 1:
            raise ValueError(key + ' must be a positive integer')
    if settings.get('workflow_role', 'standalone_minimum') != 'standalone_minimum':
        raise ValueError('recovery_unsupported: standalone minima only')


def improving(steps):
    accepted = [step for step in steps if step.accepted is True]
    if len(accepted) < 2:
        return False
    first, last = accepted[0], accepted[-1]
    # Both independent residual and energy must improve above logging precision.
    return (last.effective_energy_hartree < first.effective_energy_hartree - 1e-9
        and last.gradient_rms_hartree_per_bohr is not None
        and first.gradient_rms_hartree_per_bohr is not None
        and last.gradient_rms_hartree_per_bohr < first.gradient_rms_hartree_per_bohr * .99)


def flat_pes_evidence(steps, requested_rms_threshold, *, window=4,
                      energy_noise_hartree=1e-7, gradient_relative_span=.05):
    """Conservatively flag a plateau without relaxing requested tolerances."""
    if (type(window) is not int or window < 2
            or not math.isfinite(requested_rms_threshold) or requested_rms_threshold <= 0
            or not math.isfinite(energy_noise_hartree) or energy_noise_hartree < 0
            or not math.isfinite(gradient_relative_span) or gradient_relative_span < 0):
        raise ValueError('Flat-PES evidence thresholds must be finite and valid')
    accepted = [step for step in steps if step.accepted is True]
    recent = accepted[-window:]
    if (len(recent) < window or requested_rms_threshold <= 0
            or any(step.gradient_rms_hartree_per_bohr is None for step in recent)):
        return False
    gradients = [step.gradient_rms_hartree_per_bohr for step in recent]
    if min(gradients) <= requested_rms_threshold or min(gradients) <= 0:
        return False
    energy_span = max(step.effective_energy_hartree for step in recent) - min(
        step.effective_energy_hartree for step in recent)
    gradient_span = (max(gradients) - min(gradients)) / max(gradients)
    return energy_span <= energy_noise_hartree and gradient_span <= gradient_relative_span
