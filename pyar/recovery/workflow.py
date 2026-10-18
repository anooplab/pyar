"""Workflow-owned acceptance and stage evidence, independent of restart authority."""
from dataclasses import replace

from .models import FailureCode, GateEvidence, GateState, _exception_chain


def acceptance(result, status, *, allow_incomplete):
    """Make the legacy owner's provisional-result permission explicit.

    Acceptance is a permission to continue to workflow identity/selection gates,
    never a numerical convergence certificate. Existing rankings remain owned by
    each workflow. Unknown/diagnose-only evidence is always labelled provisional.
    """
    from pyar.optimiser import is_success, is_usable
    permitted = is_usable(status) if allow_incomplete else is_success(status)
    detail = ('Owner permits incomplete results for subsequent workflow validation'
              if allow_incomplete else 'Owner requires a legacy converged status')
    decision = dict(accepted=bool(permitted), allow_incomplete=allow_incomplete,
                    validated_converged=bool(permitted and result.converged),
                    provisional=bool(permitted and not result.converged), reason=detail)
    # Leave the final chemistry/workflow gate unknown: only the owner can pass it
    # after connectivity, shell, ranking or path checks have actually executed.
    gates = tuple(g for g in result.gates if g.axis != 'acceptance_policy') + (
        GateEvidence('acceptance_policy', GateState.PASSED if permitted else GateState.FAILED, detail),)
    return replace(result, gates=gates), decision


def stage_diagnostic(role, gates, *, exception=None, reason=None):
    """Record separate owner gates; a TS optimizer boolean cannot prove a saddle."""
    states = {key: 'passed' if value is True else 'failed' if value is False else 'unknown'
              for key, value in gates.items()}
    failed = exception is not None or any(state == 'failed' for state in states.values())
    code = (FailureCode.EXEC_FAILED.value if exception is not None else
            FailureCode.INVALID_RESULT.value if failed else
            'OWNER_GATES_PASSED' if states and all(s == 'passed' for s in states.values()) else
            FailureCode.UNKNOWN.value)
    return dict(schema_version=1, workflow_role=role, stage_gate=(
        'failed' if failed else 'passed' if code == 'OWNER_GATES_PASSED' else 'not_assessed'),
        gates=states, reason=reason or code, exception_chain=_exception_chain(exception),
        automatic_recovery_attempted=False, recovery_supported=False,
        recovery_stop_reason='No certified restart contract for this stage',
        restart_authority=False)
