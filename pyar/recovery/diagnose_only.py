"""Negative evidence for routes without a certified accepted-frame parser.

These observations never certify convergence or authorize restart. Backend
markers are negative evidence only; a normal exit or a result XYZ is insufficient.
"""
import re

from .models import (FailureCode, FailureDiagnostic, GateEvidence, GateState,
                     OptimizationAttempt, OptimizationResult, Outcome)
from .observations import safe_path


BACKEND_FAILURES = {
    'gaussian': ((r'Convergence failure|SCF has not converged', FailureCode.CALC_SCF_MAX_STEPS),
                 (r'Number of steps exceeded', FailureCode.OPT_MAX_STEPS)),
    'psi4': ((r'Could not converge SCF|SCFConvergenceError', FailureCode.CALC_SCF_MAX_STEPS),
             (r'Geometry optimization.*failed|maximum.*optimization.*steps', FailureCode.OPT_MAX_STEPS)),
    'turbomole': ((r'SCF NOT CONVERGED|scf convergence not reached', FailureCode.CALC_SCF_MAX_STEPS),),
    'xtb_turbo': ((r'SCF NOT CONVERGED|SCC.*not converged', FailureCode.CALC_SCF_MAX_STEPS),),
    'mopac': ((r'SCF FIELD WAS NOT ACHIEVED|UNABLE TO ACHIEVE SELF-CONSISTENCE', FailureCode.CALC_SCF_MAX_STEPS),
              (r'EXCESS NUMBER OF OPTIMIZATION CYCLES', FailureCode.OPT_MAX_STEPS)),
    'aimnet_2': ((r'CUDA out of memory|OutOfMemoryError', FailureCode.EXEC_FAILED),
                 (r'non-finite|NaN.*(?:force|gradient)', FailureCode.CALC_INVALID_GRADIENT)),
    'ani': ((r'CUDA out of memory|OutOfMemoryError', FailureCode.EXEC_FAILED),),
    'aiqm1_mlatom': ((r'SCF.*(?:failed|not converged)|Convergence failure', FailureCode.CALC_SCF_MAX_STEPS),),
    'obabel': (), 'xtb-aimnet2': (), 'xtb-aiqm1': (),
}


def diagnose_only(request, root, artifacts, *, execution=None, exception=None):
    """Retain artifact/exception evidence without inventing a final frame."""
    diagnostics = []
    patterns = BACKEND_FAILURES.get(request.backend, ())
    # Composite stages must not borrow a successful first-stage certificate.
    if request.backend == 'xtb-aimnet2':
        patterns = BACKEND_FAILURES['aimnet_2']
    elif request.backend == 'xtb-aiqm1':
        patterns = BACKEND_FAILURES['aiqm1_mlatom']
    for artifact in artifacts:
        if not artifact.path.endswith(('.out', '.log', '.arc')):
            continue
        with safe_path(root, artifact.path).open(errors='replace') as stream:
            for number, line in enumerate(stream, 1):
                for pattern, code in patterns:
                    if re.search(pattern, line, re.I):
                        diagnostics.append(FailureDiagnostic(code, request.backend,
                            line.strip()[:300], artifact.path, number, line.strip()[:300],
                            certainty='observed_unvalidated_context'))
                        break
                if len(diagnostics) >= 32:
                    break
        if len(diagnostics) >= 32:
            break
    failed = exception is not None or bool(execution and any(
        p.returncode not in (None, 0) for p in execution.processes))
    if failed:
        diagnostics.append(FailureDiagnostic(FailureCode.EXEC_FAILED, 'execution',
            str(exception) if exception is not None else 'Observed nonzero process exit'))
    diagnostics.append(FailureDiagnostic(FailureCode.BACKEND_UNSUPPORTED, 'capability',
        'Independent numerical validation is unavailable for the selected route/objective'))
    diagnostics.append(FailureDiagnostic(FailureCode.RECOVERY_UNSUPPORTED, 'capability',
        'Diagnose-only route: accepted-frame, numerical convergence and restart compatibility are unverified'))
    gates = tuple(GateEvidence(axis, GateState.UNKNOWN, 'Not certified on this route')
                  for axis in ('backend', 'scf', 'optimizer', 'stationarity', 'objective', 'constraints', 'workflow'))
    return OptimizationResult(Outcome.FAILED_FINAL if failed else Outcome.INVALID_RESULT,
        tuple(diagnostics), gates,
        (OptimizationAttempt(1, execution, request=request, artifacts=artifacts),), exception=exception)
