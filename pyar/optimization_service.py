"""Opt-in, single-execution diagnosis; existing workflows/CLIs keep legacy APIs.

Inspection is read-only. optimize_once runs the legacy adapter exactly once on
an independent molecule in a fresh job directory, then validates native evidence.
No recovery, cache adoption, workflow migration or chemistry change occurs here.
"""
from copy import deepcopy
from pathlib import Path
from time import perf_counter
from uuid import uuid4
import logging

from pyar.persistence import RunLock, LockError
from pyar.backend_capabilities import optimization_diagnostic_capabilities
from pyar.recovery.models import (
    ExecutionObservation, FailureCode, FailureDiagnostic, GateEvidence, GateState,
    OptimizationAttempt, OptimizationRequest, OptimizationResult, Outcome, legacy_diagnostic,
)
from pyar.recovery.observations import capture_processes, safe_path, snapshot_artifacts
from pyar.recovery.parsers import parse_observation
from pyar.recovery.validation import verify_start
from pyar.persistence import write_json


def _invalid(code, message, *, execution=None, exception=None, request=None, artifacts=()):
    return OptimizationResult(Outcome.INVALID_RESULT,
        (FailureDiagnostic(code, 'validator', message),),
        tuple(GateEvidence(axis, GateState.UNKNOWN, 'Not established by available evidence')
              for axis in ('backend', 'scf', 'optimizer', 'stationarity', 'objective', 'constraints', 'workflow')),
        (OptimizationAttempt(1, execution, request=request, artifacts=artifacts),), exception=exception)


def _unsupported_request(request):
    settings = request.settings
    try:
        supported = bool(optimization_diagnostic_capabilities(request.backend, request.optimizer))
    except ValueError:
        supported = False
    return (not supported or request.target != 'minimum' or bool(settings.get('constraints'))
            or (request.optimizer == 'native' and bool(settings.get('gamma') or settings.get('microsolvation_wall'))))


def inspect_optimization(request, directory, *, execution=None, exception=None):
    """Independently inspect pilot outputs; a receipt binds them to one execution.

    Missing provenance, unknown formats and unsupported workflow constraints fail
    closed. Returned ArtifactRefs describe observations, not committed restart
    authority. The original exception object and its chain remain accessible.
    """
    root = Path(directory).absolute()
    artifacts = ()
    try:
        artifacts = snapshot_artifacts(root)
        trusted = (execution is not None and execution.fresh_directory
                   and execution.request_fingerprint == request.request_fingerprint
                   and execution.artifacts == artifacts)
        if _unsupported_request(request):
            from pyar.recovery.diagnose_only import diagnose_only
            return diagnose_only(request, root, artifacts, execution=execution, exception=exception)
        known = {ref.path for ref in artifacts}
        def path(relative):
            result = safe_path(root, relative)
            if result.exists() and relative not in known:
                raise ValueError('Unobserved artifact: ' + relative)
            return result
        if trusted:
            verify_start(request, root, path)
        observed = parse_observation(request, path)
        if snapshot_artifacts(root) != artifacts:
            raise ValueError('Artifacts changed during observational inspection')
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return _invalid(FailureCode.INVALID_RESULT, str(exc), execution=execution,
                        exception=exception or exc, request=request, artifacts=artifacts)
    processes = execution.processes if trusted else ()
    returncode = processes[-1].returncode if processes else None
    terminated = observed.terminated and returncode == 0
    diagnostics = list(observed.diagnostics)
    if trusted and execution.legacy_status is not None and execution.legacy_status is not True and execution.legacy_status != 'converged':
        diagnostics.append(legacy_diagnostic(execution.legacy_status))
    if not trusted:
        diagnostics.append(FailureDiagnostic(FailureCode.CHECKPOINT_INCOMPATIBLE, 'validator',
            'Historical/unbound evidence: missing or mismatched fresh execution receipt', certainty='unverified_context'))
    if returncode not in (None, 0) and not observed.limited:
        diagnostics.append(FailureDiagnostic(FailureCode.EXEC_FAILED, 'environment', f'Process exited with {returncode}'))
    if observed.limited and not any(d.code is FailureCode.OPT_MAX_STEPS for d in diagnostics):
        diagnostics.append(FailureDiagnostic(FailureCode.OPT_MAX_STEPS, 'optimizer', 'Optimizer iteration limit reached'))
    if observed.frame is None or not observed.accepted:
        diagnostics.append(FailureDiagnostic(FailureCode.ARTIFACT_INCOMPLETE, 'validator', 'No independently verified accepted frame'))
    if not observed.converged and not diagnostics:
        diagnostics.append(FailureDiagnostic(FailureCode.UNKNOWN, 'optimizer', 'No genuine optimizer convergence evidence'))
    if exception is not None:
        diagnostics.append(FailureDiagnostic(FailureCode.EXEC_FAILED, 'legacy_adapter', str(exception)))
    scf_failed = any(d.code is FailureCode.CALC_SCF_MAX_STEPS for d in diagnostics)
    gates = (
        GateEvidence('backend', GateState.PASSED if terminated else GateState.UNKNOWN if returncode is None else GateState.FAILED, 'Native termination and observed process return code'),
        GateEvidence('scf', GateState.FAILED if scf_failed else GateState.UNKNOWN, 'SCF errors parsed; electronic root identity is not independently proven'),
        GateEvidence('optimizer', GateState.PASSED if observed.converged else GateState.FAILED, 'Native optimizer convergence evidence'),
        GateEvidence('stationarity', GateState.PASSED if observed.stationary else GateState.FAILED, 'Measured gradients must meet original requested thresholds'),
        GateEvidence('objective', GateState.PASSED if trusted and observed.frame is not None else GateState.UNKNOWN, 'Frame-associated physical and effective energies validated separately'),
        GateEvidence('constraints', GateState.NOT_APPLICABLE, 'No physical constraints requested'),
        GateEvidence('workflow', GateState.UNKNOWN, 'Workflow-specific identity, Hessian and pathway gates remain owner-controlled'),
    )
    if not trusted:
        status = Outcome.INVALID_RESULT
    elif observed.converged and observed.stationary and observed.accepted and terminated and not diagnostics:
        status = Outcome.CONVERGED
    elif scf_failed or any(d.code is FailureCode.USER_INPUT_INVALID for d in diagnostics) or (returncode not in (None, 0) and not observed.limited):
        status = Outcome.FAILED_FINAL
    elif observed.limited and observed.accepted and exception is None:
        status = Outcome.INCOMPLETE_USABLE
    elif observed.limited:
        status = Outcome.FAILED_RECOVERABLE
    else:
        status = Outcome.INVALID_RESULT
    if status is not Outcome.CONVERGED and not diagnostics:
        diagnostics.append(FailureDiagnostic(FailureCode.INVALID_RESULT, 'validator', 'Independent convergence validation failed'))
    attempt = OptimizationAttempt(1, execution, ((request.optimizer if request.optimizer == 'geometric' else request.backend, observed.version),), observed.steps, request, artifacts)
    return OptimizationResult(status, tuple(diagnostics), gates, (attempt,),
                              observed.frame if observed.accepted and trusted else None, observed.frame, exception)


def optimize_once(molecule, qc_params):
    """Execute one opt-in minimum optimization in the caller's current directory.

    Uses the existing optimizer's cwd behavior, never adds a new chdir layer.
    Refuses all existing job directories and obtains the existing RunLock. The
    caller molecule/settings remain untouched. No files are created on inspection.
    """
    from pyar import optimiser
    request = OptimizationRequest.from_molecule(molecule, qc_params)
    if _unsupported_request(request):
        return _invalid(FailureCode.BACKEND_UNSUPPORTED, 'R1 pilot does not validate selected route, target or physical controls', request=request)
    root = Path.cwd() / ('job_' + request.name)
    if root.exists() or root.is_symlink():
        return _invalid(FailureCode.CHECKPOINT_INCOMPATIBLE, 'Existing job directory cannot be adopted by R1', request=request)
    original = None
    try:
        with RunLock(root):
            if root.exists() or root.is_symlink():
                return _invalid(FailureCode.CHECKPOINT_INCOMPATIBLE, 'Job directory appeared before execution', request=request)
            candidate = deepcopy(molecule)
            params = request.settings
            params.pop('fragments', None)
            status, original = None, None
            with capture_processes() as processes:
                started = perf_counter()
                try:
                    status = optimiser.optimise(candidate, params)
                except (Exception, SystemExit) as exc:
                    original = exc
            try:
                execution = ExecutionObservation(request.request_fingerprint, snapshot_artifacts(root),
                                                 tuple(processes), True, status, elapsed_seconds=perf_counter()-started)
            except (OSError, ValueError) as exc:
                return _invalid(FailureCode.INVALID_RESULT, str(exc), exception=original or exc, request=request)
            return inspect_optimization(request, root, execution=execution, exception=original)
    except (OSError, LockError) as exc:
        return _invalid(FailureCode.EXEC_FAILED, str(exc), exception=original or exc, request=request)


def optimize_workflow_once(molecule, qc_params, *, workflow_role, optimizer=None,
                           allow_incomplete=False):
    """Run a workflow-owned optimization once and persist independent ORE evidence.

    The owning workflow continues to decide whether a legacy incomplete status
    is usable. This adapter records the independent result beside the workflow
    execution and does not retry, select a geometry, or alter legacy status.
    """
    from pyar import optimiser

    if not isinstance(workflow_role, str) or not workflow_role.strip():
        raise ValueError('workflow_role must be a nonempty string')
    if type(allow_incomplete) is not bool:
        raise ValueError('allow_incomplete must be an explicit Boolean')
    qc_params = dict(qc_params)
    workflow_recovery = qc_params.pop('workflow_recovery', None)
    request_error = None
    try:
        request = OptimizationRequest.from_molecule(molecule, qc_params)
        from dataclasses import replace
        request = replace(request, workflow_role=workflow_role)
    except (ValueError, TypeError, KeyError, AttributeError, OSError) as exc:
        # Preserve legacy workflow behavior for duck-typed compatibility objects;
        # their incomplete request metadata cannot be treated as validated.
        request, request_error = None, exc
    request_name = request.name if request is not None else str(getattr(molecule, 'name', 'unnamed'))
    execution_root = Path.cwd() / ('job_' + request_name)
    fresh = not execution_root.exists() and not execution_root.is_symlink()
    original = None
    status = None
    recovery_report = None
    started = perf_counter()
    try:
        with capture_processes() as processes:
            try:
                if workflow_recovery is not None:
                    from pyar.recovery.controller import run_recovery
                    from pyar.recovery.policy import RecoveryPolicy
                    if not isinstance(workflow_recovery, dict):
                        raise ValueError('workflow_recovery must contain explicit bounded policy settings')
                    controls = dict(workflow_recovery)
                    resume = controls.pop('resume', False)
                    if type(resume) is not bool:
                        raise ValueError('workflow_recovery.resume must be Boolean')
                    policy = RecoveryPolicy(**controls)
                    if policy.mode == 'off':
                        raise ValueError('Omit workflow_recovery to retain the legacy execution path')
                    recovery_report = run_recovery(molecule, dict(qc_params, workflow_role=workflow_role),
                                                   policy, resume=resume)
                    recovered = recovery_report.result
                    status = False
                    if recovered.converged and recovered.accepted_frame is not None:
                        import numpy as np
                        frame = recovered.accepted_frame
                        molecule.coordinates = np.array(frame.coordinates_angstrom)
                        molecule.energy = frame.physical_energy_hartree
                        status = True
                else:
                    status = (optimizer or optimiser.optimise)(molecule, qc_params)
            except Exception as exc:
                original = exc
        elapsed = max(0., perf_counter() - started)
        if recovery_report is not None:
            result = recovery_report.result
        elif request is None:
            result = _invalid(FailureCode.INVALID_RESULT,
                'Workflow optimization request could not be normalized: ' + str(request_error),
                exception=original or request_error)
        else:
            artifacts = snapshot_artifacts(execution_root)
            execution = ExecutionObservation(request.request_fingerprint, artifacts,
                tuple(processes), fresh, status, elapsed_seconds=elapsed)
            result = inspect_optimization(request, execution_root, execution=execution,
                                          exception=original)
    except Exception as exc:
        result = _invalid(FailureCode.INVALID_RESULT, str(exc), exception=original or exc,
                          request=request)
    from pyar.recovery.workflow import acceptance
    result, decision = acceptance(result, status, allow_incomplete=allow_incomplete)
    record = {
        'schema_version': 1,
        'workflow_role': workflow_role,
        'molecule_name': request_name,
        'request_error': None if request_error is None else str(request_error),
        'legacy_status': optimiser._status_label(status),
        'result': result.to_dict(),
        'owner_decision': decision,
        'recovery': None if recovery_report is None else recovery_report.to_dict(),
    }
    report = Path.cwd() / '.pyar_optimization_diagnostics' / f'{uuid4().hex}.json'
    try:
        write_json(report, record, overwrite=False)
    except Exception as exc:
        logging.getLogger(__name__).warning('Could not persist optimization diagnostics: %s', exc)
        report = None
    if original is not None:
        try:
            original.ore_diagnostic_path = str(report) if report is not None else None
        except Exception:
            pass
        raise original
    try:
        molecule.optimization_owner_decision = decision
    except (AttributeError, TypeError):
        pass
    return status, result, report
