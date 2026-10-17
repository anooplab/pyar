"""Bounded standalone SAFE cycle continuation, using existing persistence primitives.

Immutable intents precede execution; immutable attempt commits follow independent
validation. outcome.json is derived and is never used to authorize continuation.
A started attempt without a sealed execution receipt is never blindly replayed.
"""
from dataclasses import asdict, dataclass, replace
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time

from pyar.persistence import ArtifactRef, RunLock, canonical_json_bytes, read_json, write_json
from pyar.persistence.artifacts import fsync_artifacts
from pyar.optimization_service import inspect_optimization
from .models import (ExecutionObservation, FailureCode, FailureDiagnostic, OptimizationRequest,
                     OptimizationAttempt, OptimizationResult, Outcome, ProcessRecord, RecordedExecutionError)
from .observations import safe_path, snapshot_artifacts
from .policy import RecoveryPolicy, improving, validate_settings


@dataclass(frozen=True)
class RecoveryReport:
    result: OptimizationResult
    policy: RecoveryPolicy
    stop_reason: str
    reserved_cycles: int
    backend_evaluations: int
    elapsed_seconds: float

    def to_dict(self):
        return dict(schema_version=1, result=self.result.to_dict(), policy=asdict(self.policy),
                    stop_reason=self.stop_reason, reserved_cycles=self.reserved_cycles,
                    backend_evaluations=self.backend_evaluations, elapsed_seconds=self.elapsed_seconds)


def prepare_request(molecule, settings, policy):
    request = OptimizationRequest.from_molecule(molecule, dict(settings, _ore_safe=True))
    validate_settings(request.settings)
    if request.settings.get('opt_cycles', 100) > policy.max_cycles:
        raise ValueError('Initial opt_cycles exceeds total recovery cycle budget')
    return request


def preflight_recovery(molecules, settings, policy, *, resume=False):
    """No directories, locks, files, backend construction or subprocess execution."""
    requests = [prepare_request(mol, settings, policy) for mol in molecules]
    for request in requests:
        root = safe_path(Path.cwd(), 'job_' + request.name)
        if root.exists():
            if not resume:
                raise ValueError('Existing job directory; use explicit --recovery-resume for an ORE run')
            metadata = read_json(safe_path(root, 'recovery/run.json'))
            _validate_run(metadata, request, policy)
            try:
                _preflight_records(root, request, policy)
            except (KeyError, TypeError, IndexError, AttributeError) as exc:
                raise ValueError('Malformed recovery record') from exc
        elif resume:
            raise ValueError('Cannot resume a missing recovery run')
    return requests


def _validate_run(metadata, request, policy):
    if not isinstance(metadata, dict):
        raise ValueError('Invalid recovery metadata')
    # JSON normalization avoids tuple/list identity differences on reload.
    expected = dict(schema_version=1, request_fingerprint=request.request_fingerprint,
                    scientific_fingerprint=request.scientific_fingerprint, policy=asdict(policy))
    if any(metadata.get(k) != v for k, v in expected.items()):
        raise ValueError('Incompatible recovery request or policy; original science/budgets are immutable')
    if (type(metadata.get('created_epoch')) not in (int, float)
            or not math.isfinite(metadata['created_epoch'])
            or metadata.get('deadline_epoch') != metadata['created_epoch'] + policy.wall_seconds):
        raise ValueError('Invalid immutable walltime budget')


def _preflight_records(root, original, policy):
    recovery = safe_path(root, 'recovery')
    if read_json(safe_path(recovery, 'request.json')) != _encode(original):
        raise ValueError('Original request record changed')
    attempts = safe_path(recovery, 'attempts')
    if not attempts.exists():
        return
    children = sorted(attempts.iterdir())
    if len(children) > policy.max_attempts:
        raise ValueError('Attempt count exceeds immutable budget')
    for number, child in enumerate(children, 1):
        directory = safe_path(attempts, f'attempt_{number:03d}')
        if child != directory or not directory.is_dir():
            raise ValueError('Unrecognized/non-contiguous recovery attempt directory')
        intent = read_json(safe_path(directory, 'input.json'))
        current = _request(intent['request'])
        if (intent.get('schema_version') != 1 or intent.get('number') != number
                or intent.get('original_scientific_fingerprint') != original.scientific_fingerprint
                or current.objective_fingerprint != original.objective_fingerprint):
            raise ValueError('Incompatible attempt intent')
        result = safe_path(directory, 'result.json')
        commit = safe_path(directory, 'commit.json')
        if result.exists():
            value = read_json(result)
            if commit.exists():
                record = read_json(commit)
                if record.get('schema_version') != 1 or record.get('intent') != intent or record.get('execution_result') != value:
                    raise ValueError('Committed attempt changed')
            _validated_result(current, safe_path(directory, 'job_' + original.name), value)
        elif commit.exists():
            raise ValueError('Committed receipt is missing')


def _request(value):
    return OptimizationRequest(value['name'], tuple(value['symbols']), value['coordinates_angstrom'],
        value['backend'], value['optimizer'], value['target'], canonical_json_bytes(value['settings']).decode(),
        canonical_json_bytes(value['asset_hashes']).decode(), value['workflow_role'])


def _encode(request):
    # Reuse the schema-1 request serializer rather than a second settings schema.
    from .models import _serializable
    return json.loads(canonical_json_bytes(_serializable(request)))


def _validated_result(request, directory, value):
    if len(value['attempts']) != 1 or value['attempts'][0]['request'] != _encode(request):
        raise ValueError('Attempt execution request does not match immutable intent')
    record = value['attempts'][0]['execution']
    if record is None:
        raise ValueError('Missing execution receipt')
    execution = ExecutionObservation(record['request_fingerprint'],
        tuple(ArtifactRef(**ref) for ref in record['artifacts']),
        tuple(ProcessRecord(p['program'], p['returncode']) for p in record['processes']),
        record['fresh_directory'], record['legacy_status'], record['elapsed_seconds'])
    if execution.artifacts != snapshot_artifacts(directory):
        raise ValueError('Execution artifacts changed after receipt publication')
    exception = RecordedExecutionError(value['exception_chain']) if value['exception_chain'] else None
    return inspect_optimization(request, directory, execution=execution, exception=exception)


def _execute_attempt(directory, timeout):
    """Walltime cap covers the complete process group, including backend children."""
    with (directory / 'worker.out').open('x') as stream:
        process = subprocess.Popen([sys.executable, '-m', 'pyar.recovery.worker'], cwd=directory,
            stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            return process.wait(timeout=timeout)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise


def run_recovery(molecule, settings, policy=RecoveryPolicy(), *, resume=False, recover_lock=False):
    if policy.mode == 'off':
        from pyar.optimization_service import optimize_once
        return RecoveryReport(optimize_once(molecule, settings), policy, 'recovery_off', 0, 0, 0.)
    request = preflight_recovery([molecule], settings, policy, resume=resume)[0]
    root = safe_path(Path.cwd(), 'job_' + request.name)
    with RunLock(root, recover_stale=recover_lock):
        if root.exists() and not resume:
            raise ValueError('Job directory appeared before execution')
        recovery = safe_path(root, 'recovery')
        metadata_path = safe_path(recovery, 'run.json')
        if resume:
            metadata = read_json(metadata_path)
            _validate_run(metadata, request, policy)
        else:
            created = time.time()
            metadata = dict(schema_version=1, request_fingerprint=request.request_fingerprint,
                scientific_fingerprint=request.scientific_fingerprint, policy=asdict(policy),
                created_epoch=created, deadline_epoch=created+policy.wall_seconds)
            write_json(metadata_path, metadata, overwrite=False)
            write_json(safe_path(recovery, 'request.json'), _encode(request), overwrite=False)
        # Raw/accepted artifacts, never outcome.json, authorize replay.
        if read_json(safe_path(recovery, 'request.json')) != _encode(request):
            raise ValueError('Original request record changed')
        attempts = []
        used_cycles, evaluations = 0, 0
        start, settings_next = request.coordinates_angstrom, request.settings
        remaining = policy.max_evaluations
        signatures = set()
        last = OptimizationResult(Outcome.INVALID_RESULT, (), (), ())
        reason = 'attempt_budget_exhausted'
        for number in range(1, policy.max_attempts+1):
            directory = safe_path(recovery, f'attempts/attempt_{number:03d}')
            cycles = min(settings_next.get('opt_cycles', 100), policy.max_cycles-used_cycles)
            if cycles < 1 or remaining < 1:
                reason = 'cycle_or_evaluation_budget_exhausted'
                break
            current = replace(request, coordinates_angstrom=start,
                settings_json=canonical_json_bytes(dict(settings_next, opt_cycles=cycles,
                                                        _ore_evaluation_limit=remaining)).decode())
            if current.objective_fingerprint != request.objective_fingerprint:
                raise ValueError('Recovery attempted to change the scientific objective')
            signature = current.request_fingerprint
            if signature in signatures:
                reason = 'duplicate_action_suppressed'
                break
            signatures.add(signature)
            intent = dict(schema_version=1, number=number, parent=number-1 if number>1 else None,
                action='cycle_extension' if number>1 else 'initial',
                rationale='Accepted energy and force residual improved' if number>1 else 'Original request',
                original_scientific_fingerprint=request.scientific_fingerprint, request=_encode(current))
            intent_path = safe_path(directory, 'input.json')
            if directory.exists():
                if not resume or read_json(intent_path) != intent:
                    raise ValueError('Attempt intent missing or incompatible')
            else:
                if time.time() >= metadata['deadline_epoch']:
                    reason = 'wall_budget_exhausted'
                    break
                write_json(intent_path, intent, overwrite=False)
            used_cycles += cycles  # Reserve before launching; never refunded after a crash.
            backend = safe_path(directory, 'job_' + request.name)
            result_path = safe_path(directory, 'result.json')
            commit_path = safe_path(directory, 'commit.json')
            if commit_path.exists():
                committed = read_json(commit_path)
                if committed['intent'] != intent or read_json(result_path) != committed['execution_result']:
                    raise ValueError('Committed attempt changed')
            elif not result_path.exists():
                # A backend or worker.out means execution was started; never
                # repeat it from a possibly stale/partial file on restart.
                if backend.exists() or safe_path(directory, 'worker.out').exists():
                    reason = 'interrupted_attempt_without_receipt'
                    attempts.append(OptimizationAttempt(number, None, request=current,
                        action=intent['action'], parent_number=intent['parent'],
                        original_scientific_fingerprint=request.scientific_fingerprint))
                    evaluations += remaining  # No receipt: charge entire reserved quota.
                    break
                try:
                    code = _execute_attempt(directory, max(.001, metadata['deadline_epoch']-time.time()))
                except subprocess.TimeoutExpired:
                    reason = 'wall_budget_exhausted'
                    attempts.append(OptimizationAttempt(number, None, request=current,
                        action=intent['action'], parent_number=intent['parent'],
                        original_scientific_fingerprint=request.scientific_fingerprint))
                    evaluations += remaining
                    break
                if code != 0 or not result_path.exists():
                    reason = 'worker_failed_without_receipt'
                    attempts.append(OptimizationAttempt(number, None, request=current,
                        action=intent['action'], parent_number=intent['parent'],
                        original_scientific_fingerprint=request.scientific_fingerprint,
                        worker_returncode=code))
                    evaluations += remaining
                    break
            result_value = read_json(result_path)
            last = _validated_result(current, backend, result_value)
            if not commit_path.exists():
                fsync_artifacts([ref.resolve(backend) for ref in snapshot_artifacts(backend)], backend)
                write_json(commit_path, dict(schema_version=1, intent=intent,
                    execution_result=result_value), overwrite=False)
            counter = safe_path(backend, 'ore_evaluations.json')
            observed_evaluations = read_json(counter)['used'] if counter.exists() else remaining
            if type(observed_evaluations) is not int or not 0 <= observed_evaluations <= remaining:
                raise ValueError('Invalid evaluation accounting')
            evaluations += observed_evaluations
            remaining -= observed_evaluations
            attempt = replace(last.attempts[0], number=number, action=intent['action'],
                parent_number=intent['parent'], original_scientific_fingerprint=request.scientific_fingerprint,
                worker_returncode=0)
            if last.accepted_frame is not None:
                attempt = replace(attempt, versions=attempt.versions + (('xtb', '6.7.1'),))
            attempts.append(attempt)
            if last.converged:
                reason = 'validated_convergence'
                last = replace(last, status=Outcome.RECOVERED_CONVERGED if number>1 else Outcome.CONVERGED)
                break
            if remaining == 0:
                reason = 'evaluation_budget_exhausted' if counter.exists() else 'unverified_evaluation_cost'
                break
            if (last.status is not Outcome.INCOMPLETE_USABLE or last.accepted_frame is None
                    or {d.code for d in last.diagnostics} - {FailureCode.OPT_MAX_STEPS}
                    or not improving(attempt.steps)):
                reason = 'no_verified_safe_action'
                break
            start = last.accepted_frame.coordinates_angstrom
            if start == current.coordinates_angstrom:
                reason = 'unchanged_accepted_geometry'
                break
            settings_next = dict(current.settings, opt_cycles=cycles*2)
        if not last.converged:
            last = replace(last, diagnostics=last.diagnostics+(FailureDiagnostic(
                FailureCode.RECOVERY_STOPPED, 'recovery', reason),))
        last = replace(last, attempts=tuple(attempts))
        elapsed = max(0., time.time()-metadata['created_epoch'])
        report = RecoveryReport(last, policy, reason, used_cycles, evaluations, elapsed)
        if last.converged:
            # Derived legacy filename, never restart authority. Only verified
            # convergence may publish it; source native attempts remain immutable.
            frame = last.accepted_frame
            filename = safe_path(root, 'result_' + request.name + '.xyz')
            payload = f'{len(frame.symbols)}\n{frame.physical_energy_hartree:.16g}\n' + ''.join(
                f'{s} {x:.12f} {y:.12f} {z:.12f}\n' for s,(x,y,z) in zip(frame.symbols,frame.coordinates_angstrom))
            if filename.exists() and filename.read_text() != payload:
                raise ValueError('Existing final result differs from validated outcome')
            if not filename.exists():
                temporary = None
                try:
                    with tempfile.NamedTemporaryFile(mode='w', dir=root, prefix='.ore-final-', delete=False) as stream:
                        temporary = Path(stream.name)
                        stream.write(payload)
                        stream.flush()
                        os.fsync(stream.fileno())
                    os.replace(temporary, filename)
                    from pyar.persistence.json_store import fsync_directory
                    fsync_directory(root)
                finally:
                    if temporary is not None:
                        temporary.unlink(missing_ok=True)
        write_json(safe_path(recovery, 'outcome.json'), report.to_dict())
        return report
