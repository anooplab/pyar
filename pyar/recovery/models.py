"""Small immutable contracts for observational optimization diagnostics (schema 1)."""
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
import hashlib
import json
import math
from numbers import Integral
from pathlib import Path

from pyar.persistence import ArtifactRef, canonical_json_bytes


class Outcome(str, Enum):
    CONVERGED = 'CONVERGED'
    RECOVERED_CONVERGED = 'RECOVERED_CONVERGED'
    INCOMPLETE_USABLE = 'INCOMPLETE_USABLE'
    FAILED_RECOVERABLE = 'FAILED_RECOVERABLE'
    FAILED_FINAL = 'FAILED_FINAL'
    INVALID_RESULT = 'INVALID_RESULT'


class FailureCode(str, Enum):
    OPT_MAX_STEPS = 'OPT_MAX_STEPS'
    CALC_SCF_MAX_STEPS = 'CALC_SCF_MAX_STEPS'
    CALC_INVALID_GRADIENT = 'CALC_INVALID_GRADIENT'
    CONSTRAINT_VIOLATION = 'CONSTRAINT_VIOLATION'
    CHECKPOINT_INCOMPATIBLE = 'CHECKPOINT_INCOMPATIBLE'
    ARTIFACT_INCOMPLETE = 'ARTIFACT_INCOMPLETE'
    INVALID_RESULT = 'INVALID_RESULT'
    BACKEND_UNSUPPORTED = 'BACKEND_UNSUPPORTED'
    EXEC_FAILED = 'EXEC_FAILED'
    USER_INPUT_INVALID = 'USER_INPUT_INVALID'
    UNKNOWN = 'UNKNOWN'
    RECOVERY_UNSUPPORTED = 'RECOVERY_UNSUPPORTED'
    RECOVERY_STOPPED = 'RECOVERY_STOPPED'


class GateState(str, Enum):
    PASSED = 'passed'
    FAILED = 'failed'
    UNKNOWN = 'unknown'
    NOT_APPLICABLE = 'not_applicable'


def matrix(value, count):
    """Freeze an N by 3 finite matrix, including N=1, without scientific imports."""
    try:
        rows = tuple(tuple(float(x) for x in row) for row in value)
    except (TypeError, ValueError) as exc:
        raise ValueError('Coordinates/vectors must have shape (N, 3)') from exc
    if len(rows) != count or any(len(row) != 3 for row in rows):
        raise ValueError('Coordinates/vectors must have shape (N, 3)')
    if not all(math.isfinite(x) for row in rows for x in row):
        raise ValueError('Coordinates/vectors must be finite')
    return rows


def fingerprint(value):
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


@dataclass(frozen=True)
class OptimizationRequest:
    name: str
    symbols: tuple[str, ...]
    coordinates_angstrom: tuple[tuple[float, ...], ...]
    backend: str
    optimizer: str
    target: str
    settings_json: str
    asset_hashes_json: str | None = None
    workflow_role: str = 'standalone_minimum'

    def __post_init__(self):
        if (not self.name or self.name in {'.', '..'} or any(c in self.name for c in '/\\\x00')):
            raise ValueError('Optimization name must be a filename stem')
        object.__setattr__(self, 'symbols', tuple(self.symbols))
        from pyar.data.new_atomic_data import atomic_number
        if not self.symbols or any(not isinstance(s, str) or s not in atomic_number for s in self.symbols):
            raise ValueError('Ordered atom symbols are required')
        object.__setattr__(self, 'coordinates_angstrom', matrix(self.coordinates_angstrom, len(self.symbols)))
        settings = json.loads(self.settings_json)
        if not isinstance(settings, dict):
            raise ValueError('Settings must be a JSON object')
        object.__setattr__(self, 'settings_json', canonical_json_bytes(settings).decode())
        assets = {}
        if self.asset_hashes_json is None:
            from pyar.persistence import sha256_file
            for key in ('model', 'gxtb_wrapper'):
                if settings.get(key) is not None:
                    assets[key] = sha256_file(Path(settings[key]).expanduser())
        else:
            assets = json.loads(self.asset_hashes_json)
        object.__setattr__(self, 'asset_hashes_json', canonical_json_bytes(assets).decode())

    @property
    def settings(self):
        return json.loads(self.settings_json)

    @classmethod
    def from_molecule(cls, molecule, settings):
        from pyar.backend_capabilities import normalize_backend_name
        values = dict(settings)
        for key in ('charge', 'multiplicity'):
            value = getattr(molecule, key)
            if isinstance(value, bool) or not isinstance(value, Integral) or (key == 'multiplicity' and value < 1):
                raise ValueError('Charge must be integer and multiplicity must be positive integer')
            if key in values and values[key] != value:
                raise ValueError('Conflicting molecule and QC electronic state: ' + key)
        if 'scftype' in values and values['scftype'] != molecule.scftype:
            raise ValueError('Conflicting molecule and QC reference state')
        backend = normalize_backend_name(values['software'])
        values.update(software=backend, charge=int(molecule.charge), multiplicity=int(molecule.multiplicity),
                      scftype=molecule.scftype, fragments=[list(f) for f in getattr(molecule, 'fragments', [])])
        if backend == 'xtb':
            from pyar.backends.xtb_utils import canonical_xtb_model
            values['xtb_model'] = canonical_xtb_model(values.get('xtb_model'))
            unpaired = int(molecule.multiplicity) - 1
            if values.get('xtb_unpaired_electrons', unpaired) != unpaired:
                raise ValueError('Conflicting xTB unpaired-electron count')
            values['xtb_unpaired_electrons'] = unpaired
        return cls(molecule.name, tuple(molecule.atoms_list), molecule.coordinates, backend,
                   values.get('geometry_optimizer', 'native'), values.get('opt_target', 'minimum'),
                   canonical_json_bytes(values).decode(), workflow_role=values.get('workflow_role', 'standalone_minimum'))

    def _identity(self, *, scientific=False):
        values = self.settings
        # Only execution budgets/parallelism vary independently of requested science.
        # Requested tolerances and unknown options stay in the scientific identity.
        if scientific:
            for key in ('nprocs', 'opt_cycles', 'scf_cycles', '_ore_evaluation_limit', '_ore_safe'):
                values.pop(key, None)
        # Known file-backed models/wrappers use content identity, not host paths.
        for key in ('model', 'gxtb_wrapper'):
            if values.get(key) is not None:
                values[key] = {'sha256': json.loads(self.asset_hashes_json)[key]}
        return dict(symbols=self.symbols, coordinates_angstrom=self.coordinates_angstrom,
                    backend=self.backend, optimizer=self.optimizer, target=self.target, workflow_role=self.workflow_role, settings=values)

    @property
    def scientific_fingerprint(self):
        return fingerprint(self._identity(scientific=True))

    @property
    def request_fingerprint(self):
        return fingerprint(self._identity())

    @property
    def objective_fingerprint(self):
        """Science independent of the explicitly recorded starting geometry."""
        identity = self._identity(scientific=True)
        identity.pop('coordinates_angstrom')
        return fingerprint(identity)


@dataclass(frozen=True)
class FailureDiagnostic:
    code: FailureCode
    source: str
    message: str
    path: str | None = None
    line: int | None = None
    excerpt: str | None = None
    certainty: str = 'observed'
    severity: str = 'error'


@dataclass(frozen=True)
class GateEvidence:
    axis: str
    state: GateState
    detail: str


@dataclass(frozen=True)
class ProcessRecord:
    program: str
    returncode: int | None
    exception: BaseException | None = None

    def __post_init__(self):
        if self.returncode is not None and type(self.returncode) is not int:
            raise ValueError('Process return code must be an integer or unknown')


@dataclass(frozen=True)
class ExecutionObservation:
    """Caller-owned receipt, not persisted restart authority or a workflow commit."""
    request_fingerprint: str
    artifacts: tuple[ArtifactRef, ...]
    processes: tuple[ProcessRecord, ...]
    fresh_directory: bool
    legacy_status: bool | str | None = None
    elapsed_seconds: float | None = None

    def __post_init__(self):
        if self.elapsed_seconds is not None and (not math.isfinite(self.elapsed_seconds) or self.elapsed_seconds < 0):
            raise ValueError('Observed execution time must be finite and nonnegative')
        if type(self.fresh_directory) is not bool:
            raise ValueError('Fresh execution attestation must be explicit Boolean')
        if (not isinstance(self.request_fingerprint, str) or len(self.request_fingerprint) != 64
                or any(c not in '0123456789abcdef' for c in self.request_fingerprint)):
            raise ValueError('Invalid execution request fingerprint')
        object.__setattr__(self, 'artifacts', tuple(self.artifacts))
        object.__setattr__(self, 'processes', tuple(self.processes))


@dataclass(frozen=True)
class GeometryFrame:
    symbols: tuple[str, ...]
    coordinates_angstrom: tuple[tuple[float, ...], ...]
    physical_energy_hartree: float
    effective_energy_hartree: float | None = None
    effective_forces_ev_per_angstrom: tuple[tuple[float, ...], ...] | None = None
    source_paths: tuple[str, ...] = ()

    def __post_init__(self):
        object.__setattr__(self, 'symbols', tuple(self.symbols))
        object.__setattr__(self, 'coordinates_angstrom', matrix(self.coordinates_angstrom, len(self.symbols)))
        if self.physical_energy_hartree is None:
            raise ValueError('Physical energy is required')
        for value in (self.physical_energy_hartree, self.effective_energy_hartree):
            if value is not None and not math.isfinite(float(value)):
                raise ValueError('Energy must be a finite scalar')
        if self.effective_forces_ev_per_angstrom is not None:
            object.__setattr__(self, 'effective_forces_ev_per_angstrom', matrix(self.effective_forces_ev_per_angstrom, len(self.symbols)))


@dataclass(frozen=True)
class StepObservation:
    index: int
    effective_energy_hartree: float
    gradient_rms_hartree_per_bohr: float | None = None
    gradient_max_hartree_per_bohr: float | None = None
    accepted: bool | None = None
    gradient_norm_hartree_per_bohr: float | None = None

    def __post_init__(self):
        if type(self.index) is not int or self.index < 0:
            raise ValueError('Step index must be a nonnegative integer')
        if self.accepted is not None and type(self.accepted) is not bool:
            raise ValueError('Step acceptance must be explicit Boolean or unknown')
        if not math.isfinite(self.effective_energy_hartree):
            raise ValueError('Step energy must be finite')
        for value in (self.gradient_rms_hartree_per_bohr, self.gradient_max_hartree_per_bohr,
                      self.gradient_norm_hartree_per_bohr):
            if value is not None and (not math.isfinite(value) or value < 0):
                raise ValueError('Gradient residuals must be finite and nonnegative')


@dataclass(frozen=True)
class OptimizationAttempt:
    number: int
    execution: ExecutionObservation | None
    versions: tuple[tuple[str, str], ...] = ()
    steps: tuple[StepObservation, ...] = ()
    request: OptimizationRequest | None = None
    artifacts: tuple[ArtifactRef, ...] = ()
    action: str = 'initial'
    parent_number: int | None = None
    original_scientific_fingerprint: str | None = None
    worker_returncode: int | None = None


def _exception_chain(exception):
    if isinstance(exception, RecordedExecutionError):
        return list(exception.recorded_chain)
    chain, seen, exc = [], set(), exception
    while exc is not None and id(exc) not in seen:
        seen.add(id(exc))
        chain.append({'type': type(exc).__name__, 'message': str(exc)})
        exc = exc.__cause__ or (None if exc.__suppress_context__ else exc.__context__)
    return chain


class RecordedExecutionError(RuntimeError):
    """Persisted evidence; does not pretend to recreate a live exception type."""
    def __init__(self, chain):
        self.recorded_chain = tuple(dict(item) for item in chain)
        super().__init__(self.recorded_chain[0]['message'] if chain else 'Recorded execution failure')


def _serializable(value):
    if isinstance(value, OptimizationRequest):
        return dict(name=value.name, symbols=value.symbols, coordinates_angstrom=value.coordinates_angstrom,
                    backend=value.backend, optimizer=value.optimizer, target=value.target, workflow_role=value.workflow_role,
                    settings=value.settings, asset_hashes=json.loads(value.asset_hashes_json),
                    scientific_fingerprint=value.scientific_fingerprint, request_fingerprint=value.request_fingerprint)
    if isinstance(value, ProcessRecord):
        return dict(program=value.program, returncode=value.returncode, exception_chain=_exception_chain(value.exception))
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return {f.name: _serializable(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, (tuple, list)):
        return [_serializable(v) for v in value]
    if isinstance(value, dict):
        return {k: _serializable(v) for k, v in value.items()}
    return value


@dataclass(frozen=True)
class OptimizationResult:
    status: Outcome
    diagnostics: tuple[FailureDiagnostic, ...]
    gates: tuple[GateEvidence, ...]
    attempts: tuple[OptimizationAttempt, ...]
    accepted_frame: GeometryFrame | None = None
    observed_frame: GeometryFrame | None = None
    exception: BaseException | None = None

    @property
    def converged(self):
        return self.status in {Outcome.CONVERGED, Outcome.RECOVERED_CONVERGED}

    @property
    def workflow_accepted(self):
        return any(g.axis == 'workflow' and g.state is GateState.PASSED for g in self.gates)

    def to_dict(self):
        result = {f.name: _serializable(getattr(self, f.name)) for f in fields(self) if f.name != 'exception'}
        chain = _exception_chain(self.exception)
        result.update(schema_version=1, reason=('VALIDATED_NUMERICAL_CONVERGENCE' if self.converged
                      else self.diagnostics[0].code.value if self.diagnostics else 'UNKNOWN'),
                      exception_chain=chain, converged=self.converged,
                      workflow_accepted=self.workflow_accepted)
        return result

    def to_json(self):
        return canonical_json_bytes(self.to_dict()).decode()


def legacy_diagnostic(status):
    """Translate status evidence; a legacy success is never a certificate."""
    if status is True or status == 'converged':
        code, message = FailureCode.ARTIFACT_INCOMPLETE, 'Legacy success requires independent validation'
    elif status in ('CycleExceeded', 'cycle_exceeded', 'cycleexceeded'):
        code, message = FailureCode.OPT_MAX_STEPS, 'Optimizer reported a cycle limit'
    elif status in ('SCFFailed', 'scf_failed'):
        code, message = FailureCode.CALC_SCF_MAX_STEPS, 'Calculator reported SCF failure'
    else:
        code, message = FailureCode.UNKNOWN, f'Unclassified legacy status: {status!r}'
    return FailureDiagnostic(code, 'legacy_adapter', message)
