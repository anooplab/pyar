"""Failures that prevent a selected chemistry backend from completing a job."""


class BackendExecutionError(RuntimeError):
    """A backend failed or produced no usable result for a requested job."""


class BackendSCFError(BackendExecutionError):
    """The electronic structure calculation failed to converge."""


class BackendGradientError(BackendExecutionError):
    """The backend did not produce a valid gradient."""


class BackendGeometryError(BackendExecutionError):
    """Geometry setup or coordinate updating failed."""


class BackendOutputError(BackendExecutionError):
    """Required output is absent, malformed or inconsistent."""


class BackendConvergenceError(BackendExecutionError):
    """The optimization budget ended without convergence.

    usable_geometry indicates only the historical incomplete-frame contract;
    the optimizer and owning chemistry workflow must still validate the frame.
    """
    def __init__(self, message, *, usable_geometry=False):
        super().__init__(message)
        self.usable_geometry = usable_geometry


from enum import Enum


class OptimizationOutcome(str, Enum):
    """Nonfatal workflow observations; serialized values preserve run contracts."""
    CYCLE_EXCEEDED = 'CycleExceeded'
    RELEASE_CANDIDATE = 'ReleaseCandidate'

    def __str__(self):
        return self.value


def backend_operation(operation):
    """Keep unexpected adapter failures inside the typed backend boundary."""
    from functools import wraps

    @wraps(operation)
    def execute(*args, **kwargs):
        try:
            result = operation(*args, **kwargs)
            if result is None or result is False or (isinstance(result, str) and not isinstance(result, OptimizationOutcome)):
                raise BackendOutputError(f'{operation.__qualname__} returned no typed result')
            return result
        except BackendExecutionError:
            raise
        except (Exception, SystemExit) as exc:
            raise BackendExecutionError(f'{operation.__qualname__} failed: {exc}') from exc
    return execute
