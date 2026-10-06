"""Execution observations shared by modern adapters and the CLI dispatcher.

Adapters publish the request they already resolved; this module never resolves
chemistry or runs workflows. A context-local observer keeps nested invocations
independent and leaves direct Python/legacy callers unaffected.
"""

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class RunObservation:
    prepared: dict | None = None
    result: dict | None = None
    started: bool = False
    on_prepared: object = None
    options: dict = field(default_factory=dict)


_current = ContextVar('pyar_cli_run', default=None)


@contextmanager
def observe_run(on_prepared=None):
    observation = RunObservation(on_prepared=on_prepared)
    token = _current.set(observation)
    try:
        yield observation
    finally:
        _current.reset(token)


def prepared_run(args, *, request, backend, molecules, requirements, outputs,
                 state_lists=True, options=None):
    """Publish validated domain settings after every preflight/restart check."""
    observer = _current.get()
    if observer is None:
        return
    values = dict(vars(args))
    values.update(options or {})
    from pyar.optimization_request import OPTIMIZATION_OPTIONS, unsupported_optimization_options
    for key, value in backend.items():
        target = 'backend' if key == 'software' and 'backend' in values else key
        if target in values:
            if (values[target] is None and key in OPTIMIZATION_OPTIONS
                    and unsupported_optimization_options(backend.get('software'),
                                                         backend.get('geometry_optimizer'), {key})):
                continue
            values[target] = value
    states = []
    for molecule in molecules:
        states.append({'name': molecule.name, 'charge': molecule.charge,
                       'multiplicity': molecule.multiplicity, 'scftype': molecule.scftype})
    for key in ('charge', 'multiplicity', 'scftype'):
        if key in values and states:
            if state_lists:
                values[key] = [state[key] for state in states]
            elif len({state[key] for state in states}) == 1:
                values[key] = states[0][key]
    observer.options = values
    observer.prepared = {
        'request': request, 'backend': dict(backend), 'electronic_states': states,
        'requirements': [{'name': item, 'status': 'satisfied'} for item in requirements],
        'outputs': [str(Path(path).expanduser().resolve()) for path in outputs],
    }
    if observer.on_prepared is not None:
        observer.on_prepared(observer)


def started_run():
    observer = _current.get()
    if observer is not None:
        observer.started = True


def finished_run(result):
    observer = _current.get()
    if observer is not None:
        observer.result = result.to_dict() if hasattr(result, 'to_dict') else result
