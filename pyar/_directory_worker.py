"""Private one-shot worker; its working directory is set by subprocess creation."""

import importlib
import logging
from pathlib import Path
import pickle
import random
import sys
import traceback

import numpy as np
from pyar.backend_errors import BackendExecutionError, BackendOutputError
from pyar.recovery.observations import capture_processes


def main():
    source, target = map(Path, sys.argv[1:])
    module, name, args, kwargs, python_state, numpy_state = pickle.loads(
        source.read_bytes()
    )
    operation = importlib.import_module(module)
    for component in name.split("."):
        operation = getattr(operation, component)
    random.setstate(python_state)
    np.random.set_state(numpy_state)
    value, failure, remote_traceback = None, None, ""
    logs = []

    class CaptureLog(logging.Handler):
        def emit(self, record):
            logs.append((record.name, record.levelno, self.format(record)))

    logging.getLogger("pyar").addHandler(CaptureLog())
    logging.getLogger("pyar").setLevel(logging.DEBUG)
    with capture_processes() as processes:
        try:
            value = operation(*args, **kwargs)
        except BaseException as exc:
            failure = (
                BackendExecutionError(str(exc)) if isinstance(exc, SystemExit) else exc
            )
            remote_traceback = traceback.format_exc()
    response = dict(
        result=value,
        exception=failure,
        processes=processes,
        logs=logs,
        traceback=remote_traceback,
        cause=failure.__cause__ if failure is not None else None,
        random_state=random.getstate(),
        numpy_state=np.random.get_state(),
    )
    try:
        encoded = pickle.dumps(response, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as exc:
        response.update(
            result=None,
            cause=None,
            exception=BackendOutputError(
                f"Could not transfer directory operation result: {exc}"
            ),
        )
        encoded = pickle.dumps(response, protocol=pickle.HIGHEST_PROTOCOL)
    target.write_bytes(encoded)


if __name__ == "__main__":
    main()
