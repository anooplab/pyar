"""Isolate relative-file third-party operations using subprocess ``cwd``.

The private, one-shot IPC files live in a mode-0700 temporary directory, never
in workflow authority/cache records. They are generated and consumed by PyAR
within one call, not a public input, restart format or reusable checkpoint.
"""

import logging
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys
import tempfile

import numpy as np

from pyar.backend_errors import BackendExecutionError


def run_in_directory(operation, *args, directory, **kwargs):
    """Execute an importable PyAR operation without changing the caller's cwd."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if (
        not operation.__module__.startswith("pyar.")
        or "<locals>" in operation.__qualname__
    ):
        raise ValueError("Directory operations must be importable PyAR functions")
    request = (
        operation.__module__,
        operation.__qualname__,
        args,
        kwargs,
        random.getstate(),
        np.random.get_state(),
    )
    with tempfile.TemporaryDirectory(prefix="pyar-directory-ipc-") as temporary:
        source, target = Path(temporary) / "request", Path(temporary) / "response"
        source.write_bytes(pickle.dumps(request, protocol=pickle.HIGHEST_PROTOCOL))
        environment = worker_environment()
        result = subprocess.run(
            [sys.executable, "-m", "pyar._directory_worker", str(source), str(target)],
            cwd=directory,
            env=environment,
            capture_output=True,
            check=False,
        )
        if result.returncode or not target.is_file():
            raise BackendExecutionError(
                f"Directory operation failed in {directory}: "
                + result.stderr.decode(errors="replace")
            )
        if result.stdout:
            sys.stdout.write(result.stdout.decode(errors="replace"))
        if result.stderr:
            sys.stderr.write(result.stderr.decode(errors="replace"))
        try:
            response = pickle.loads(target.read_bytes())
        except Exception as exc:
            raise BackendExecutionError(
                "Could not read directory operation response"
            ) from exc
    for name, level, message in response["logs"]:
        logger = logging.getLogger(name)
        if logger.isEnabledFor(level):
            logger.log(level, message)
    random.setstate(response["random_state"])
    np.random.set_state(response["numpy_state"])
    from pyar.recovery.observations import extend_process_records

    extend_process_records(response["processes"])
    if response["exception"] is not None:
        failure = response["exception"]
        failure.add_note(response["traceback"])
        failure.__cause__ = response["cause"]
        raise failure
    return response["result"]


def worker_environment():
    """Use this installed/check-out package and resolve externally owned assets."""
    environment = dict(os.environ)
    for key in ("PYAR_AIMNET2_MODEL", "PYAR_ORCA_AIQM1_PROGEXT"):
        if environment.get(key):
            environment[key] = str(Path(environment[key]).expanduser().resolve())
    package_parent = str(Path(__file__).resolve().parents[1])
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, (package_parent, environment.get("PYTHONPATH")))
    )
    return environment
