"""Shared subprocess helpers for PyAR backend wrappers."""

from contextlib import ExitStack
import logging
import os
import subprocess as subp

logger = logging.getLogger("pyar.backends.subprocess")


def _normalize_command(command):
    """Return a subprocess-safe command list with string arguments."""
    if isinstance(command, (str, bytes, os.PathLike)):
        return [os.fspath(command)]
    return [
        os.fspath(part) if isinstance(part, os.PathLike) else str(part)
        for part in command
    ]


def _open_shared(stack, cache, path, mode):
    if path is None:
        return None
    if path not in cache:
        cache[path] = stack.enter_context(open(path, mode))
    return cache[path]


def run_command(command, stdout_path=None, stderr_path=None, stdin_path=None):
    """Run a command and return its exit status.

    When stdout and stderr should go to the same file, pass the same path for
    both arguments.
    """
    command = _normalize_command(command)
    cwd = os.getcwd()
    logger.debug(
        "running command=%s cwd=%s stdout=%s stderr=%s stdin=%s",
        command, cwd, stdout_path, stderr_path, stdin_path,
    )
    with ExitStack() as stack:
        files = {}
        stdout = _open_shared(stack, files, stdout_path, "w")
        stderr = _open_shared(stack, files, stderr_path, "w")
        stdin = _open_shared(stack, files, stdin_path, "r")
        proc = subp.Popen(command, stdout=stdout, stderr=stderr, stdin=stdin)
        proc.communicate()
        if proc.returncode != 0:
            logger.error(
                "command failed command=%s returncode=%s cwd=%s stdout=%s stderr=%s stdin=%s",
                command, proc.returncode, cwd, stdout_path, stderr_path, stdin_path,
            )
        else:
            logger.debug("command finished command=%s returncode=0 cwd=%s", command, cwd)
        return proc.returncode


def run_output(command, stderr_path=None):
    """Run a command and return its stdout bytes.

    Raises:
        subprocess.CalledProcessError: If the command exits with a non-zero
            status.
    """
    command = _normalize_command(command)
    cwd = os.getcwd()
    logger.debug("running command for output=%s cwd=%s stderr=%s", command, cwd, stderr_path)
    with ExitStack() as stack:
        files = {}
        stderr = _open_shared(stack, files, stderr_path, "w")
        proc = subp.Popen(command, stdout=subp.PIPE, stderr=stderr)
        stdout, _ = proc.communicate()
        if proc.returncode != 0:
            logger.error(
                "output command failed command=%s returncode=%s cwd=%s stderr=%s",
                command, proc.returncode, cwd, stderr_path,
            )
            raise subp.CalledProcessError(proc.returncode, command, output=stdout)
        logger.debug("output command finished command=%s returncode=0 cwd=%s", command, cwd)
        return stdout
