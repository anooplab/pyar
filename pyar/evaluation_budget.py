"""Persistent physical-backend evaluation budgets for optimizer subprocesses."""

from __future__ import annotations

import json
import os
from pathlib import Path
import time


class EvaluationBudgetExhausted(RuntimeError):
    """A geometry optimizer reached its configured backend evaluation limit."""


def reserve_backend_evaluation(path):
    """Durably reserve one backend call or raise before starting an over-budget call.

    The microsolvation workflow runs one geometry optimizer at a time under its
    run lock. Atomic replacement protects the counter from torn writes if the
    calculator subprocess exits unexpectedly.
    """
    counter_path = Path(path).expanduser().resolve()
    data = json.loads(counter_path.read_text(encoding="utf-8"))
    if time.time() >= float(data["deadline_epoch"]):
        raise EvaluationBudgetExhausted("microsolvation resampling time limit reached")
    used = int(data["used"])
    limit = int(data["limit"])
    if used >= limit:
        raise EvaluationBudgetExhausted("microsolvation resampling evaluation budget exhausted")
    data["used"] = used + 1
    temporary = counter_path.with_name(f"{counter_path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(data, stream, sort_keys=True)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, counter_path)
    descriptor = os.open(counter_path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return data["used"]
