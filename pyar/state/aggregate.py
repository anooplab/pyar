"""Versioned, inspectable restart state for aggregation workflows."""

from __future__ import annotations


import numpy as np


STATE_VERSION = 1
STATE_FILENAME = "state.json"


class AggregateStateError(RuntimeError):
    """Raised when an aggregation run cannot be safely resumed."""


def _json_value(value):
    """Return ``value`` in a JSON-serializable representation."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)
