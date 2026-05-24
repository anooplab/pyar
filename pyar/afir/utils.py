"""Shared AFIR helper functions."""

import logging
import math

afir_logger = logging.getLogger("pyar.afir")


def resolve_gamma(value, fallback=100.0):
    """Return a numeric AFIR gamma value.

    The current workflow keeps a backward-compatible fallback so older call
    sites that do not supply gamma still behave as before.
    """
    if value is None:
        return float(fallback)
    try:
        gamma = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid AFIR gamma value: {value!r}") from None
    if not math.isfinite(gamma):
        raise ValueError(f"Invalid AFIR gamma value: {value!r}")
    return gamma
