"""Failures that prevent a selected chemistry backend from completing a job."""


class BackendExecutionError(RuntimeError):
    """A backend failed or produced no usable result for a requested job."""
