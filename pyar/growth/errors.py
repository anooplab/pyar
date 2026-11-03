"""Control-flow outcomes owned by aggregation and fixed-addend growth."""


class GrowthStopped(RuntimeError):
    """A stop marker requests interruption before the next growth operation."""
