"""Compatibility shim for the legacy ``pyar.orientation_sampling`` path.

Use :mod:`pyar.sampling.sphere`, :mod:`pyar.sampling.rotation`, and
:mod:`pyar.sampling.metrics` instead.
"""

from pyar.sampling.metrics import (
    RotationCoverageMetrics,
    SphereCoverageMetrics,
    quaternion_coverage_metrics,
    sphere_coverage_metrics,
)
from pyar.sampling.rotation import (
    canonicalize_quaternion,
    halton_quaternions,
    quaternions_to_euler_zxz,
    random_quaternions,
)
from pyar.sampling.sphere import (
    fibonacci_sphere,
    generate_directions,
    latin_hypercube_sphere,
    maximin_select,
    random_sphere,
)

__all__ = [
    "RotationCoverageMetrics",
    "SphereCoverageMetrics",
    "canonicalize_quaternion",
    "fibonacci_sphere",
    "generate_directions",
    "halton_quaternions",
    "latin_hypercube_sphere",
    "maximin_select",
    "quaternion_coverage_metrics",
    "quaternions_to_euler_zxz",
    "random_quaternions",
    "random_sphere",
    "sphere_coverage_metrics",
]
