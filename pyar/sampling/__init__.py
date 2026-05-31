"""Public sampling entrypoints for PyAR 2.0."""

from pyar.sampling.metrics import (
    RotationCoverageMetrics,
    SphereCoverageMetrics,
    quaternion_coverage_metrics,
    sphere_coverage_metrics,
)
from pyar.sampling.rotation import halton_quaternions, random_quaternions
from pyar.sampling.sphere import (
    fibonacci_sphere,
    generate_directions,
    latin_hypercube_sphere,
    maximin_select,
    random_sphere,
)
from pyar.sampling.trial_generator import generate_trial_vectors

__all__ = [
    "RotationCoverageMetrics",
    "SphereCoverageMetrics",
    "fibonacci_sphere",
    "generate_directions",
    "generate_trial_vectors",
    "halton_quaternions",
    "latin_hypercube_sphere",
    "maximin_select",
    "quaternion_coverage_metrics",
    "random_quaternions",
    "random_sphere",
    "sphere_coverage_metrics",
]
