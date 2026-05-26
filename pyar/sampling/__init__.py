"""Sampling entrypoints for the target package layout."""

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
from pyar.sampling.trial_generator import (
    create_composite_molecule,
    create_trial_geometries,
    generate_points,
    generate_trial_vectors,
    merge_two_molecules,
    write_trial_vectors,
)

__all__ = [
    "RotationCoverageMetrics",
    "SphereCoverageMetrics",
    "canonicalize_quaternion",
    "fibonacci_sphere",
    "generate_directions",
    "generate_trial_vectors",
    "halton_quaternions",
    "latin_hypercube_sphere",
    "maximin_select",
    "quaternion_coverage_metrics",
    "quaternions_to_euler_zxz",
    "random_quaternions",
    "random_sphere",
    "sphere_coverage_metrics",
    "create_composite_molecule",
    "create_trial_geometries",
    "generate_points",
    "merge_two_molecules",
    "write_trial_vectors",
]
