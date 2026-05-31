"""Compatibility shim for the legacy ``pyar.trial_generation`` path.

Use :mod:`pyar.sampling.trial_generator` instead.
"""

from pyar.sampling.trial_generator import (
    broken,
    check_close_contact,
    create_composite_molecule,
    create_trial_geometries,
    ellipsoid_wall_potential,
    generate_points,
    generate_trial_vectors,
    get_connectivity,
    merge_two_molecules,
    minimum_separation,
    plot_points,
    polar_to_cartesian,
    spherical_to_cartesian,
    write_trial_vectors,
)

__all__ = [
    "broken",
    "check_close_contact",
    "create_composite_molecule",
    "create_trial_geometries",
    "ellipsoid_wall_potential",
    "generate_points",
    "generate_trial_vectors",
    "get_connectivity",
    "merge_two_molecules",
    "minimum_separation",
    "plot_points",
    "polar_to_cartesian",
    "spherical_to_cartesian",
    "write_trial_vectors",
]
