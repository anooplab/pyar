"""Geometry helpers for :mod:`pyar.Molecule`."""

from __future__ import annotations

from math import cos, sin

import numpy as np

import pyar.property

__all__ = [
    "rotation_matrix",
    "moments_of_inertia_tensor",
    "translate",
    "move_to_origin",
    "move_to_centre_of_mass",
    "rotate_3d",
    "align",
]


def _require_coordinates(molecule) -> np.ndarray:
    if molecule.coordinates is None:
        raise ValueError("molecule has no coordinates")
    return molecule.coordinates


def rotation_matrix(angles) -> np.ndarray:
    """Build the Z-X-Z' Euler rotation matrix."""
    phi, theta, psi = angles
    matrix_d = np.array(((cos(phi), sin(phi), 0.0), (-sin(phi), cos(phi), 0.0), (0.0, 0.0, 1.0)))
    matrix_c = np.array(((1.0, 0.0, 0.0), (0.0, cos(theta), sin(theta)), (0.0, -sin(theta), cos(theta))))
    matrix_b = np.array(((cos(psi), sin(psi), 0.0), (-sin(psi), cos(psi), 0.0), (0.0, 0.0, 1.0)))
    return np.dot(matrix_b, np.dot(matrix_c, matrix_d))


def moments_of_inertia_tensor(molecule):
    """Return the inertia tensor without mutating the molecule."""
    coordinates = _require_coordinates(molecule)
    mass = np.asarray(molecule.atomic_mass, dtype=float)
    centered = coordinates - pyar.property.get_centre_of_mass(coordinates, mass)
    x = centered[:, 0]
    y = centered[:, 1]
    z = centered[:, 2]
    i_xx = np.sum(mass * (y * y + z * z))
    i_yy = np.sum(mass * (x * x + z * z))
    i_zz = np.sum(mass * (x * x + y * y))
    i_xy = -np.sum(mass * (x * y))
    i_yz = -np.sum(mass * (y * z))
    i_xz = -np.sum(mass * (x * z))
    return np.array([[i_xx, i_xy, i_xz], [i_xy, i_yy, i_yz], [i_xz, i_yz, i_zz]]) / molecule.number_of_atoms


def translate(molecule, magnitude):
    """Translate the molecule in place by ``magnitude``."""
    molecule.coordinates = _require_coordinates(molecule) - magnitude
    return molecule


def move_to_origin(molecule):
    """Shift the molecule so its centroid lies at the origin."""
    return translate(molecule, pyar.property.get_centroid(_require_coordinates(molecule)))


def move_to_centre_of_mass(molecule):
    """Shift the molecule so its center of mass lies at the origin."""
    return translate(
        molecule,
        pyar.property.get_centre_of_mass(_require_coordinates(molecule), molecule.atomic_mass),
    )


def rotate_3d(molecule, angles):
    """Rotate the molecule in place using Z-X-Z' Euler angles."""
    matrix_a = rotation_matrix(angles)
    coordinates = _require_coordinates(molecule)
    centroid = molecule.centroid
    molecule.coordinates = np.dot(coordinates - centroid, matrix_a.T) + centroid
    return molecule


def align(molecule):
    """Align the molecule with its principal axes in place."""
    move_to_centre_of_mass(molecule)
    moi = moments_of_inertia_tensor(molecule)
    _, eigen_vectors = np.linalg.eig(moi)
    transformed_coordinates = np.array(np.dot(molecule.coordinates, eigen_vectors))
    order = [0, 1, 2]
    for p in range(3):
        for q in range(p + 1, 3):
            if moi.item(p, p) < moi.item(q, q):
                order[p], order[q] = order[q], order[p]
    move_axes = np.zeros((3, 3))
    for p in range(3):
        move_axes[p][order[p]] = 1.0
    molecule.coordinates = np.dot(transformed_coordinates, move_axes)
    return molecule
