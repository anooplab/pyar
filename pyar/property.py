# -*- coding: utf-8 -*-
"""
Calculate or retrieve properties of the input molecules.


"""
from pyar.constants import BOND_DISTANCE_FACTOR
from math import pi

import numpy as np
from scipy.spatial import distance as scipy_distance


def distance(coords_a, coords_b):
    """
    Calculates distance between atoms a and b.

    :type coords_a: ndarray(float)
    :param coords_a: Coordinate of atom a
    :type coords_b: ndarray(float)
    :param coords_b: Coordinate of atom b
    :return: The distance between atom a and b
    :rtype: float
    """
    return np.linalg.norm(coords_a - coords_b)


def get_centroid(coordinates):
    return coordinates.mean(0)


def get_centre_of_mass(coordinates, atomic_mass):
    return np.average(coordinates, axis=0, weights=atomic_mass)


def get_average_radius(coordinates, centroid):
    return np.mean(np.array([distance(centroid, coord_i)
                             for coord_i in coordinates]))


def get_std_of_radius(coordinates, centroid):
    return np.std([distance(centroid, coord_i)
                   for coord_i in coordinates])


def get_principal_axes(moments_of_inertia_tensor):
    return np.linalg.eigvals(moments_of_inertia_tensor)


def get_distance_list(coordinates):
    return scipy_distance.pdist(coordinates)


def get_distance_matrix(coordinates):
    dm = np.zeros((len(coordinates), len(coordinates)))
    for i, c in enumerate(coordinates):
        for j, d in enumerate(coordinates):
            dm[i, j] = np.linalg.norm(c - d)
    return dm


def get_bond_matrix(coordinates, covalent_radius):
    """return bond matrix"""
    dm = get_distance_matrix(coordinates)
    bm = np.zeros((len(coordinates), len(coordinates)), dtype=int)
    for i, c in enumerate(coordinates):
        for j, d in enumerate(coordinates):
            sum_of_covalent_radii = (covalent_radius[i] + covalent_radius[j]) * BOND_DISTANCE_FACTOR
            if i != j and dm[i, j] < sum_of_covalent_radii:
                bm[i, j] = 1
    return bm


def get_connectivity(coordinates, covalent_radius):
    """return connection graph"""
    import collections
    dm = get_distance_matrix(coordinates)
    bond_graph = collections.defaultdict(list)
    for i, c in enumerate(coordinates):
        for j, d in enumerate(coordinates):
            sum_of_covalent_radii = (covalent_radius[i] + covalent_radius[j]) * BOND_DISTANCE_FACTOR
            if i != j and dm[i, j] < sum_of_covalent_radii:
                bond_graph[j].append(i)
    return bond_graph


def calculate_angle(a1, b1, c1):
    """Return the A-B-C angle in degrees; reject a zero-length angle ray."""
    v1 = a1 - b1
    v2 = c1 - b1
    norm1, norm2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if norm1 == 0.0 or norm2 == 0.0:
        raise ValueError("Bond angle is undefined for a zero-length angle vector")
    cosine = np.dot(v1, v2) / (norm1 * norm2)
    # Roundoff can put a normalized dot product just outside arccos' domain.
    return np.arccos(np.clip(cosine, -1.0, 1.0)) * 180 / pi


def hydrogen_bond_analysis(coordinates, covalent_radius, atomic_number, atoms_list):
    """Return directed H-to-acceptor distances using the existing contact cutoffs."""
    dm = get_distance_matrix(coordinates)
    bm = get_bond_matrix(coordinates, covalent_radius)
    hbm = np.zeros((len(coordinates), len(coordinates)), dtype=float)
    for i, c in enumerate(coordinates):
        for j, d in enumerate(coordinates):
            if i != j and atomic_number[i] == 1 and 1.5 < dm[i, j] < 2.5:
                bnd = int(np.where(bm[i, :])[0].item())
                angle = calculate_angle(coordinates[bnd], coordinates[i],
                                        coordinates[j])
                if angle > 160.0:
                    print(f"{atoms_list[i]}({i + 1}) - {atoms_list[j]}({j + 1}) = {dm[i, j]}")
                    hbm[i, j] = dm[i, j]
    return hbm


def calculate_dihedral(p0, p1, p2, p3):
    """Return the existing signed torsion in degrees; reject undefined planes."""
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    norm = np.linalg.norm(b1)
    if norm == 0.0:
        raise ValueError("Torsion is undefined for a zero-length central torsion bond")
    if np.linalg.norm(np.cross(b0, b1)) == 0.0 or np.linalg.norm(np.cross(b1, b2)) == 0.0:
        raise ValueError("Undefined torsion: collinear geometry does not define both planes")
    b1 /= norm

    # vector projections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(v, b1), w)
    return np.degrees(np.arctan2(y, x))
