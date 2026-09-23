"""Independent energy/force consistency tests for the AFIR bias."""

import numpy as np

from pyar.biases.afir import isotropic
from pyar.data.units import angstrom2bohr


def test_afir_bias_force_matches_central_energy_difference():
    fragments = [[0, 2], [1, 3]]
    symbols = ["C", "H", "O", "H"]
    coordinates = angstrom2bohr(np.asarray([
        [-1.2, 0.1, 0.0],
        [1.1, 0.0, 0.2],
        [-1.1, 1.0, 0.3],
        [1.3, 1.2, -0.1],
    ]))
    gamma = 125.0
    energy, forces = isotropic(fragments, symbols, coordinates, gamma)
    displacement = 2.0e-5
    numerical_forces = np.zeros_like(coordinates)
    for atom in range(len(symbols)):
        for axis in range(3):
            forward = coordinates.copy()
            backward = coordinates.copy()
            forward[atom, axis] += displacement
            backward[atom, axis] -= displacement
            forward_energy, _ = isotropic(fragments, symbols, forward, gamma)
            backward_energy, _ = isotropic(fragments, symbols, backward, gamma)
            numerical_forces[atom, axis] = -(forward_energy - backward_energy) / (2 * displacement)

    assert np.isfinite(energy)
    np.testing.assert_allclose(forces, numerical_forces, rtol=2e-7, atol=2e-9)
