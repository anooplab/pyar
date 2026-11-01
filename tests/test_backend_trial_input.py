"""Trial geometry must belong to the current attempt."""
import numpy as np
from pyar.backends import SF
from pyar.core.molecule import Molecule
from pyar.io.xyz import parse_xyz


def test_repeated_attempt_rewrites_same_trial_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = Molecule(['H', 'H'], [[0, 0, 0], [0, 0, .8]], name='retry')
    SF(molecule)
    molecule.coordinates[1, 2] = 1.2
    SF(molecule)
    _, coordinates, *_ = parse_xyz('trial_retry.xyz')
    np.testing.assert_allclose(coordinates, molecule.coordinates, atol=1e-5)


