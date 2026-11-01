"""Electronic state and attempt geometry at the optimizer/backend boundary."""
from pathlib import Path
from types import SimpleNamespace
import shutil

import numpy as np
import pytest

from pyar import optimiser
from pyar.core.molecule import Molecule
from pyar.io.xyz import parse_xyz


def test_optimise_uses_molecule_state_without_mutating_parameters(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = Molecule(['H', 'H'], [[0, 0, 0], [0, 0, .8]], name='ion', charge=1, multiplicity=2)
    supplied = dict(software='xtb', charge=0, multiplicity=1, xtb_unpaired_electrons=0)
    observed = {}
    monkeypatch.setattr(optimiser, '_backend_cache_provenance', lambda params: None)

    def build(molecule, parameters):
        observed.update(parameters)
        return SimpleNamespace(optimize=lambda: False)

    monkeypatch.setattr(optimiser, 'build_geometry', build)
    from pyar.backend_errors import BackendExecutionError
    with pytest.raises(BackendExecutionError, match="optimization failed"):
        optimiser.optimise(molecule, supplied)
    assert observed['charge'] == 1
    assert observed['multiplicity'] == 2
    assert observed['xtb_unpaired_electrons'] == 1
    assert supplied == dict(software='xtb', charge=0, multiplicity=1, xtb_unpaired_electrons=0)


@pytest.mark.external_backend
@pytest.mark.skipif(shutil.which('xtb') is None, reason='native xTB executable required')
def test_charged_aggregate_native_xtb(tmp_path, monkeypatch):
    """A sodium/water aggregate must run with its +1 state, without caller injection."""
    from pyar.workflows.aggregate import aggregate
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    monkeypatch.setenv('MKL_NUM_THREADS', '1')
    sodium = Molecule(['Na'], [[0, 0, 0]], name='sodium', charge=1, multiplicity=1)
    water = Molecule(['O', 'H', 'H'], [[0, 0, 0], [.758, 0, .504], [-.758, 0, .504]], name='water')
    result = aggregate([sodium, water], [1, 1], 1,
                       dict(software='xtb', xtb_model='gfn2', opt_threshold='normal',
                            nprocs=1, gamma=None), 1, 0, 1, None,
                       connectivity_policy='off')
    outputs = list(Path('aggregates').rglob('xtb.out'))
    assert outputs, result
    texts = [p.read_text() for p in outputs]
    assert all('-chrg 1' in text for text in texts)
    assert any('GEOMETRY OPTIMIZATION CONVERGED' in text for text in texts)
    frames = list(Path('aggregates').rglob('result_*.xyz'))
    assert frames
    for path in frames:
        symbols, coordinates, _, _, energy = parse_xyz(path)
        assert sorted(symbols) == ['H', 'H', 'Na', 'O']
        assert np.isfinite(coordinates).all()
        assert energy is not None and np.isfinite(energy)
