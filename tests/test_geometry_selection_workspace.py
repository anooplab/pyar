"""Energy-free selection diagnostics remain inside their owning run."""
from pathlib import Path

from pyar.core.molecule import Molecule
from pyar.growth import service
from pyar.growth.policy import GrowthPolicy


def test_geometry_only_growth_diagnostics_belong_to_explicit_directory(tmp_path):
    caller = Path.cwd()
    before = set(caller.iterdir())
    seed = Molecule(['He'], [[0,0,0]], name='seed')
    monomer = Molecule(['Ne'], [[0,0,0]], name='monomer')
    directory = tmp_path/'addition'
    pool = service.add_one('addition', [seed], monomer, GrowthPolicy(1, {}, 1, connectivity_policy='off'), directory=directory)
    assert pool
    assert (directory/'selection_diagnostics.json').is_file()
    assert (directory/'selected').is_dir()
    assert set(caller.iterdir()) == before
    assert Path.cwd() == caller
