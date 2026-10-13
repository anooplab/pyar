"""Installation paths must survive workflow directories and preserve XYZ input."""

from pathlib import Path
import runpy
import sys
from types import SimpleNamespace

import pytest


def test_relative_model_path_survives_workflow_directory_change(tmp_path, monkeypatch):
    from pyar.backends import aimnet2_assets

    module_path = Path(aimnet2_assets.__file__)
    model = tmp_path / 'model.jpt'
    model.write_text('model placeholder')
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('PYAR_AIMNET2_MODEL', 'model.jpt')
    assets = runpy.run_path(str(module_path))
    job = tmp_path / 'job'
    job.mkdir()
    monkeypatch.chdir(job)
    assert assets['MODEL_PATH'] == str(model)
    assets['validate_aimnet2_runtime_assets'](include_script=False)


def test_explicit_charge_xyz_optimizer_preserves_all_frames(tmp_path, monkeypatch):
    pytest.importorskip('ase')
    from ase.calculators.calculator import Calculator
    from ase.io import read
    from pyar.AIMNet2.calculators import aimnet2_ase_opt

    source, output = tmp_path / 'input.xyz', tmp_path / 'output.xyz'
    source.write_text('1\nfirst\nH 0 0 0\n1\nsecond\nHe 0 0 2\n')
    optimized, charges = [], []

    def calculator(model, charge):
        charges.append(charge)
        return Calculator()

    def optimize(atoms, **kwargs):
        optimized.append(atoms.get_chemical_formula())
        atoms.positions[:, 2] += 1

    torch = SimpleNamespace(
        backends=SimpleNamespace(cuda=SimpleNamespace(matmul=SimpleNamespace()),
                                 cudnn=SimpleNamespace()),
        cuda=SimpleNamespace(is_available=lambda: False),
        device=lambda name: name,
        jit=SimpleNamespace(load=lambda *args, **kwargs: object()),
    )
    monkeypatch.setitem(sys.modules, 'torch', torch)
    monkeypatch.setitem(sys.modules,
                        'pyar.AIMNet2.calculators.aimnet2ase',
                        SimpleNamespace(AIMNet2Calculator=calculator))
    monkeypatch.setattr(aimnet2_ase_opt, 'optimize', optimize)
    monkeypatch.setattr('sys.argv', ['pyar-aimnet2-ase-opt', '--charge', '1',
                                    'model.jpt', str(source), str(output)])
    aimnet2_ase_opt.main()
    frames = read(output, index=':')
    assert optimized == ['H', 'He']
    assert charges == [1, 1]
    assert [atoms.get_chemical_formula() for atoms in frames] == ['H', 'He']
    assert [atoms.positions[0, 2] for atoms in frames] == [1, 3]
