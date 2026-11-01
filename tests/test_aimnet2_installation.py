"""Installation paths must survive workflow directories and preserve XYZ input."""

from pathlib import Path
import runpy
import sys
from types import SimpleNamespace

import pytest


@pytest.mark.parametrize('configured', [None, ''])
def test_model_environment_is_required_even_with_local_model(tmp_path, monkeypatch, configured):
    from pyar.backends import aimnet2_assets

    # Neither a working-directory model nor the former package location is authority.
    (tmp_path / 'aimnet2_wb97m-d3_0.jpt').write_text('unused placeholder')
    monkeypatch.chdir(tmp_path)
    if configured is None:
        monkeypatch.delenv('PYAR_AIMNET2_MODEL', raising=False)
    else:
        monkeypatch.setenv('PYAR_AIMNET2_MODEL', configured)
    assets = runpy.run_path(aimnet2_assets.__file__)
    assert assets['MODEL_PATH'] == ''
    assert Path(assets['SCRIPT_PATH']).is_file()
    assert Path(assets['SCRIPT_PATH']).parent.name == 'backends'
    assert assets['missing_aimnet2_assets']() == ['PYAR_AIMNET2_MODEL (unset)']
    with pytest.raises(FileNotFoundError, match='Set PYAR_AIMNET2_MODEL'):
        assets['validate_aimnet2_runtime_assets']()


def test_missing_configured_model_is_rejected(tmp_path, monkeypatch):
    from pyar.backends import aimnet2_assets

    missing = tmp_path / 'missing.jpt'
    monkeypatch.setenv('PYAR_AIMNET2_MODEL', str(missing))
    assets = runpy.run_path(aimnet2_assets.__file__)
    assert assets['missing_aimnet2_assets']() == [str(missing)]
    with pytest.raises(FileNotFoundError, match='missing.jpt'):
        assets['validate_aimnet2_runtime_assets']()


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
    assets['validate_aimnet2_runtime_assets']()


def test_explicit_charge_xyz_optimizer_preserves_all_frames(tmp_path, monkeypatch):
    pytest.importorskip('ase')
    from ase.calculators.calculator import Calculator
    from ase.io import read
    from pyar.backends import aimnet2_ase_opt

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
                        'pyar.backends.aimnet2ase',
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
