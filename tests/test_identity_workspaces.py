"""Identity conversion owns diagnostic paths without changing chemistry gates."""
from pathlib import Path
import shutil

import pytest

from pyar.reaction_identity import molecule_identity_from_xyz


@pytest.mark.external_backend
@pytest.mark.skipif(shutil.which('obabel') is None, reason='OpenBabel executable required')
def test_native_identity_log_stays_in_explicit_workspace(tmp_path, monkeypatch):
    import os
    root = tmp_path/'identity'; root.mkdir()
    (root/'water.xyz').write_text('3\nwater\nO 0 0 0\nH .9572 0 0\nH -.239 .927 0\n')
    caller = Path.cwd()
    def forbidden(*args):
        pytest.fail('identity changed process cwd')
    monkeypatch.setattr(os, 'chdir', forbidden)
    result = molecule_identity_from_xyz('water.xyz', directory=root)
    assert result['inchi'].startswith('InChI=1S/H2O/')
    assert result['smiles'] == 'O'
    assert (root/'OBabel.log').is_file()
    assert Path.cwd() == caller
