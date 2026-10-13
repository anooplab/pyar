"""Run from outside the checkout after a normal (non-editable) installation."""

import argparse
import importlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def run(*arguments, expected=0, env=None, cwd=None):
    result = subprocess.run(arguments, text=True, capture_output=True, env=env, cwd=cwd)
    assert result.returncode == expected, (arguments, result.stdout, result.stderr)
    assert 'Traceback' not in result.stderr, result.stderr
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['core', 'aggregate', 'reaction', 'aimnet2', 'neb', 'ml', 'analysis', 'microsolvate', 'grow', 'solvate'])
    mode = parser.parse_args().mode
    import pyar
    from pyar.dependencies import require_workflow

    assert 'site-packages' in pyar.__file__, pyar.__file__
    run(sys.executable, '-m', 'pyar.modern_cli', '--help')
    run(sys.executable, '-m', 'pyar.cli', '--help')
    doctor = run(sys.executable, '-m', 'pyar.modern_cli', 'doctor', '--json')
    assert json.loads(doctor.stdout)['ok']
    if mode == 'core':
        for module in ('ase', 'torch', 'rdkit', 'dscribe', 'geometric', 'sella'):
            assert importlib.util.find_spec(module) is None, module
        missing = run(sys.executable, '-m', 'pyar.modern_cli', 'aggregate', '--formula', 'Au2', '--check', expected=2)
        assert 'pyar-chem[aggregate]' in missing.stderr
        print('Core-only imports, CLI help and doctor passed.')
        return
    if mode == 'ml':
        import mlatom
        import torchani
        from pyar.backends import ani
        from pyar import mlatom as compatibility
        assert 'site-packages/mlatom' in mlatom.__file__
        assert compatibility.data is mlatom.data
        print('MLatom and TorchANI runtime imports passed.')
        return
    if mode == 'analysis':
        import MDAnalysis
        import pandas
        import matplotlib
        from pyar.scripts.descriptor import _require_descriptor_dependencies
        _require_descriptor_dependencies()
        print('Analysis capability imports passed.')
        return
    workflow = 'conformer' if mode == 'aimnet2' else mode
    require_workflow(workflow)
    modules = {
        'microsolvate': 'microsolvation', 'solvate': 'microsolvation',
        'grow': 'grow', 'reaction': 'reaction', 'aimnet2': 'conformer',
        'aggregate': 'aggregate', 'neb': 'scan_path',
    }
    module_name = modules[mode]
    importlib.import_module('pyar.workflows.' + module_name)
    if mode == 'aggregate':
        import os
        from pyar import Molecule
        from pyar.selection.features import compute_feature_matrix
        from pyar.selection.clusterers import cluster_molecules

        molecules = [Molecule(['H', 'H'], [[0, 0, 0], [0, 0, distance]], energy=-1)
                     for distance in (0.7, 0.8, 0.9, 1.0)]
        features = compute_feature_matrix(molecules, feature='soap')
        assert features.name == 'soap', features
        clustered = cluster_molecules(molecules, feature='soap', algorithm='hdbscan')
        assert clustered is not None
        assert importlib.util.find_spec('torch') is None
        assert importlib.util.find_spec('geometric') is None
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'h2.xyz'
            source.write_text('2\nH2\nH 0 0 0\nH 0 0 0.74\n')
            for check in ([], ['--check']):
                result = run(sys.executable, '-m', 'pyar.modern_cli', 'aggregate',
                             str(source), '--size', '2', '--backend', 'xtb', *check,
                             expected=2, env=dict(os.environ, PATH=''), cwd=directory)
                assert "Missing executable 'xtb'" in result.stderr
                assert 'mamba install -c conda-forge xtb' in result.stderr
            assert list(Path(directory).iterdir()) == [source]
    elif mode == 'aimnet2':
        from pyar.AIMNet2.calculators.aimnet2ase import AIMNet2Calculator
        from pyar.backends.aimnet_2 import Aimnet2
        import rdkit
        import torch
        assert AIMNet2Calculator and Aimnet2 and rdkit and torch
        assert importlib.util.find_spec('openbabel') is None
        assert importlib.util.find_spec('geometric') is None
        run(sys.executable, '-m', 'pyar.modern_cli', 'conformer', 'CCO',
                '--input-format', 'smiles', '--check')
    elif mode == 'reaction':
        import os
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / 'a.xyz', Path(directory) / 'b.xyz'
            first.write_text('2\nH2\nH 0 0 0\nH 0 0 0.74\n')
            second.write_text(first.read_text())
            env = dict(os.environ, PATH='')
            for check in ([], ['--check']):
                result = run(sys.executable, '-m', 'pyar.modern_cli', 'react',
                             str(first), str(second), '--backend', 'orca',
                             '--method', 'BP86', '--basis', 'def2-SVP',
                             '--bias-max', '100', *check, expected=2, env=env,
                             cwd=directory)
                assert "Missing executable 'orca'" in result.stderr, result.stderr
                assert 'https://www.faccts.de/' in result.stderr
            assert sorted(path.name for path in Path(directory).iterdir()) == ['a.xyz', 'b.xyz']
    elif mode == 'neb':
        importlib.import_module('pyar.neb')
        run(sys.executable, '-m', 'pyar.neb', '--help')
    print(f'{mode}: installation smoke passed.')


if __name__ == '__main__':
    main()
