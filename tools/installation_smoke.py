"""Run from outside the checkout after a normal (non-editable) installation."""

import argparse
import importlib
import importlib.util
import importlib.metadata
import shutil
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
    parser.add_argument('mode', choices=['core', 'aggregate', 'reaction', 'aimnet2', 'neb', 'ml', 'analysis', 'microsolvate', 'grow', 'solvate', 'conformer'])
    mode = parser.parse_args().mode
    import pyar
    from pyar.dependencies import require_workflow

    assert 'site-packages' in pyar.__file__, pyar.__file__
    assert pyar.__version__ == importlib.metadata.version('pyar-chem')
    for command in ('pyar', 'pyar-cli', 'pyar-react'):
        executable = shutil.which(command)
        assert executable, command
        run(executable, '--help')
    distribution = importlib.metadata.distribution('pyar-chem')
    for entry in distribution.entry_points:
        if entry.group == 'console_scripts':
            assert shutil.which(entry.name), entry.name
    for name in ('pyar.persistence', 'pyar.growth.layout_v2', 'pyar.microsolvation.layout_v2',
                 'pyar.aggregation.layout_v2', 'pyar.state.reaction_v2',
                 'pyar.workflows.scan_layout', 'pyar.workflows.path_layout', 'pyar.conformer.layout_v2'):
        importlib.import_module(name)
    run(sys.executable, '-m', 'pyar.modern_cli', '--help')
    run(sys.executable, '-m', 'pyar.cli', '--help')
    doctor = run(sys.executable, '-m', 'pyar.modern_cli', 'doctor', '--json')
    assert json.loads(doctor.stdout)['ok']
    status = run(sys.executable, '-m', 'pyar.modern_cli', 'install', 'status', '--json')
    inventory = json.loads(status.stdout)
    assert inventory['ok'] and inventory['executable'] == sys.executable
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
        'aggregate': 'aggregate', 'neb': 'scan_path', 'conformer': 'conformer',
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
            result = run(sys.executable, '-m', 'pyar.modern_cli', 'aggregate', str(source),
                         '--size', '2', '--orientations', '2', '--maximum-number-of-seeds', '1',
                         '--selection-feature', 'distance-histogram', '--connectivity-policy', 'off', cwd=directory)
            root = Path(directory) / 'aggregates'
            assert json.loads((root / 'layout.json').read_text())['layout_version'] == 2
            manifest = json.loads((root / 'results/manifest.json').read_text())
            assert manifest['artifacts']
            from pyar.persistence import ArtifactRef
            for ref in manifest['artifacts']:
                ArtifactRef.from_dict({key: ref[key] for key in ('artifact_id', 'role', 'path', 'sha256')}).verify(root)
            growth = Path(directory) / 'growth'
            growth.mkdir()
            run(sys.executable, '-m', 'pyar.modern_cli', 'grow', str(source), str(source),
                '--count', '1', '--orientations', '2', '--maximum-number-of-seeds', '1',
                '--selection-feature', 'distance-histogram', '--connectivity-policy', 'off', cwd=growth)
            assert json.loads((growth / 'grow/layout.json').read_text())['layout_version'] == 2
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
    elif mode == 'conformer':
        from pyar.conformer.layout_v2 import selected_artifacts
        with tempfile.TemporaryDirectory() as directory:
            run(sys.executable, '-m', 'pyar.modern_cli', 'conformer', 'CCCC',
                '--num-conformers', '4', '--top-n', '1', '--num-seeds', '1',
                '--no-torsion-kicks', '--seed', '13', '--num-threads', '1', cwd=directory)
            root = Path(directory) / 'conformers'
            assert selected_artifacts(root)
            assert json.loads((root / 'state.json').read_text())['workflow_state_schema_version'] == 4
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
    elif mode == 'microsolvate':
        with tempfile.TemporaryDirectory() as directory:
            source, solvent = Path(directory) / 'methane.xyz', Path(directory) / 'water.xyz'
            source.write_text('5\nmethane\nC 0 0 0\nH .629 .629 .629\nH -.629 -.629 .629\nH -.629 .629 -.629\nH .629 -.629 -.629\n')
            solvent.write_text('3\nwater\nO 0 0 0\nH .9572 0 0\nH -.239 .927 0\n')
            run(sys.executable, '-m', 'pyar.modern_cli', 'microsolvate', str(source), str(solvent),
                '--count', '1', '--orientations', '2', '--maximum-number-of-seeds', '1',
                '--surface-points-per-atom', '32', '--shell-tolerance', '3.5', cwd=directory)
            root = Path(directory) / 'microsolvation'
            assert json.loads((root / 'layout.json').read_text())['layout_version'] == 2
            assert json.loads((root / 'results/manifest.json').read_text())['artifacts']
    elif mode == 'neb':
        importlib.import_module('pyar.neb')
        run(sys.executable, '-m', 'pyar.neb', '--help')
    print(f'{mode}: installation smoke passed.')


if __name__ == '__main__':
    main()
