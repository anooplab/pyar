"""Fail-fast modern bulk optimization and legacy isolation."""

from pathlib import Path
from unittest import mock

import pytest

from pyar import modern_cli, optimiser
from pyar.scripts import optimize
from pyar import optimization_request as request


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for name, atom in [('a', 'He'), ('b', 'He'), ('radical', 'H')]:
        Path(name + '.xyz').write_text(f'1\n{name}\n{atom} 0 0 0\n')
    return ['a.xyz', 'b.xyz']


@pytest.fixture
def dependencies():
    with mock.patch('pyar.cli.shutil.which', return_value='/bin/backend'):
        yield


def test_help(capsys):
    for argv in [['--help'], ['optimize', '--help']]:
        with pytest.raises(SystemExit) as exc:
            modern_cli.main(argv)
        assert exc.value.code == 0
        text = capsys.readouterr().out
        assert 'optimize' in text
        if argv == ['--help']:
            assert 'clustering' in text
        else:
            assert '--backend' in text and '--check' in text
            assert '--software' not in text


def test_minimal_and_inferred_spin(inputs, dependencies):
    with mock.patch.object(optimiser, 'bulk_optimize', side_effect=lambda mols, qc: mols) as bulk:
        modern_cli.main(['optimize', *inputs, 'radical.xyz', '--backend', 'xtb'])
    molecules, settings = bulk.call_args.args
    assert [mol.name for mol in molecules] == ['a', 'b', 'radical']
    assert [mol.charge for mol in molecules] == [0, 0, 0]
    assert [mol.multiplicity for mol in molecules] == [1, 1, 2]
    assert settings['software'] == 'xtb'
    assert settings['xtb_model'] == 'gfn2'
    assert settings['method'] is None and settings['basis'] is None
    assert settings['nprocs'] == 1


def test_overrides_and_engine_electronic_state(inputs, dependencies):
    calls = []
    def run(mol, qc):
        calls.append((mol, qc))
        return True
    with mock.patch.object(optimiser, 'optimise', side_effect=run):
        modern_cli.main(['optimize', *inputs, '--backend', 'xtb', '--charge', '-1',
                         '--multiplicity', '2', '--opt-threshold', 'tight', '--nprocs', '4'])
    assert len(calls) == 2
    assert calls[0][1] == calls[1][1]
    for mol, settings in calls:
        assert mol.charge == settings['charge'] == -1
        assert mol.multiplicity == settings['multiplicity'] == 2
        assert settings['xtb_unpaired_electrons'] == 1
        assert settings['opt_threshold'] == 'tight' and settings['nprocs'] == 4


def test_check_no_mutations(inputs, dependencies, capsys):
    before = set(Path('.').iterdir())
    with mock.patch.object(optimiser, 'bulk_optimize') as bulk, mock.patch.object(
        optimize, 'preflight', wraps=request.preflight
    ) as preflight:
        modern_cli.main(['optimize', *inputs, '--backend', 'xtb', '--check'])
    bulk.assert_not_called()
    preflight.assert_called_once()
    assert set(Path('.').iterdir()) == before
    assert 'Ready to run' in capsys.readouterr().out


@pytest.mark.parametrize('options', [
    [], ['--backend', 'unknown'], ['--backend', 'ani'],
    ['--backend', 'xtb', '--multiplicity', '2'],
    ['--backend', 'xtb', '--multiplicity', '0'],
    ['--backend', 'xtb', '--basis', 'def2-SVP'],
    ['--backend', 'xtb', '--nprocs', '0'],
    ['--backend', 'xtb', '--opt-target', 'ts'],
    ['--backend', 'psi4', '--geometry-optimizer', 'geometric'],
    ['--backend', 'orca', '--method', 'r2scan-3c', '--geometry-optimizer', 'geometric'],
    ['--backend', 'xtb', '--custom-keywords', 'Opt'],
    ['--backend', 'gaussian', '--geometry-optimizer', 'native'],
    ['--backend', 'orca', '--method', ''],
])
def test_invalid_settings_no_execution(inputs, dependencies, options):
    with mock.patch.object(optimiser, 'bulk_optimize') as bulk:
        with pytest.raises(SystemExit) as exc:
            modern_cli.main(['optimize', *inputs, *options])
    assert exc.value.code != 0
    bulk.assert_not_called()


@pytest.mark.parametrize('contents', [None, 'bad XYZ', '1\ninvalid\nHe nan 0 0\n'])
def test_all_inputs_validated_before_execution(inputs, dependencies, contents):
    if contents is not None:
        Path('bad.xyz').write_text(contents)
    with mock.patch.object(optimiser, 'bulk_optimize') as bulk:
        with pytest.raises(SystemExit):
            modern_cli.main(['optimize', *inputs, 'bad.xyz', '--backend', 'xtb'])
    bulk.assert_not_called()


@pytest.mark.parametrize('check', [[], ['--check']])
def test_missing_executable(inputs, check, capsys):
    with mock.patch('pyar.cli.shutil.which', return_value=None), mock.patch.object(
        optimiser, 'bulk_optimize'
    ) as bulk:
        with pytest.raises(SystemExit):
            modern_cli.main(['optimize', *inputs, '--backend', 'xtb', *check])
    bulk.assert_not_called()
    text = capsys.readouterr().err
    assert 'xtb-docs.readthedocs.io' in text and 'No calculations' in text


def test_import_broken_dependency(inputs, dependencies):
    with mock.patch.object(request, 'import_module', side_effect=ImportError('broken torch')):
        with mock.patch('pyar.backends.aimnet2_assets.validate_aimnet2_runtime_assets'):
            with mock.patch.object(optimiser, 'bulk_optimize') as bulk:
                with pytest.raises(SystemExit):
                    modern_cli.main(['optimize', *inputs, '--backend', 'aimnet_2', '--check'])
    bulk.assert_not_called()


def test_orca_composite_and_overrides(inputs, dependencies):
    from pyar.backends.orca_methods import orca_method_keywords
    with mock.patch.object(optimiser, 'bulk_optimize', side_effect=lambda mols, qc: mols) as bulk:
        modern_cli.main(['optimize', *inputs, '--backend', 'orca', '--method', 'r2scan-3c',
                         '--nprocs', '8', '--opt-cycles', '200'])
    settings = bulk.call_args.args[1]
    assert settings['basis'] is None
    assert orca_method_keywords(settings, 'Opt')[0] == '! r2scan-3c Opt'
    assert settings['opt_cycles'] == 200 and settings['nprocs'] == 8


def test_legacy_bulk_parameters_unchanged():
    settings = {'software': 'xtb'}
    mol = mock.Mock()
    with mock.patch.object(optimiser, 'optimise', return_value=True) as single:
        assert optimiser.bulk_optimize([mol], settings) == [mol]
    assert single.call_args.args[1] is settings


def test_xtb_spin_mapping_is_opt_in():
    from pyar.backends.xtb_utils import build_xtb_command
    settings = {'multiplicity': 2}
    assert build_xtb_command('xtb', 'a.xyz', settings)[-2:] == ['-uhf', '2']
    settings['xtb_unpaired_electrons'] = 1
    assert build_xtb_command('xtb', 'a.xyz', settings)[-2:] == ['-uhf', '1']


def test_missing_python_requirement(inputs, dependencies, capsys):
    with mock.patch('pyar.cli.importlib.util.find_spec', return_value=None), mock.patch.object(
        optimiser, 'bulk_optimize'
    ) as bulk:
        with pytest.raises(SystemExit):
            modern_cli.main(['optimize', *inputs, '--backend', 'xtb',
                             '--geometry-optimizer', 'geometric', '--check'])
    bulk.assert_not_called()
    assert 'geometric' in capsys.readouterr().err


def test_per_molecule_inference_reaches_engine(inputs, dependencies):
    with mock.patch.object(optimiser, 'optimise', return_value=True) as single:
        modern_cli.main(['optimize', 'a.xyz', 'radical.xyz', '--backend', 'xtb'])
    assert [call.args[1]['multiplicity'] for call in single.call_args_list] == [1, 2]


def test_backend_alias(inputs, dependencies):
    with mock.patch.object(optimiser, 'bulk_optimize', side_effect=lambda mols, qc: mols) as bulk:
        modern_cli.main(['optimize', *inputs, '--backend', 'g16'])
    assert bulk.call_args.args[1]['software'] == 'gaussian'
    assert bulk.call_args.args[1]['method'] == 'BP86'
    assert bulk.call_args.args[1]['geometry_optimizer'] == 'geometric'


def test_duplicate_names_fail(inputs, dependencies, tmp_path):
    (tmp_path / 'other').mkdir()
    (tmp_path / 'other' / 'a.xyz').write_text(Path('a.xyz').read_text())
    with mock.patch.object(optimiser, 'bulk_optimize') as bulk:
        with pytest.raises(SystemExit):
            modern_cli.main(['optimize', 'a.xyz', 'other/a.xyz', '--backend', 'xtb'])
    bulk.assert_not_called()


def test_legacy_optimizer_contract(inputs, capsys):
    # Import only inside the temporary work directory: legacy logging stays legacy.
    from pyar.scripts import optimiser as legacy
    try:
        with mock.patch('sys.argv', ['pyar-optimiser', *inputs, '--software', 'xtb']):
            with pytest.raises(SystemExit) as exc:
                legacy.main()
        assert exc.value.code == 2
        with mock.patch('sys.argv', ['pyar-optimiser', *inputs, '--software', 'xtb',
                                    '--charge', '0', '--multiplicity', '1']):
            with mock.patch.object(optimiser, 'bulk_optimize') as bulk:
                legacy.main()
        molecules, settings = bulk.call_args.args
        assert len(molecules) == 2
        assert settings['software'] == 'xtb'
        assert settings['method'] == 'BP86' and settings['basis'] == 'def2-SVP'
        assert '_electronic_state_from_molecule' not in settings
    finally:
        legacy.logger.removeHandler(legacy.handler)
        legacy.handler.close()
