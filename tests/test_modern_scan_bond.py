"""Modern scan request validation, stage-aware preflight and legacy isolation."""

from pathlib import Path
from unittest import mock

import pytest

from pyar import modern_cli
from pyar.scripts import scan_bond


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('a.xyz').write_text('2\nA\nH 0 0 0\nH 0 0 0.74\n')
    Path('b.xyz').write_text('2\nB\nH 0 0 0\nH 0 0 0.85\n')
    return ['scan-bond', 'a.xyz', 'b.xyz', '--atoms', '0', '1']


@pytest.fixture
def dependencies():
    with mock.patch('pyar.cli.shutil.which', return_value='/bin/backend'):
        yield


def test_help(capsys):
    for argv in [['--help'], ['scan-bond', '--help']]:
        with pytest.raises(SystemExit) as exc:
            modern_cli.main(argv)
        assert exc.value.code == 0
        text = capsys.readouterr().out
        assert 'scan-bond' in text
        if len(argv) == 1:
            assert 'optimize' in text and 'clustering' in text
        else:
            assert '0-based' in text and 'local' in text
            assert '--backend' in text and '--software' not in text
            for option in ['--through', '--orientations', '--check', '--gxtb-wrapper',
                           '--interpolation', '--sella-internal-coordinates', '--irc-max-cycles']:
                assert option in text


@pytest.mark.parametrize('override, count', [([], 8), (['--orientations', '20'], 20), (['-N', '20'], 20)])
def test_minimal_and_orientations(inputs, dependencies, override, count):
    with mock.patch.object(scan_bond, 'run_scan_bond') as run:
        modern_cli.main([*inputs, '--backend', 'xtb', *override])
    args = run.call_args.args
    assert args[:4] == ('a.xyz', 'b.xyz', [0, 1], count)
    assert args[4]['software'] == 'xtb' and args[4]['xtb_model'] == 'gfn2'
    assert args[4]['method'] is None and args[4]['basis'] is None
    assert args[4]['charge_a'] == args[4]['charge_b'] == 0
    assert args[4]['multiplicity_a'] == args[4]['multiplicity_b'] == 1
    assert args[9] == 'scan'
    assert args[6:9] == (None, None, None)


@pytest.mark.parametrize('atoms', [['-1', '0'], ['2', '0'], ['0', '-1'], ['0', '2']])
def test_invalid_atoms(inputs, dependencies, atoms, capsys):
    with mock.patch.object(scan_bond, 'run_scan_bond') as run:
        with pytest.raises(SystemExit):
            modern_cli.main([*inputs[:3], '--atoms', *atoms, '--backend', 'xtb'])
    run.assert_not_called()
    assert 'out of range for fragment' in capsys.readouterr().err


@pytest.mark.parametrize('options', [
    [], ['--backend', 'invalid'],
    ['--backend', 'orca', '--method', 'BP86'],
    ['--backend', 'orca', '--method', 'XTB2', '--basis', 'def2-SVP'],
    ['--backend', 'orca', '--method', 'g-xTB'],
    ['--backend', 'gaussian'],
    ['--backend', 'xtb', '--method', 'BP86'],
    ['--backend', 'xtb', '--basis', 'def2-SVP'],
    ['--backend', 'xtb', '--multiplicity', '2'],
    ['--backend', 'xtb', '--multiplicity', '0'],
    ['--backend', 'xtb', '--charge', '0', '0', '0'],
    ['--backend', 'xtb', '--scan-step', '0.1', '--scan-points', '10'],
    ['--backend', 'xtb', '--scan-step', 'nan'],
    ['--backend', 'xtb', '--scan-points', '1'],
    ['--backend', 'xtb', '--through', 'all', '--images', '4'],
    ['--backend', 'orca', '--method', 'XTB2', '--xtb-model', 'gxtb'],
])
def test_invalid_request(inputs, dependencies, options):
    with mock.patch.object(scan_bond, 'run_scan_bond') as run:
        with pytest.raises(SystemExit) as exc:
            modern_cli.main([*inputs, *options])
    assert exc.value.code != 0
    run.assert_not_called()


@pytest.mark.parametrize('contents', [None, 'bad XYZ', '1\nB\nH nan 0 0\n'])
def test_invalid_input(inputs, dependencies, contents):
    Path('b.xyz').unlink()
    if contents is not None:
        Path('b.xyz').write_text(contents)
    with mock.patch.object(scan_bond, 'run_scan_bond') as run:
        with pytest.raises(SystemExit):
            modern_cli.main([*inputs, '--backend', 'xtb'])
    run.assert_not_called()


def test_electronic_inference_and_pair_overrides(inputs, dependencies):
    Path('b.xyz').write_text('1\nB\nH 0 0 0\n')
    base = [*inputs[:3], '--atoms', '1', '0', '--backend', 'xtb']
    with mock.patch.object(scan_bond, 'run_scan_bond') as run:
        modern_cli.main(base)
        params = run.call_args.args[4]
        assert params['multiplicity_a'] == 1 and params['multiplicity_b'] == 2
        modern_cli.main([*base, '--charge', '0', '-1'])
        params = run.call_args.args[4]
        assert params['charge_b'] == -1 and params['multiplicity_b'] == 1
        modern_cli.main([*base, '--multiplicity', '1', '2'])


def test_two_radical_fragments_keep_scan_scf_default(inputs, dependencies):
    for path in ['a.xyz', 'b.xyz']:
        Path(path).write_text('1\nH radical\nH 0 0 0\n')
    with mock.patch.object(scan_bond, 'run_scan_bond') as run:
        modern_cli.main([*inputs[:3], '--atoms', '0', '0', '--backend', 'xtb'])
    params = run.call_args.args[4]
    assert params['multiplicity_a'] == params['multiplicity_b'] == 2
    assert params['scftype_a'] == params['scftype_b'] == 'rhf'


@pytest.mark.parametrize('through', ['scan', 'neb', 'ts', 'frequency', 'irc', 'endpoints', 'all'])
def test_check_no_calculation_or_mutations(inputs, dependencies, through, capsys):
    before = set(Path('.').iterdir())
    with mock.patch.object(scan_bond, 'run_scan_bond') as run, mock.patch(
        'subprocess.Popen', side_effect=AssertionError('check must not launch programs')
    ):
        modern_cli.main([*inputs, '--backend', 'xtb', '--through', through, '--check'])
    run.assert_not_called()
    assert set(Path('.').iterdir()) == before
    assert 'Ready to run' in capsys.readouterr().out


def test_stage_optional_dependencies(inputs, dependencies):
    with mock.patch('pyar.neb._sella_api', side_effect=RuntimeError('install pyar-chem[sella]')), mock.patch(
        'pyar.neb._geodesic_api', side_effect=RuntimeError('install pyar-chem[geodesic]')
    ), mock.patch.object(scan_bond, 'run_scan_bond') as run:
        modern_cli.main([*inputs, '--backend', 'xtb', '--ts-optimizer', 'sella',
                         '--interpolation', 'geodesic', '--check'])
        for options in [['--through', 'ts', '--ts-optimizer', 'sella'],
                        ['--through', 'neb', '--interpolation', 'geodesic']]:
            with pytest.raises(SystemExit):
                modern_cli.main([*inputs, '--backend', 'xtb', *options, '--check'])
    run.assert_not_called()


def test_missing_backend(inputs, capsys):
    with mock.patch('pyar.cli.shutil.which', return_value=None), mock.patch.object(scan_bond, 'run_scan_bond') as run:
        with pytest.raises(SystemExit):
            modern_cli.main([*inputs, '--backend', 'xtb', '--check'])
    run.assert_not_called()
    assert 'xtb-docs.readthedocs.io' in capsys.readouterr().err


def test_orca_scan_only_does_not_need_geometric(inputs, dependencies):
    with mock.patch('pyar.cli.importlib.util.find_spec', return_value=None):
        modern_cli.main([*inputs, '--backend', 'orca', '--method', 'XTB2', '--check'])
        with pytest.raises(SystemExit):
            modern_cli.main([*inputs, '--backend', 'orca', '--method', 'XTB2', '--through', 'neb', '--check'])


def test_wrapper_validation(inputs, dependencies, tmp_path):
    base = [*inputs, '--backend', 'orca', '--method', 'g-xTB', '--gxtb-wrapper']
    wrapper = tmp_path / 'wrapper'
    for path in [tmp_path / 'missing', wrapper]:
        if path == wrapper:
            wrapper.write_text('#!/bin/sh\nexit 0\n')
        with pytest.raises(SystemExit):
            modern_cli.main([*base, str(path), '--check'])
    wrapper.chmod(0o755)
    modern_cli.main([*base, str(wrapper), '--check'])
    with pytest.raises(SystemExit):
        modern_cli.main([*base, str(wrapper), '--through', 'all', '--check'])
    with pytest.raises(SystemExit):
        modern_cli.main([*inputs, '--backend', 'orca', '--method', 'XTB2', '--gxtb-wrapper', str(wrapper), '--check'])


def test_legacy_parser_contract(inputs):
    parser = scan_bond.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(inputs[1:])
    args = parser.parse_args([*inputs[1:], '-N', '1'])
    assert args.software == 'orca' and args.multiplicity == [1]
    assert not hasattr(args, 'check')


def test_legacy_and_modern_share_scientific_arguments(inputs, dependencies):
    with mock.patch.object(scan_bond, 'run_scan_bond') as run:
        scan_bond.main([*inputs[1:], '--software', 'xtb', '-N', '8', '--through', 'neb',
                        '--images', '13', '--scan-points', '5', '--output', 'my_scan'])
        legacy_args = run.call_args.args
        modern_cli.main([*inputs, '--backend', 'xtb', '--through', 'neb',
                         '--images', '13', '--scan-points', '5', '--output', 'my_scan'])
        assert legacy_args == run.call_args.args


def test_structured_failure_stays_nonzero(inputs, dependencies):
    with mock.patch.object(scan_bond, 'run_scan_bond', return_value={'status': 'failed', 'output_dir': 'scan_bond'}):
        with pytest.raises(SystemExit) as exc:
            modern_cli.main([*inputs, '--backend', 'xtb'])
    assert exc.value.code == 1


def test_successful_workflow_is_not_an_exit_value(inputs, dependencies):
    with mock.patch.object(scan_bond, 'run_scan_bond', return_value={'status': 'complete', 'output_dir': 'scan_bond'}):
        assert modern_cli.main([*inputs, '--backend', 'xtb']) is None
