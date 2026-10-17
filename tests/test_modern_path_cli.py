"""Modern path tasks must stop at their named stage and validate before writes."""
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from pyar import modern_cli, neb
from pyar.scripts import modern_path_stage as cli
from pyar.scripts import scan_bond


@pytest.fixture
def endpoints(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for name in ('a', 'b', 'guess'):
        Path(name + '.xyz').write_text('2\nH2\nH 0 0 0\nH 0 0 0.74\n')
    return ['neb', 'a.xyz', 'b.xyz', '--ts-guess', 'guess.xyz', '--backend', 'xtb']


def test_neb_runs_relax_and_neb_only(endpoints):
    relaxed = {key: True for key in neb._STAGE_GATES['relax']}
    with patch.object(cli, 'preflight', return_value=['xtb', 'geometric']), patch.object(cli, 'run_neb', side_effect=[relaxed, {'converged': True, 'interior_maximum': True}]) as engine:
        modern_cli.main(endpoints)
    assert [call.kwargs['stage'] for call in engine.call_args_list] == ['relax', 'neb']
    assert all(call.kwargs['reuse'] is True for call in engine.call_args_list)
    assert engine.call_args.kwargs['method'] is None
    assert engine.call_args.kwargs['xtb_model'] == 'gfn2'


def test_cycle_defaults_are_raised_across_path_entry_points():
    neb_args = cli.build_parser('neb').parse_args([
        'a.xyz', 'b.xyz', '--ts-guess', 'guess.xyz', '--backend', 'xtb',
    ])
    ts_args = cli.build_parser('ts').parse_args(['guess.xyz', '--backend', 'xtb'])
    irc_args = cli.build_parser('irc').parse_args(['run', '--backend', 'xtb'])
    scan_args = scan_bond.build_parser(modern=True).parse_args([
        'a.xyz', 'b.xyz', '--atoms', '0', '0', '--backend', 'xtb',
    ])
    assert neb_args.max_cycles == 1000
    assert ts_args.ts_max_cycles == 1000
    assert irc_args.irc_max_cycles == 10000
    assert scan_args.opt_cycles == 1000
    assert scan_args.max_cycles == 1000
    assert scan_args.product_relaxation_max_steps == 1000
    assert scan_args.ts_max_cycles == 1000
    assert scan_args.irc_max_cycles == 10000
    assert scan_args.endpoint_max_cycles == 1000


def test_neb_stops_on_failed_endpoint_gate(endpoints):
    with patch.object(cli, 'preflight', return_value=[]), patch.object(cli, 'run_neb', return_value={}) as engine, pytest.raises(SystemExit) as exc:
        modern_cli.main(endpoints)
    assert exc.value.code == 1
    assert engine.call_count == 1
    assert engine.call_args.kwargs['stage'] == 'relax'


def test_neb_path_seed_is_an_alternative_initializer_for_check(endpoints, capsys):
    seed = Path('seed.xyz')
    seed.write_text(''.join(
        f'2\nframe={index}\nH 0 0 0\nH 0 0 {distance}\n'
        for index, distance in enumerate((0.74, 0.75, 0.76))
    ))
    arguments = ['neb', 'a.xyz', 'b.xyz', '--path-seed', str(seed), '--backend', 'xtb', '--check']
    with patch.object(cli, 'preflight', return_value=['xtb', 'geometric']), \
            patch.object(cli, 'run_neb') as engine:
        modern_cli.main(arguments)
    assert 'Ready to run' in capsys.readouterr().out
    engine.assert_not_called()


def test_neb_initializer_options_are_validated_before_execution(endpoints):
    with patch.object(cli, 'preflight', return_value=[]), patch.object(cli, 'run_neb') as engine:
        with pytest.raises(SystemExit):
            modern_cli.main(['neb', 'a.xyz', 'b.xyz', '--backend', 'xtb', '--check'])
        with pytest.raises(SystemExit):
            modern_cli.main(endpoints + ['--path-seed', 'guess.xyz', '--check'])
        engine.assert_not_called()


@pytest.mark.parametrize('options', [['--images', '4'], ['--max-cycles', '0'], ['--spring', 'nan'], ['--multiplicity', '2']])
def test_check_validates_without_mutation(endpoints, options):
    before = set(Path('.').iterdir())
    with patch.object(cli, 'preflight', return_value=[]), patch.object(cli, 'run_neb') as engine, pytest.raises(SystemExit) as exc:
        modern_cli.main(endpoints + options + ['--check'])
    assert exc.value.code == 2
    engine.assert_not_called()
    assert set(Path('.').iterdir()) == before


def test_valid_check_and_output_alias(endpoints, capsys):
    with patch.object(cli, 'preflight', return_value=['xtb']), patch.object(cli, 'run_neb') as engine:
        modern_cli.main(endpoints + ['--check', '--json'])
        assert json.loads(capsys.readouterr().out)['status'] == 'ready'
        Path('ts_guess.xyz').write_text(Path('guess.xyz').read_text())
        arguments = [*endpoints]
        arguments[arguments.index('guess.xyz')] = 'ts_guess.xyz'
        with pytest.raises(SystemExit):
            modern_cli.main(arguments + ['--output', '.', '--check'])
        engine.assert_not_called()
    assert not Path('neb_run').exists()


def test_irc_check_verifies_artifact_hashes(endpoints):
    root = Path('prior'); root.mkdir()
    geometry = root / 'frequency_geometry.xyz'
    geometry.write_text(Path('guess.xyz').read_text())
    from types import SimpleNamespace
    calculator = SimpleNamespace(qc_params=dict(software='xtb', method=None, basis=None,
                                               charge=0, multiplicity=1, nprocs=1,
                                               gamma=0.0, xtb_model='gfn2'))
    neb._save_stage(root, 'frequency', {'first_order_saddle_confirmed': True}, calculator,
                    inputs=[], artifacts=[geometry])
    before = {p: p.read_bytes() for p in root.iterdir()}
    with patch.object(cli, 'preflight', return_value=[]), patch.object(cli, 'run_neb') as engine:
        modern_cli.main(['irc', str(root), '--backend', 'xtb', '--check'])
        assert before == {p: p.read_bytes() for p in root.iterdir()}
        geometry.write_text(geometry.read_text().replace('0.74', '0.75'))
        with pytest.raises(SystemExit):
            modern_cli.main(['irc', str(root), '--backend', 'xtb', '--check'])
        engine.assert_not_called()
