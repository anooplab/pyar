"""Modern path tasks must stop at their named stage and validate before writes."""
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from pyar import modern_cli, neb
from pyar.scripts import modern_path_stage as cli


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
    assert engine.call_args.kwargs['method'] is None
    assert engine.call_args.kwargs['xtb_model'] == 'gfn2'


def test_neb_stops_on_failed_endpoint_gate(endpoints):
    with patch.object(cli, 'preflight', return_value=[]), patch.object(cli, 'run_neb', return_value={}) as engine, pytest.raises(SystemExit) as exc:
        modern_cli.main(endpoints)
    assert exc.value.code == 1
    assert engine.call_count == 1
    assert engine.call_args.kwargs['stage'] == 'relax'


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
