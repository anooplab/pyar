"""Behavioral regressions found in the modern CLI review."""
import argparse
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from pyar import modern_cli
from pyar.run_config import (PROFILE_KIND, RUN_CONFIG_KIND, arguments_from_profile,
                             load_toml_document, resolve_run_spec, write_toml_document)
from pyar.scripts import modern_grow
from pyar.workflow_results import GrowResult


@pytest.fixture
def fragments(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for name in ('seed', 'addend'):
        Path(name + '.xyz').write_text('1\nfixture\nHe 0 0 0\n')
    return ['seed.xyz', 'addend.xyz']


def result(request, *, output):
    return GrowResult('grow', 'completed', str(Path(output).resolve()),
                      selected_paths=('selected.xyz',), metadata={'completed_additions': request.count})


def test_saved_check_config_executes_and_cli_overrides_edited_values(fragments, capsys):
    modern_cli.main(['grow', *fragments, '--count', '1', '--check', '--write-config', 'run.toml'])
    config = load_toml_document('run.toml', RUN_CONFIG_KIND)
    assert '--check' not in config['arguments']
    assert config['effective']['multiplicity'] == [1, 1]
    config['effective']['count'] = 2
    write_toml_document('edited.toml', kind=RUN_CONFIG_KIND, workflow='grow',
                        arguments=config['arguments'], effective=config['effective'])
    with patch.object(modern_grow, 'grow', side_effect=result) as engine:
        modern_cli.main(['run', 'edited.toml'])
        assert engine.call_args.args[0].count == 2
        modern_cli.main(['run', 'edited.toml', '--count', '3'])
        assert engine.call_args.args[0].count == 3
    records = list(Path('.pyar/runs').glob('*/pyar-run.toml'))
    assert len(records) == 2
    data = [load_toml_document(path, 'pyar-run-record') for path in records]
    overridden = next(r for r in data if r['effective']['count'] == 3)
    assert overridden['metadata']['value_sources']['count'] == 'cli'
    assert overridden['metadata']['result']['metadata']['completed_additions'] == 3
    assert overridden['metadata']['resolved_request']['electronic_states'][0]['multiplicity'] == 1


def test_dry_run_json_has_actual_resolved_backend_and_no_writes(fragments, capsys):
    before = set(Path('.').iterdir())
    with patch.object(modern_grow, 'preflight', return_value=['xtb']), patch.object(modern_grow, 'grow') as engine:
        modern_cli.main(['grow', *fragments, '--count', '2', '--backend', 'xtb', '--dry-run', '--json'])
    value = json.loads(capsys.readouterr().out)
    backend = value['resolved']['scientific']['backend']
    assert backend['software'] == 'xtb'
    assert backend['xtb_model'] == 'gfn2'
    assert backend['basis'] is None
    assert value['plan']['outputs'] == [str(Path('grow').resolve())]
    assert value['status'] == 'ready'
    engine.assert_not_called()
    assert set(Path('.').iterdir()) == before


def test_xtb_config_replays_without_injecting_unsupported_settings(fragments, capsys):
    with patch.object(modern_grow, 'preflight', return_value=['xtb']):
        modern_cli.main(['grow', *fragments, '--count', '2', '--backend', 'xtb', '--check', '--write-config', 'xtb.toml'])
        modern_cli.main(['run', 'xtb.toml', '--check', '--json'])
    assert json.loads(capsys.readouterr().out.split('Ready to run.\n--check specified; no calculations were performed.\n')[-1])['status'] == 'ready'


def test_preflight_failure_writes_neither_config_nor_run_record(fragments):
    with patch.object(modern_grow, 'preflight', side_effect=ValueError('missing xtb')), pytest.raises(SystemExit):
        modern_cli.main(['grow', *fragments, '--count', '1', '--backend', 'xtb', '--write-config', 'invalid.toml'])
    assert not Path('invalid.toml').exists()
    assert not Path('.pyar').exists()
    assert not Path('grow').exists()


def test_failure_after_start_is_recorded_even_for_parser_exit(fragments):
    with patch.object(modern_grow, 'grow', side_effect=ValueError('backend failed')), pytest.raises(SystemExit):
        modern_cli.main(['grow', *fragments, '--count', '1'])
    path = next(Path('.pyar/runs').glob('*/pyar-run.toml'))
    record = load_toml_document(path, 'pyar-run-record')
    assert record['metadata']['status'] == 'failed'
    assert record['metadata']['resolved_request']['request']['count'] == 1


def test_boolean_profile_required_options_and_cli_precedence():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('--backend', required=True, choices=['xtb', 'orca'])
    parser.add_argument('--torsion-kicks', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('-N', '--orientations', type=int, default=8)
    profile = {'backend': 'xtb', 'torsion_kicks': False, 'orientations': 20}
    _, resolved = resolve_run_spec('conformer', parser, ['CCO', '-N8'], profile_values=profile)
    assert resolved.values['backend'] == 'xtb'
    assert resolved.values['torsion_kicks'] is False
    assert resolved.values['orientations'] == 8
    assert resolved.value_sources['orientations'] == 'cli'
    vector = arguments_from_profile(profile, parser, workflow='conformer', skip_defaults=False)
    assert '--no-torsion-kicks' in vector
    assert 'False' not in vector
    with pytest.raises(SystemExit):
        resolve_run_spec('conformer', parser, ['CCO'], profile_values={'backend': 'unknown'})


def test_option_values_after_delimiter_are_preserved(fragments, capsys):
    Path('--verbose').write_text('1\nenergy=-1\nHe 0 0 0\n')
    modern_cli.main(['energies', '--json', '--', '--verbose'])
    assert json.loads(capsys.readouterr().out)['minimum'] == '--verbose'


def test_profile_list_is_read_only_and_creation_refuses_overwrite(fragments, monkeypatch):
    monkeypatch.setenv('XDG_CONFIG_HOME', str(Path('settings').resolve()))
    modern_cli.main(['profile', 'list'])
    assert not Path('settings').exists()
    record = 'record.toml'
    write_toml_document(record, kind='pyar-run-record', workflow='grow', arguments=[],
                        effective={'count': 2, 'orientations': 8})
    modern_cli.main(['profile', 'create', 'careful', '--from', record])
    target = Path('settings/pyar/profiles/careful.toml')
    before = target.read_bytes()
    with pytest.raises(SystemExit):
        modern_cli.main(['profile', 'create', 'careful', '--from', record])
    assert target.read_bytes() == before


def test_no_abbreviation_can_override_profile_by_accident():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('--orientations', type=int, default=8)
    with pytest.raises(SystemExit):
        resolve_run_spec('grow', parser, ['--orient', '20'], profile_values={'orientations': 10})


def test_json_run_contains_result_without_console_exit_object(fragments, capsys):
    with patch.object(modern_grow, 'grow', side_effect=result):
        assert modern_cli.main(['grow', *fragments, '--count', '2', '--json']) is None
    streams = capsys.readouterr()
    value = json.loads(streams.out)
    assert value['result']['status'] == 'completed'
    assert 'Growth completed' in streams.err


def test_verbosity_does_not_leak_between_calls(fragments):
    import logging
    logger = logging.getLogger('pyar')
    old_level, old_handlers = logger.level, list(logger.handlers)
    modern_cli.main(['-vv', 'grow', *fragments, '--count', '1', '--check'])
    assert logger.level == old_level and logger.handlers == old_handlers


def test_editable_toml_and_validation_reject_scientific_errors(fragments, capsys):
    modern_cli.main(['grow', *fragments, '--count', '1', '--check', '--write-config', 'run.toml'])
    text = Path('run.toml').read_text()
    assert '[settings]' in text
    Path('run.toml').write_text(text.replace('"count" = 1', '"count" = 0'))
    with patch.object(modern_grow, 'grow') as engine, pytest.raises(SystemExit):
        modern_cli.main(['config', 'validate', 'run.toml'])
    engine.assert_not_called()
    assert not Path('grow').exists()


def test_latest_run_record_and_backend_detail(tmp_path, capsys):
    from pyar.scripts.cli_meta import _load_run_record
    from pyar.run_config import RUN_RECORD_KIND
    import os
    for name, status in [('pyar-run.toml', 'failed'), ('pyar-run-new.toml', 'completed')]:
        write_toml_document(tmp_path / name, kind=RUN_RECORD_KIND, workflow='grow', arguments=[], metadata={'status': status})
    os.utime(tmp_path / 'pyar-run.toml', ns=(1, 1))
    assert _load_run_record(tmp_path)['metadata']['status'] == 'completed'
    modern_cli.main(['backends', 'xtb', '--json'])
    value = json.loads(capsys.readouterr().out)
    assert list(value['backends']) == ['xtb']
    assert 'availability' in value['backends']['xtb']


def test_quiet_run_config_does_not_restore_progress(fragments, capsys):
    modern_cli.main(['grow', *fragments, '--count', '1', '--check', '--write-config', 'run.toml'])
    capsys.readouterr()
    with patch.object(modern_grow, 'grow', side_effect=result):
        modern_cli.main(['run', 'run.toml', '--quiet'])
    assert capsys.readouterr().out == ''


def test_formula_profile_can_be_reused_with_new_building_blocks(tmp_path, monkeypatch, capsys):
    from pyar.scripts import modern_aggregate
    from pyar.run_config import RUN_RECORD_KIND, profile_path
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('XDG_CONFIG_HOME', str(tmp_path / 'config'))
    original = vars(modern_aggregate.build_parser().parse_args(['--formula', 'He2']))
    write_toml_document('record.toml', kind=RUN_RECORD_KIND, workflow='aggregate',
                        arguments=['--formula', 'He2'], effective=original)
    modern_cli.main(['profile', 'create', 'geometry', '--from', 'record.toml'])
    profile = load_toml_document(profile_path('geometry'), PROFILE_KIND)
    assert 'formula' not in profile['effective']
    with patch.object(modern_aggregate, 'aggregate') as engine:
        modern_cli.main(['aggregate', 'Ne', '--size', '2', '--profile', 'geometry', '--check'])
    engine.assert_not_called()
    assert 'Ready to run.' in capsys.readouterr().out
    assert not Path('aggregates').exists()


def test_profiles_exclude_analysis_artifact_paths():
    from pyar.run_config import profile_overrides_from_record
    from pyar.scripts import clustering
    values = {'first': 'a.xyz', 'second': 'b.xyz', 'path': 'old-run',
              'labels_output': 'old-labels.csv', 'report_output': 'old-report.json',
              'structure_report': 'old-structure.json', 'plot_directory': 'old-plots',
              'algorithm': 'auto'}
    assert profile_overrides_from_record({'workflow': 'clustering', 'effective': values}) == {'algorithm': 'auto'}
    assert arguments_from_profile(values, clustering.build_parser(), workflow='clustering', skip_defaults=False) == ['--algorithm=auto']


def test_saved_reaction_config_freezes_resolved_protocol(fragments, capsys):
    from pyar.scripts import modern_react
    with patch.object(modern_react, 'preflight_reaction', return_value=[]):
        modern_cli.main(['react', *fragments, '--backend', 'xtb', '--bias-max', '100',
                         '--check', '--write-config', 'reaction.toml'])
        saved = load_toml_document('reaction.toml', RUN_CONFIG_KIND)
        effective = saved['effective']
        assert effective['geometry_optimizer'] == 'geometric'
        assert effective['bias_controller'] == 'adaptive'
        assert effective['bias_potential'] == 'afir'
        assert effective['release_retry_limit'] is not None
        # Explicit frozen defaults must be accepted by the canonical resolver.
        modern_cli.main(['run', 'reaction.toml', '--check'])
    assert not Path('reaction').exists()


def test_saved_growth_config_freezes_xtb_model_and_optimizer(fragments):
    with patch.object(modern_grow, 'preflight', return_value=[]):
        modern_cli.main(['grow', *fragments, '--count', '1', '--backend', 'xtb',
                         '--check', '--write-config', 'grow.toml'])
        effective = load_toml_document('grow.toml', RUN_CONFIG_KIND)['effective']
        assert effective['geometry_optimizer'] == 'native'
        assert effective['xtb_model'] == 'gfn2'
        assert effective['scf_threshold'] is None
        modern_cli.main(['run', 'grow.toml', '--check'])
