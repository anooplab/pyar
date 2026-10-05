import json
from pathlib import Path
from unittest.mock import patch

import pytest

from pyar.modern_cli import main
from pyar.reaction_request import resolve_reaction_request, validate_restart, preflight_reaction
from pyar.state.reaction import ReactionRunState, ReactionStateError
from pyar.workflow_results import ReactionResult


def xyz(tmp_path, name, symbol='H'):
    path = tmp_path/name
    path.write_text(f'1\nreactant\n{symbol} 0 0 0\n')
    return str(path)


def inputs(tmp_path):
    return [xyz(tmp_path, 'a.xyz'), xyz(tmp_path, 'b.xyz')]


def resolve(files, **kwargs):
    return resolve_reaction_request(files, dict(backend='xtb', bias_max=100, **kwargs))


@pytest.fixture
def workflow():
    with patch('pyar.scripts.modern_react.preflight_reaction', return_value=['xtb', 'geometric', 'obabel']), patch(
            'pyar.workflows.reaction.react', return_value=ReactionResult(
                'reaction', 'completed_no_products', 'reaction', 'reaction/state.json')) as run:
        yield run


def test_minimal_defaults_and_native_state(tmp_path, monkeypatch, workflow, capsys):
    files = inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    main(['react', *files, '--backend', 'xtb', '--bias-max', '100'])
    args = workflow.call_args.args
    settings = args[5]
    assert args[2:5] == (None, 100, 8)
    assert args[6:] == (None, 2.3)
    assert settings['bias_controller'] == 'adaptive' and settings['bias_potential'] == 'afir'
    assert settings['geometry_optimizer'] == 'geometric' and settings['software'] == 'xtb'
    assert settings['method'] is None and settings['basis'] is None
    assert settings['xtb_model'] == 'gfn2'
    assert settings['charge'] == 0 and settings['multiplicity'] == 1  # Two doublets combine to a singlet.
    assert settings['trace_enabled']
    assert all(m.multiplicity == 2 for m in args[:2])
    assert 'Products found: 0' in capsys.readouterr().out


@pytest.mark.parametrize('count', [0, 1, 3])
def test_exactly_two(tmp_path, count, workflow):
    files = inputs(tmp_path)
    with pytest.raises(SystemExit) as exc:
        main(['react', *(files+[files[0]])[:count], '--backend', 'xtb', '--bias-max', '100'])
    assert exc.value.code == 2 and not workflow.called


def test_missing_ceiling_and_backend(tmp_path, workflow, capsys):
    files = inputs(tmp_path)
    for options in [['--backend', 'xtb'], ['--bias-max', '100']]:
        with pytest.raises(SystemExit):
            main(['react', *files, *options])
        assert not workflow.called
    captured = capsys.readouterr().err
    assert '--bias-max' in captured and 'backend' in captured


@pytest.mark.parametrize('options', [dict(bias_controller='fixed'), dict(bias_controller='scheduled'),
                                    dict(bias_alpha_margin=0), dict(softmin_beta=-1),
                                    dict(site=[-1, 0]), dict(site=[0, 1]), dict(bias_min=float('nan')),
                                    dict(geometry_optimizer='native'), dict(multiplicity=[1]),
                                    dict(release_margin_factor=1), dict(release_distance_fraction=2),
                                    dict(method='BP86')])
def test_invalid_request_no_work(tmp_path, options):
    files = inputs(tmp_path)
    with pytest.raises(ValueError):
        resolve(files, **options)
    assert not (tmp_path/'reaction').exists()


def test_controllers_softmin_site_and_spin(tmp_path):
    files = inputs(tmp_path)
    default = resolve(files)
    explicit = resolve(files, bias_controller='adaptive')
    assert default.restart_request == explicit.restart_request
    fixed = resolve(files, bias_controller='fixed', bias_min=10)
    assert fixed.qc_params['bias_controller'] == 'fixed' and fixed.bias_min == 10
    assert fixed.restart_request['gamma_schedule'][0] == 10
    assert fixed.restart_request['gamma_schedule'][-1] == 100
    scheduled = resolve(files, bias_controller='scheduled', bias_min=10, bias_scheduled_alpha=.002)
    assert scheduled.qc_params['bias_scheduled_alpha'] == .002
    tuned = resolve(files, bias_potential='softmin', softmin_beta=2., site=[0, 0], charge=[-1, 0],
                    multiplicity=[1, 2])
    assert tuned.site == [0, 1]
    assert tuned.qc_params['softmin_beta'] == 2 and tuned.qc_params['bias_potential'] == 'softmin'
    assert tuned.qc_params['charge'] == -1 and tuned.qc_params['multiplicity'] == 2
    assert tuned.qc_params['xtb_unpaired_electrons'] == 1
    neutral_even = resolve([xyz(tmp_path, 'he.xyz', 'He'), files[0]])
    assert [m.multiplicity for m in neutral_even.reactants] == [1, 2]


def test_one_charge_applies_to_both(tmp_path):
    resolved = resolve(inputs(tmp_path), charge=[-1])
    assert [m.charge for m in resolved.reactants] == [-1, -1]
    assert [m.multiplicity for m in resolved.reactants] == [1, 1]
    assert resolved.qc_params['charge'] == -2


def test_capabilities_and_qc_settings(tmp_path):
    files = inputs(tmp_path)
    with pytest.raises(ValueError, match='Cartesian'):
        resolve_reaction_request(files, dict(backend='psi4', bias_max=100))
    with pytest.raises(ValueError, match='--method'):
        resolve_reaction_request(files, dict(backend='orca', bias_max=100))
    with pytest.raises(ValueError, match='basis'):
        resolve_reaction_request(files, dict(backend='orca', bias_max=100, method='BP86'))
    orca = resolve_reaction_request(files, dict(backend='orca16', bias_max=100, method='BP86', basis='def2-SVP'))
    assert orca.qc_params['software'] == 'orca'
    aimnet = resolve_reaction_request(files, dict(backend='aimnet_2', bias_max=100))
    assert aimnet.qc_params['method'] is None and 'model' not in aimnet.qc_params
    with pytest.raises(ValueError, match='native Gaussian'):
        resolve_reaction_request(files, dict(backend='gaussian', bias_max=100, method='B3LYP', basis='6-31G'))
    from dataclasses import replace
    from pyar.backend_capabilities import BACKEND_CAPABILITIES
    with patch.dict(BACKEND_CAPABILITIES, xtb=replace(BACKEND_CAPABILITIES['xtb'], native_optimization=False)):
        with pytest.raises(ValueError, match='native'):
            resolve(files)


def test_check_no_mutations(tmp_path, workflow, monkeypatch, capsys):
    files = inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    main(['react', *files, '--backend', 'xtb', '--bias-max', '100', '--check'])
    assert not workflow.called
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    assert 'adaptive' in capsys.readouterr().out


@pytest.mark.parametrize('failure', [ValueError("Missing executable 'xtb'; official installation"),
                                     ImportError('geometric missing; pip install "pyar-chem[xtb]"'),
                                     ValueError('native relaxation unavailable')])
def test_preflight_failure_no_work(tmp_path, monkeypatch, workflow, failure, capsys):
    monkeypatch.chdir(tmp_path)
    with patch('pyar.scripts.modern_react.preflight_reaction', side_effect=failure):
        with pytest.raises(SystemExit) as exc:
            main(['react', *inputs(tmp_path), '--backend', 'xtb', '--bias-max', '100', '--check'])
    assert exc.value.code == 2 and not workflow.called
    assert not (tmp_path/'reaction').exists()
    assert 'No calculations were started' in capsys.readouterr().err


def test_two_phase_preflight_and_identity(tmp_path):
    request = resolve(inputs(tmp_path))
    with patch('pyar.reaction_request.preflight', return_value=['xtb']) as check, patch(
            'pyar.cli._workflow_requirement_messages', return_value=[]), patch(
            'pyar.reaction_request.import_module'):
        requirements = preflight_reaction(request)
    assert check.call_count == 2
    biased, native = [call.args[0] for call in check.call_args_list]
    assert biased['geometry_optimizer'] == 'geometric'
    assert native['geometry_optimizer'] == 'native' and native['gamma'] == 0
    assert 'bias_controller' not in native and native['trace_enabled'] is False
    assert {'obabel', 'native unbiased relaxation'} <= set(requirements)


def test_real_missing_executable_requirement(tmp_path):
    request = resolve(inputs(tmp_path))
    import shutil
    actual_which = shutil.which
    with patch('pyar.cli.shutil.which', side_effect=lambda name: None if name == 'xtb' else actual_which(name)):
        with pytest.raises(ValueError, match='xtb-docs'):
            preflight_reaction(request)


def test_restart_compatible_readonly_and_mismatch(tmp_path):
    request = resolve(inputs(tmp_path))
    state = ReactionRunState.create(tmp_path, request.restart_request, [], request.reactants)
    before = state.state_file.read_bytes()
    assert validate_restart(request, tmp_path) is not None
    assert state.state_file.read_bytes() == before
    changed = resolve(request.reactants and [m.relative_path for m in request.reactants], orientations=20)
    with pytest.raises(ReactionStateError):
        validate_restart(changed, tmp_path)
    assert state.state_file.read_bytes() == before


def test_restart_check_and_resume_call(tmp_path, monkeypatch, workflow):
    files = inputs(tmp_path)
    request = resolve(files)
    state = ReactionRunState.create(tmp_path, request.restart_request, [], request.reactants)
    before = state.state_file.read_bytes()
    monkeypatch.chdir(tmp_path)
    main(['react', *files, '--backend', 'xtb', '--bias-max', '100', '--check'])
    assert not workflow.called and state.state_file.read_bytes() == before
    main(['react', *files, '--backend', 'xtb', '--bias-max', '100'])
    assert workflow.called


def test_help_and_legacy_controller_default(capsys):
    from pyar.biases.controller import resolve_controller_policy
    assert resolve_controller_policy(None) == 'fixed'
    with pytest.raises(SystemExit) as exc:
        main(['react', '--help'])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    assert '--bias-max' in output and '--backend' in output and 'adaptive' in output


def test_missing_geometric_has_pyar_extra_guidance(tmp_path):
    request = resolve(inputs(tmp_path))
    import importlib.util
    find = importlib.util.find_spec
    with patch('pyar.cli.importlib.util.find_spec', side_effect=lambda name: None if name == 'geometric' else find(name)):
        with pytest.raises(ValueError, match=r'pyar-chem\[xtb\]'):
            preflight_reaction(request)


def test_missing_identity_executable(tmp_path):
    request = resolve(inputs(tmp_path))
    with patch('pyar.reaction_request.preflight', return_value=[]), patch(
            'pyar.cli.shutil.which', return_value=None):
        with pytest.raises(ValueError, match='openbabel.org'):
            preflight_reaction(request)


def test_import_failure_prevents_execution(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with patch('pyar.reaction_request.preflight', return_value=[]), patch(
            'pyar.cli._workflow_requirement_messages', return_value=[]), patch(
            'pyar.reaction_request.import_module', side_effect=ImportError('provider unavailable')), patch(
            'pyar.workflows.reaction.react') as run:
        with pytest.raises(SystemExit):
            main(['react', *inputs(tmp_path), '--backend', 'xtb', '--bias-max', '100'])
    assert not run.called and not (tmp_path/'reaction').exists()
    assert 'provider unavailable' in capsys.readouterr().err


def test_missing_restart_snapshot_readonly(tmp_path):
    request = resolve(inputs(tmp_path))
    state = ReactionRunState.create(tmp_path, request.restart_request, [request.reactants[0]], request.reactants)
    reference = state.data['pending_orientations'][0]['path']
    (tmp_path/'reaction'/reference).unlink()
    before = state.state_file.read_bytes()
    with pytest.raises(ReactionStateError, match='snapshot'):
        validate_restart(request, tmp_path)
    assert state.state_file.read_bytes() == before


@pytest.mark.parametrize('options', [dict(opt_threshold='invalid'), dict(bias_scheduled_alpha=.01),
                                    dict(bias_controller='fixed', bias_min=200), dict(charge=[0, 0, 0])])
def test_review_invalid_options(tmp_path, options):
    with pytest.raises(ValueError):
        resolve(inputs(tmp_path), **options)


def test_invalid_xyz_before_directory_creation(tmp_path, workflow, monkeypatch):
    files = inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    Path(files[1]).write_text('2\nmalformed\nH 0 0 0\n')
    with pytest.raises(SystemExit):
        main(['react', *files, '--backend', 'xtb', '--bias-max', '100'])
    assert not workflow.called and not (tmp_path/'reaction').exists()


def test_legacy_standalone_contract_unchanged(tmp_path, monkeypatch):
    from importlib import import_module
    import sys
    monkeypatch.chdir(tmp_path)
    legacy = import_module('pyar.scripts.react')
    files = inputs(tmp_path)
    with patch.object(sys, 'argv', ['pyar-react', *files, '--software', 'xtb', '--bias-max', '100']):
        with pytest.raises(SystemExit):
            legacy.argument_parse()  # -N is still required.
    with patch.object(sys, 'argv', ['pyar-react', *files, '-N', '8', '--software', 'xtb', '--bias-max', '100']):
        args = legacy.argument_parse()
        assert args.bias_controller is None
        with pytest.raises(SystemExit, match='bias-min'):
            legacy.main()  # No controller means legacy fixed, still requiring min.
    with patch.object(sys, 'argv', ['pyar-react', *files, '-N', '8', '--software', 'xtb', '--bias-min', '10', '--bias-max', '100']), patch(
            'pyar.workflows.reaction.react') as run:
        legacy.main()
    assert run.call_args.args[2:5] == (10, 100, 8)
    assert 'bias_controller' not in run.call_args.args[5]
