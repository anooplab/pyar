"""Saved settings preserve argparse action semantics and exact replay tokens."""

import argparse

import pytest

from pyar.run_config import arguments_from_values


@pytest.mark.parametrize(('settings', 'expected'), [
    ({'enabled': True}, ['--enabled']),
    ({'enabled': False}, ['--no-enabled']),
    ({'force': True}, ['--force']),
    ({'force': False}, []),
    ({'cache': False}, ['--disable-cache']),
    ({'cache': True}, []),
    ({'verbosity_count': 3}, ['-v', '-v', '-v']),
    ({'verbosity_count': 0}, []),
    ({'label': ['alpha', '-beta']}, ['--label=alpha', '--label=-beta']),
    ({'pair': [['H', 'He'], ['C', 'N']]}, ['--pair', 'H', 'He', '--pair', 'C', 'N']),
    ({'atoms': ['H', 'He']}, ['--atoms', 'H', 'He']),
    ({'atoms': []}, []),
    ({'optional_atoms': []}, []),
    ({'optional_atoms': ['H']}, ['--optional-atoms', 'H']),
    ({'fixed': [1, 2]}, ['--fixed', '1', '2']),
    ({'keyword': '-unrestricted'}, ['--keyword=-unrestricted']),
    ({'short': 'value'}, ['-s', 'value']),
    ({'keyword': None}, []),
    ({'inputs': ['--backend', 'space name.xyz']}, ['--', '--backend', 'space name.xyz']),
])
def test_saved_action_tokens_and_parser_round_trip(settings, expected):
    parser = action_parser()
    tokens = arguments_from_values(settings, parser)
    assert tokens == expected
    parsed = vars(parser.parse_args(tokens))
    for key, value in settings.items():
        if value is not None:
            assert parsed[key] == value


def action_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--enabled', action=argparse.BooleanOptionalAction)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--disable-cache', dest='cache', action='store_false')
    parser.add_argument('-v', dest='verbosity_count', action='count', default=0)
    parser.add_argument('--label', action='append')
    parser.add_argument('--pair', action='append', nargs=2)
    parser.add_argument('--atoms', nargs='+', default=[])
    parser.add_argument('--optional-atoms', nargs='*', default=[])
    parser.add_argument('--fixed', nargs=2, type=int)
    parser.add_argument('--keyword')
    parser.add_argument('-s', dest='short')
    parser.add_argument('inputs', nargs='*')
    return parser


@pytest.mark.parametrize(('settings', 'message'), [
    ({'enabled': 1}, '--enabled must be a Boolean'),
    ({'force': 'false'}, '--force must be a Boolean'),
    ({'cache': 0}, '--disable-cache must be a Boolean'),
    ({'verbosity_count': -1}, '-v must be a nonnegative integer'),
    ({'verbosity_count': 1.5}, '-v must be a nonnegative integer'),
    ({'label': 'H'}, '--label requires a list'),
    ({'pair': 'H'}, '--pair requires a list'),
    ({'atoms': 'H'}, '--atoms requires a list'),
    ({'fixed': 2}, '--fixed requires a list'),
    ({'keyword': ['a']}, '--keyword requires a scalar'),
    ({'keyword': {'method': 'HF'}}, '--keyword requires a scalar'),
    ({'keyword': False}, '--keyword requires a scalar'),
])
def test_invalid_saved_action_values_fail_before_replay(settings, message):
    with pytest.raises(ValueError, match=message):
        arguments_from_values(settings, action_parser())


def test_store_constant_rejects_unrepresentable_saved_default():
    parser = argparse.ArgumentParser()
    parser.add_argument('--force', action='store_true', default=None)
    with pytest.raises(ValueError, match='cannot represent the saved value False'):
        arguments_from_values({'force': False}, parser)


def test_replay_preserves_option_order_and_excludes_execution_permissions():
    parser = action_parser()
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--recover-stale-lock', action='store_true')
    values = {'inputs': ['-input.xyz'], 'keyword': '-x', 'force': True,
              'recover_stale_lock': True, 'check': True}
    assert arguments_from_values(values, parser) == ['--force', '--keyword=-x', '--', '-input.xyz']
