"""The installation inventory is read-only and works with a minimal core."""

from importlib.metadata import PackageNotFoundError
import json
from unittest.mock import patch

import pytest

from pyar import modern_cli
from pyar.dependencies import CORE_MODULES
from pyar.scripts import installation


@pytest.mark.parametrize('python_version', [(3, 10), (3, 13)])
def test_core_only_status_lists_missing_optional_requirements(tmp_path, monkeypatch, capsys, python_version):
    monkeypatch.chdir(tmp_path)
    available = set(CORE_MODULES) | {'tomli'}

    def installed_version(package):
        if package in available:
            return '1.0'
        raise PackageNotFoundError(package)

    with patch.object(installation.sys, 'version_info', python_version), \
         patch.object(installation, 'version', side_effect=installed_version), \
         patch.object(installation, 'module_available', side_effect=lambda name: name in available), \
         patch.object(installation.shutil, 'which', return_value=None), \
         patch('subprocess.run', side_effect=AssertionError('status must not launch programs')):
        modern_cli.main(['install', 'status', '--json'])
    result = json.loads(capsys.readouterr().out)
    assert result['ok']
    assert result['executable'] == installation.sys.executable
    assert result['environment'] == installation.sys.prefix
    assert result['packages']['torch']['status'] == 'missing'
    assert result['executables']['orca']['status'] == 'missing'
    assert not result['workflows']['aggregate']['available']
    assert result['packages']['numpy']['version'] == '1.0'
    assert ('tomli' in result['packages']) == (python_version < (3, 11))
    assert list(tmp_path.iterdir()) == []


def test_status_reports_installed_but_unavailable_package(capsys):
    with patch.object(installation, 'version', return_value='1.0'), \
         patch.object(installation, 'module_available', side_effect=lambda name: name != 'torch'):
        modern_cli.main(['install', 'status'])
    output = capsys.readouterr().out
    assert 'Interpreter:' in output and 'Environment:' in output
    assert 'BROKEN     torch' in output
    assert 'AIMNet2 model:' in output


def test_missing_core_requirement_makes_status_fail(capsys):
    with patch.object(installation, 'module_available', side_effect=lambda name: name != 'scipy'):
        with pytest.raises(SystemExit) as error:
            modern_cli.main(['install', 'status', '--json'])
    assert error.value.code == 1
    assert not json.loads(capsys.readouterr().out)['ok']


def test_only_status_is_supported():
    with patch.object(installation, 'installation_status') as inventory:
        with pytest.raises(SystemExit) as error:
            modern_cli.main(['install', 'geometric'])
    assert error.value.code == 2
    inventory.assert_not_called()
