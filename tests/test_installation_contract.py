"""User installation contracts independent of the developer environment."""

import importlib.util
from pathlib import Path
from unittest.mock import patch

import pytest
from packaging.requirements import Requirement

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from pyar import dependencies


def resolved_extra(extra):
    metadata = tomllib.loads(Path('pyproject.toml').read_text())
    extras = metadata['project']['optional-dependencies']
    result = set()

    def collect(name):
        for text in extras[name]:
            requirement = Requirement(text)
            if requirement.name == 'pyar-chem':
                for nested in requirement.extras:
                    collect(nested)
            else:
                result.add(requirement.name.lower())

    collect(extra)
    return result


@pytest.mark.parametrize('workflow', ['aggregate', 'grow', 'microsolvate', 'solvate', 'conformer', 'reaction', 'neb'])
def test_workflow_installs_every_registered_python_requirement(workflow):
    requirements = resolved_extra(workflow)
    for module in dependencies.WORKFLOW_MODULES[workflow]:
        assert dependencies.PYTHON_REQUIREMENTS[module].package.lower() in requirements


def test_aggregate_does_not_install_unrelated_stacks():
    assert not resolved_extra('aggregate') & {'torch', 'mlatom', 'rdkit', 'geometric', 'sella', 'pytest', 'sphinx'}


def test_conformer_aimnet_composition_avoids_reaction_stack():
    requirements = resolved_extra('conformer') | resolved_extra('aimnet2')
    assert {'rdkit', 'ase', 'torch'} <= requirements
    assert not requirements & {'geometric', 'dscribe', 'hdbscan', 'mlatom'}


def test_all_excludes_developer_tools_and_fragile_openbabel_requirement():
    assert not resolved_extra('all') & {'pytest', 'pytest-cov', 'build', 'twine', 'sphinx', 'sphinx-rtd-theme', 'openbabel'}


def test_missing_workflow_stack_gives_one_complete_install_command():
    with patch.object(dependencies, 'module_available', return_value=False):
        with pytest.raises(ImportError, match=r'pyar-chem\[reaction\]') as caught:
            dependencies.require_workflow('reaction')
    assert "'dscribe'" in str(caught.value)
    assert "'geometric'" in str(caught.value)


def test_missing_orca_explains_separate_installation_and_discovery():
    message = dependencies.missing_executable_message('orca')
    assert 'ORCA' in message and 'PATH' in message
    assert 'https://www.faccts.de/' in message
    assert 'pip install' not in message


def test_doctor_optional_absence_is_informational(capsys):
    from pyar.scripts import cli_meta

    with patch.object(cli_meta, '_module_available', side_effect=lambda name: name in dependencies.CORE_MODULES), \
         patch.object(cli_meta.shutil, 'which', return_value=None):
        cli_meta.doctor(['--all-backends'])
    output = capsys.readouterr().out
    assert 'healthy' in output
    assert 'MISSING' in output
    assert 'only workflows' in output
