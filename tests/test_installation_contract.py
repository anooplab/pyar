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


def test_backend_python_requirements_are_installed_by_their_capability_extra():
    from pyar.backend_capabilities import BACKEND_CAPABILITIES

    for name, capability in BACKEND_CAPABILITIES.items():
        for module in capability.required_python_modules:
            requirement = dependencies.PYTHON_REQUIREMENTS[module]
            assert requirement.package.lower() in resolved_extra(capability.optional_extra_hint), name


def test_backend_executables_have_authoritative_installation_metadata():
    from pyar.backend_capabilities import BACKEND_CAPABILITIES

    for capability in BACKEND_CAPABILITIES.values():
        for executable in capability.required_executables:
            requirement = dependencies.EXTERNAL_PROGRAMS[executable]
            assert requirement.installation_url.startswith('https://')
            assert not requirement.pip_installable


def test_aiqm1_preflights_gaussian_before_starting_its_subprocess():
    from pyar.backend_capabilities import get_backend_capabilities

    assert 'g16' in get_backend_capabilities('aiqm1_mlatom').required_executables
    assert {'g16', 'xtb'} <= get_backend_capabilities('xtb-aiqm1').required_executables


def test_doctor_optional_absence_is_informational(capsys):
    from pyar.scripts import cli_meta

    with patch.object(cli_meta, '_module_available', side_effect=lambda name: name in dependencies.CORE_MODULES), \
         patch.object(cli_meta.shutil, 'which', return_value=None):
        cli_meta.doctor(['--all-backends'])
    output = capsys.readouterr().out
    assert 'healthy' in output
    assert 'MISSING' in output
    assert 'only workflows' in output


def test_broken_workflow_dependency_is_reported_before_execution():
    with patch.object(dependencies, 'module_available', return_value=True), \
         patch.object(dependencies, 'import_module', side_effect=ImportError('missing shared dependency')):
        with pytest.raises(ImportError, match=r'pyar-chem\[aggregate\]') as caught:
            dependencies.require_workflow('aggregate')
    assert 'missing shared dependency' in str(caught.value)


@pytest.mark.parametrize('source', ['input.xyz', 'input.sdf'])
def test_direct_aimnet_optimizer_missing_torch_is_an_actionable_exit(source):
    import builtins
    from pyar.AIMNet2.calculators import aimnet2_ase_opt

    original_import = builtins.__import__

    def without_torch(name, *args, **kwargs):
        if name == 'torch':
            raise ModuleNotFoundError("No module named 'torch'", name='torch')
        return original_import(name, *args, **kwargs)

    with patch('sys.argv', ['pyar-aimnet2-ase-opt', '--charge', '0',
                            'model.jpt', source, 'output.xyz']), \
         patch('builtins.__import__', side_effect=without_torch):
        with pytest.raises(SystemExit, match=r'pyar-chem\[aimnet2\]'):
            aimnet2_ase_opt.main()
