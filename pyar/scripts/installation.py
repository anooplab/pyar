"""Read-only inventory of the environment running PyAR."""

import argparse
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import shutil
import sys

from pyar import __version__
from pyar.dependencies import (
    CORE_MODULES, EXTERNAL_PROGRAMS, PYTHON_REQUIREMENTS, WORKFLOW_MODULES,
    module_available,
)


def installation_status():
    core = (*CORE_MODULES, 'tomli') if sys.version_info < (3, 11) else CORE_MODULES
    packages = {}
    for module, requirement in PYTHON_REQUIREMENTS.items():
        if module == 'tomli' and module not in core:
            continue
        try:
            installed_version = version(requirement.package)
        except PackageNotFoundError:
            installed_version = None
        available = module_available(module)
        status = ('installed' if installed_version else 'available') if available else (
            'broken' if installed_version else 'missing')
        packages[module] = dict(package=requirement.package, version=installed_version,
                                status=status, core=module in core,
                                install=requirement.install_hint)
    executables = {}
    for executable, requirement in EXTERNAL_PROGRAMS.items():
        path = shutil.which(executable)
        executables[executable] = dict(status='found' if path else 'missing', path=path,
                                      installation_url=requirement.installation_url,
                                      install=requirement.install_hint)
    from pyar.backends.aimnet2_assets import MODEL_PATH
    model = dict(path=MODEL_PATH, status='found' if Path(MODEL_PATH).is_file() else 'missing')
    workflows = {
        workflow: dict(available=all(packages[module]['status'] in {'installed', 'available'}
                                     for module in modules),
                       missing=[module for module in modules
                                if packages[module]['status'] not in {'installed', 'available'}])
        for workflow, modules in WORKFLOW_MODULES.items()
    }
    return dict(version=__version__, python=sys.version.split()[0],
                executable=sys.executable, environment=sys.prefix,
                ok=all(packages[module]['status'] in {'installed', 'available'} for module in core),
                packages=packages, workflows=workflows, executables=executables,
                resources={'aimnet2_model': model})


def main(argv=None, *, prog='pyar install'):
    parser = argparse.ArgumentParser(prog=prog, description='Inspect the current PyAR installation')
    commands = parser.add_subparsers(dest='action', required=True)
    status_parser = commands.add_parser('status', help='List installed and missing requirements')
    status_parser.add_argument('--json', action='store_true', help='Output a machine-readable inventory')
    args = parser.parse_args(argv)
    result = installation_status()
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"PyAR {result['version']} | Python {result['python']}")
        print(f"Interpreter: {result['executable']}\nEnvironment: {result['environment']}")
        for label, is_core in (('Core Python packages', True), ('Optional Python packages', False)):
            print(f'\n{label}:')
            for package in result['packages'].values():
                if package['core'] == is_core:
                    print(f"  {package['status'].upper():<10} {package['package']:<22} {package['version'] or '-'}")
        print('\nWorkflow Python requirements:')
        for workflow, state in result['workflows'].items():
            detail = 'available' if state['available'] else 'missing: ' + ', '.join(state['missing'])
            print(f'  {workflow:<14} {detail}')
        print('\nExternal programs on PATH:')
        for executable, state in result['executables'].items():
            print(f"  {state['status'].upper():<10} {executable:<22} {state['path'] or '-'}")
        print(f"\nAIMNet2 model: {result['resources']['aimnet2_model']['status']}")
        print(f"  {result['resources']['aimnet2_model']['path']}")
        print('\nStatus checks package metadata and discovery; use workflow --check to validate a run.')
        print('Missing optional requirements do not affect core installation health.')
    if not result['ok']:
        raise SystemExit(1)
