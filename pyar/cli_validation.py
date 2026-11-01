"""Shared electronic-state and backend input validation, independent of CLI parsers."""


import importlib.util




import shutil

import sys




from pyar.backend_errors import BackendExecutionError




from pyar.backend_capabilities import (
    get_backend_capabilities,
    normalize_backend_name,
)

QC_OPTION_ALIASES = {
    "basis": "--basis (basis set)",
    "method": "--method (functional/method)",
    "opt_threshold": "--opt-threshold",
    "opt_cycles": "--opt-cycles",
    "scf_threshold": "--scf-threshold",
    "scf_cycles": "--scf-cycles",
    "nprocs": "--nprocs",
    "custom_keywords": "--custom-keywords",
    "gamma": "--bias-min/--bias-max",
    "model": "--model",
}

QC_PARAMETER_KEYS = {
    "basis",
    "method",
    "opt_cycles",
    "opt_threshold",
    "scf_cycles",
    "scf_threshold",
    "nprocs",
    "gamma",
    "custom_keywords",
    "custom_keyword",
    "model",
}

def _normalize_parameter_list(values, default_value, number_of_inputs, label):
    """Expand a single provided value or validate a full per-input list."""
    resolved = values or [default_value for _ in range(number_of_inputs)]
    if len(resolved) == 1 and number_of_inputs > 1:
        return resolved * number_of_inputs
    if len(resolved) != number_of_inputs:
        sys.exit(f'{label} are not specified for all input files')
    return resolved

def _missing_module_message(module_name, install_hint):
    """Return a short install hint for a missing Python module."""
    if install_hint.startswith("http"):
        install_text = f"See {install_hint} for installation details."
    else:
        install_text = f"Install it with {install_hint}."
    return (
        f"Missing Python package '{module_name}'. "
        f"{install_text}"
    )

def _missing_executable_message(program, friendly_name, install_hint):
    """Return a short install hint for a missing executable."""
    if install_hint.startswith("http"):
        install_text = f"See {install_hint} for installation details."
    else:
        install_text = f"Install it with {install_hint}."
    return (
        f"Missing executable '{program}' for {friendly_name}. "
        f"{install_text} Make sure '{program}' is on PATH."
    )

from pyar.dependencies import EXTERNAL_PROGRAMS, PYTHON_REQUIREMENTS

_BACKEND_EXECUTABLE_REQUIREMENT_HINTS = {
    name: (requirement.name, requirement.installation_url)
    for name, requirement in EXTERNAL_PROGRAMS.items()
}

_BACKEND_MODULE_REQUIREMENT_HINTS = {
    name: requirement.install_hint for name, requirement in PYTHON_REQUIREMENTS.items()
}

def _workflow_requirement_messages(run_mode, software, geometry_optimizer):
    """Return user-facing messages for missing requirements on the active path."""
    from pyar.dependencies import (WORKFLOW_EXECUTABLES, workflow_requirement_messages,
                                   missing_executable_message)

    missing = workflow_requirement_messages(run_mode)
    capabilities = get_backend_capabilities(software)

    if geometry_optimizer == "geometric" and software is not None:
        if importlib.util.find_spec("ase") is None:
            missing.append(
                _missing_module_message(
                    "ase",
                    _BACKEND_MODULE_REQUIREMENT_HINTS["ase"],
                )
            )
        if importlib.util.find_spec("geometric") is None:
            missing.append(
                _missing_module_message(
                    "geometric",
                    _BACKEND_MODULE_REQUIREMENT_HINTS["geometric"],
                )
            )
        if shutil.which("geometric-optimize") is None:
            missing.append(
                _missing_executable_message(
                    "geometric-optimize",
                    "geomeTRIC",
                    "https://geometric.readthedocs.io/en/1.1/install.html",
                )
            )

    canonical_software = normalize_backend_name(software)
    executables = set(capabilities.required_executables) | set(WORKFLOW_EXECUTABLES.get(run_mode, ()))
    for executable in sorted(executables):
        if shutil.which(executable) is None:
            friendly_name, install_hint = _BACKEND_EXECUTABLE_REQUIREMENT_HINTS.get(
                executable,
                (canonical_software, 'python -m pip install "pyar-chem[all]"'),
            )
            missing.append(
                missing_executable_message(executable) if executable in EXTERNAL_PROGRAMS
                else _missing_executable_message(executable, friendly_name, install_hint)
            )

    for module_name in capabilities.required_python_modules:
        if importlib.util.find_spec(module_name) is None:
                missing.append(
                    _missing_module_message(
                        module_name,
                        _BACKEND_MODULE_REQUIREMENT_HINTS.get(module_name, "https://pypi.org/"),
                    )
                )

    return list(dict.fromkeys(missing))

def _preflight_cli_requirements(run_mode, software, geometry_optimizer):
    """Fail fast when the selected workflow path is missing requirements."""
    missing = _workflow_requirement_messages(run_mode, software, geometry_optimizer)
    if normalize_backend_name(software) in {"aimnet_2", "xtb-aimnet2"}:
        from pyar.backends.aimnet2_assets import validate_aimnet2_runtime_assets
        try:
            validate_aimnet2_runtime_assets(include_script=geometry_optimizer == "native")
        except FileNotFoundError as exc:
            missing.append(str(exc))
    if missing:
        raise SystemExit(
            "Missing requirements for the selected workflow:\n- " + "\n- ".join(missing)
        )

def _preflight_xtb_model(qc_params):
    """Reject a requested xTB Hamiltonian the installed executable cannot run."""
    if qc_params.get("software") != "xtb" or qc_params.get("xtb_model") != "gxtb":
        return
    from pyar.backends import require_executable
    from pyar.backends.xtb import xtb_supports_gxtb

    executable = require_executable("xtb", "xTB")
    if not xtb_supports_gxtb(executable):
        raise BackendExecutionError(
            f"xTB executable {executable!r} does not support --gxtb; "
            "select --xtb-model gfn2 or install a g-xTB-capable executable"
        )

def _infer_default_multiplicities(input_molecules, charges):
    """Infer spin multiplicity from electron parity when the CLI omits it."""
    multiplicities = []
    for mol, charge in zip(input_molecules, charges):
        n_electrons = sum(mol.atomic_number) - charge
        multiplicity = 1 if n_electrons % 2 == 0 else 2
        multiplicities.append(multiplicity)
        mol.multiplicity = multiplicity
    return multiplicities

def _validate_backend_spin_inputs(input_molecules):
    """Reject charge/multiplicity combinations with impossible electron parity."""
    for mol in input_molecules:
        n_electrons = sum(mol.atomic_number) - mol.charge
        if n_electrons % 2 == 0:
            if mol.multiplicity % 2 != 1:
                sys.exit(f"{n_electrons} (even) electrons and multiplicty {mol.multiplicity} (odd) is not pssible for {mol.name}")
        else:
            if mol.multiplicity % 2 == 1:
                sys.exit(f"{n_electrons} (odd) electrons and multiplicty {mol.multiplicity} (even) is not pssible for {mol.name}")

