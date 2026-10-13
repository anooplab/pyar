"""Installation metadata shared by workflow preflight and diagnostics.

Keep package requirements in pyproject.toml; tests check this registry against
that metadata. Probes are descriptive and are never launched during discovery.
"""

from dataclasses import dataclass
import importlib.util
import shutil


@dataclass(frozen=True)
class PythonRequirement:
    package: str
    extra: str

    @property
    def install_hint(self):
        spec = f"pyar-chem[{self.extra}]" if self.extra else "pyar-chem"
        return f'python -m pip install "{spec}"'


PYTHON_REQUIREMENTS = {
    "numpy": PythonRequirement("numpy", ""),
    "scipy": PythonRequirement("scipy", ""),
    "networkx": PythonRequirement("networkx", ""),
    "tomli": PythonRequirement("tomli", ""),
    "ase": PythonRequirement("ase", "geometric"),
    "dscribe": PythonRequirement("dscribe", "aggregate"),
    "hdbscan": PythonRequirement("hdbscan", "aggregate"),
    "sklearn": PythonRequirement("scikit-learn", "aggregate"),
    "rdkit": PythonRequirement("rdkit", "conformer"),
    "torch": PythonRequirement("torch", "aimnet2"),
    "torchani": PythonRequirement("torchani", "ml"),
    "mlatom": PythonRequirement("mlatom", "ml"),
    "geometric": PythonRequirement("geometric", "geometric"),
    "sella": PythonRequirement("sella", "sella"),
    "geodesic_interpolate": PythonRequirement("geodesic-interpolate", "geodesic"),
    "irmsd": PythonRequirement("irmsd", "structure-comparison"),
    "pandas": PythonRequirement("pandas", "analysis"),
    "MDAnalysis": PythonRequirement("MDAnalysis", "analysis"),
    "matplotlib": PythonRequirement("matplotlib", "analysis"),
}

CORE_MODULES = ("numpy", "scipy", "networkx")
_SELECTION = ("ase", "dscribe", "hdbscan", "sklearn")
WORKFLOW_MODULES = {
    "aggregate": _SELECTION,
    "grow": _SELECTION,
    "microsolvate": _SELECTION,
    "solvate": _SELECTION,
    "conformer": ("rdkit",),
    "reaction": (*_SELECTION, "geometric"),
    "react": (*_SELECTION, "geometric"),
    "neb": ("ase", "geometric"),
    "ts": ("ase", "geometric"),
    "irc": ("ase", "geometric"),
}


@dataclass(frozen=True)
class ExecutableRequirement:
    name: str
    installation_url: str
    install_hint: str
    version_probe: tuple[str, ...]
    pip_installable: bool = False


EXTERNAL_PROGRAMS = {
    "xtb": ExecutableRequirement("xTB", "https://xtb-docs.readthedocs.io/en/latest/setup.html",
        "mamba install -c conda-forge xtb", ("xtb", "--version")),
    "orca": ExecutableRequirement("ORCA", "https://www.faccts.de/docs/orca/6.1/tutorials/first_steps/install.html",
        "Download the platform installer from the official ORCA distribution (registration required)", ("orca", "--version")),
    "g16": ExecutableRequirement("Gaussian", "https://gaussian.com/installation/",
        "Install licensed Gaussian 16 using the vendor instructions", ("g16",)),
    "psi4": ExecutableRequirement("Psi4", "https://psicode.org/installs/",
        "Install Psi4 using its official platform instructions", ("psi4", "--version")),
    "mopac": ExecutableRequirement("MOPAC", "https://openmopac.net/",
        "Install the official MOPAC distribution", ("mopac",)),
}
for _executable in ("define", "ridft", "rdgrad", "statpt", "jobex"):
    EXTERNAL_PROGRAMS[_executable] = ExecutableRequirement("Turbomole",
        "https://www.turbomole.org/turbomole/turbomole-documentation/",
        "Install licensed Turbomole and configure its executable directory on PATH", (_executable,))
for _executable in ("obabel", "obminimize", "obenergy"):
    EXTERNAL_PROGRAMS[_executable] = ExecutableRequirement("OpenBabel",
        "https://openbabel.org/docs/Installation/install.html",
        "mamba install -c conda-forge openbabel", ("obabel", "-V"))
EXTERNAL_PROGRAMS["geometric-optimize"] = ExecutableRequirement("geomeTRIC",
    "https://geometric.readthedocs.io/en/latest/install.html",
    'python -m pip install "pyar-chem[geometric]"', ("geometric-optimize", "--help"), True)


def module_available(module):
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError, AttributeError):
        return False


def missing_python_message(module, *, extra=None):
    requirement = PYTHON_REQUIREMENTS.get(module)
    hint = (f'python -m pip install "pyar-chem[{extra}]"' if extra else
            requirement.install_hint if requirement else f'python -m pip install "{module}"')
    return f"Missing Python package '{module}'. Install it with {hint}."


def missing_executable_message(executable):
    requirement = EXTERNAL_PROGRAMS[executable]
    probe = " ".join(requirement.version_probe)
    return (f"Missing executable '{executable}' for {requirement.name}. "
            f"{requirement.install_hint}. Make sure '{executable}' is on PATH. "
            f"Verify discovery with: command -v {executable}. "
            f"Program check: {probe}. Official installation instructions: {requirement.installation_url}")


def workflow_requirement_messages(workflow):
    extra = "reaction" if workflow == "react" else "neb" if workflow in {"ts", "irc"} else workflow
    return [missing_python_message(module, extra=extra)
            for module in WORKFLOW_MODULES.get(workflow, ()) if not module_available(module)]


def require_workflow(workflow):
    """Validate a workflow's Python stack before execution creates artifacts."""
    messages = workflow_requirement_messages(workflow)
    if messages:
        raise ImportError("\n".join(messages))
    return list(WORKFLOW_MODULES.get(workflow, ()))
