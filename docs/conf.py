"""Sphinx configuration for PyAR documentation."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
DOCS = Path(__file__).resolve().parent
PACKAGE_ROOT = ROOT / "pyar"
sys.path.insert(0, str(ROOT))

project = "PyAR"
author = "Anakuthil Anoop"
copyright = "2010, AnoopLab"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
]

autosummary_generate = True
autosectionlabel_prefix_document = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
}

autodoc_mock_imports = [
    "ase",
    "autograd",
    "dscribe",
    "h5py",
    "hdbscan",
    "matplotlib",
    "MDAnalysis",
    "networkx",
    "numpy",
    "mlatom",
    "pandas",
    "pyh5md",
    "scipy",
    "sklearn",
    "torch",
    "torchani",
]

API_EXCLUDED_MODULES = {
    "pyar",
    "pyar.tests",
}

API_EXCLUDED_PREFIXES = (
    "pyar.tests.",
    "pyar.AIMNet2.",
    "pyar.backends.",
    "pyar.checkvalidity",
    "pyar.interface.",
    "pyar.mlatom.",
)

templates_path = ["_templates"]
exclude_patterns = ["_build", "_autosummary"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

from pyar import __version__  # noqa: E402

version = __version__
release = __version__


def _module_name_from_path(path):
    """Return the importable PyAR module name for a Python source file."""
    relative = path.relative_to(ROOT).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _discover_pyar_modules():
    """Discover PyAR modules for generated API documentation."""
    modules = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        if any(part.startswith(".") for part in path.relative_to(ROOT).parts):
            continue
        module_name = _module_name_from_path(path)
        if not module_name or module_name in API_EXCLUDED_MODULES:
            continue
        if "-" in module_name:
            continue
        if any(module_name.startswith(prefix) for prefix in API_EXCLUDED_PREFIXES):
            continue
        modules.append(module_name)
    return modules


def _write_generated_api(_app):
    """Create a docstring-driven API page every time Sphinx builds.

    Chemistry-facing pages are hand-written. This generated page is for
    developers and follows the current Python package layout automatically, so
    newly added or refactored modules appear in the documentation when their
    docstrings are present and the docs are rebuilt.
    """
    modules = _discover_pyar_modules()
    output = DOCS / "generated_api.rst"
    lines = [
        "Generated API from Docstrings",
        "==============================",
        "",
        ".. note::",
        "",
        "   This page is generated during the Sphinx build by scanning the ``pyar``",
        "   package. Add or update Python docstrings when adding, moving, or",
        "   refactoring modules; the developer API documentation will then update",
        "   automatically on the next documentation build.",
        "",
    ]
    current_group = None
    for module in modules:
        parts = module.split(".")
        group = parts[1] if len(parts) > 1 else "top-level modules"
        if group != current_group:
            current_group = group
            title = group.replace("_", " ").title()
            lines.extend([title, "-" * len(title), ""])
        lines.extend([
            module,
            "~" * len(module),
            "",
            f".. automodule:: {module}",
            "   :no-index:",
            "   :members:",
            "   :undoc-members:",
            "   :show-inheritance:",
            "",
        ])
    output.write_text("\n".join(lines), encoding="utf-8")


def setup(app):
    """Register build-time documentation generation hooks."""
    app.connect("builder-inited", _write_generated_api)
