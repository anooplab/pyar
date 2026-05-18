"""Sphinx configuration for PyAR documentation."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

project = "PyAR"
author = "Anakuthil Anoop"
copyright = "2010, AnoopLab"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
]

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
    "pandas",
    "pyh5md",
    "scipy",
    "sklearn",
    "torch",
    "torchani",
]

templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

from pyar import __version__  # noqa: E402

version = __version__
release = __version__
