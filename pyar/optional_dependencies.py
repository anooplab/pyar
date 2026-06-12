"""Helpers for optional dependency error messages."""

from __future__ import annotations


OPTIONAL_EXTRA_HINTS = {
    "ase": ("selection", "pyar-chem[selection]"),
    "dscribe": ("selection", "pyar-chem[selection]"),
    "hdbscan": ("selection", "pyar-chem[selection]"),
    "pandas": ("selection", "pyar-chem[selection]"),
    "sklearn": ("selection", "pyar-chem[selection]"),
    "MDAnalysis": ("selection", "pyar-chem[selection]"),
    "torch": ("aimnet2", "pyar-chem[aimnet2]"),
    "torchani": ("ml", "pyar-chem[ml]"),
    "mlatom": ("ml", "pyar-chem[ml]"),
    "pyh5md": ("ml", "pyar-chem[ml]"),
    "h5py": ("ml", "pyar-chem[ml]"),
    "geometric": ("xtb", "pyar-chem[xtb]"),
    "openbabel": ("openbabel", "pyar-chem[openbabel]"),
    "rdkit": ("conformer", "pyar-chem[conformer]"),
}


def optional_dependency_error(module_name, feature=None, extra=None):
    """Return an ImportError with a consistent optional-extra install hint."""
    default_extra, package_spec = OPTIONAL_EXTRA_HINTS.get(
        module_name,
        ("all", "pyar-chem[all]"),
    )
    extra_name = extra or default_extra
    install_spec = package_spec if extra is None else f"pyar-chem[{extra}]"
    feature_name = feature or module_name
    return ImportError(
        f"{feature_name} requires the `{extra_name}` extra. "
        f"Install with `python -m pip install \"{install_spec}\"`."
    )


def require_optional(module_name, feature=None, extra=None):
    """Import an optional module or raise a clear extra-specific ImportError."""
    try:
        return __import__(module_name)
    except ImportError as exc:
        raise optional_dependency_error(module_name, feature=feature, extra=extra) from exc
