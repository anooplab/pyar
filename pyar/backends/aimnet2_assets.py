"""AIMNet2 runtime asset paths and validation."""

from importlib import resources
from pathlib import Path


MODEL_PATH = str(
    resources.files("pyar").joinpath("AIMNet2/models/aimnet2_wb97m-d3_0.jpt")
)
SCRIPT_PATH = str(
    resources.files("pyar").joinpath("AIMNet2/calculators/aimnet2_ase_opt.py")
)


def missing_aimnet2_assets(include_script=True):
    """Return missing AIMNet2 asset paths required by PyAR wrappers."""
    paths = [MODEL_PATH]
    if include_script:
        paths.append(SCRIPT_PATH)
    return [path for path in paths if not Path(path).is_file()]


def aimnet2_asset_error(missing):
    """Build the standard missing-asset error message."""
    return (
        "AIMNet2 model assets are not bundled in the pyar-chem wheel. "
        "Install or configure the model files separately, for example from the "
        "upstream project, a separate model package, or a future model-download "
        "command. Missing files: "
        + ", ".join(missing)
    )


def validate_aimnet2_runtime_assets(include_script=True):
    """Raise a clear error if required AIMNet2 runtime assets are absent."""
    missing = missing_aimnet2_assets(include_script=include_script)
    if missing:
        raise FileNotFoundError(aimnet2_asset_error(missing))


__all__ = [
    "MODEL_PATH",
    "SCRIPT_PATH",
    "aimnet2_asset_error",
    "missing_aimnet2_assets",
    "validate_aimnet2_runtime_assets",
]
