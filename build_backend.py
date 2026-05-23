"""Minimal PEP 517 backend for PyAR.

This backend keeps installation self-contained and offline-friendly while
preserving editable installs via a generated .pth file.
"""

from __future__ import annotations

import base64
import csv
import hashlib
import io
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
import zipfile


NAME = "pyar"
VERSION = "1.1.0"
SUMMARY = "A Python Code for Aggregation and Reaction"
HOME_PAGE = "https://github.com/anooplab/pyar"
AUTHOR = "Anoop et al"
AUTHOR_EMAIL = "anoop@chem.iitkgp.ac.in"
LICENSE = "GPL v3"
REQUIRES_PYTHON = ">=3.6"
INSTALL_REQUIRES = [
    "numpy>=1.18.4",
    "autograd>=1.3",
    "ase",
    "torch",
    "torchani==2.0",
    "MDAnalysis",
    "pandas>=1.0.5",
    "scipy>=1.5.2",
    "scikit-learn>=0.23.2",
    "dscribe",
    "pyh5md",
    "h5py",
    "networkx",
    "matplotlib",
    "hdbscan",
    "DBCV @ git+https://github.com/christopherjenness/DBCV.git",
    "openbabel-wheel",
]
ENTRY_POINTS = [
    "pyar-cli=pyar.cli:main",
    "pyar-react=pyar.scripts.react:main",
    "pyar-explore=pyar.scripts.explore:main",
    "pyar-optimiser=pyar.scripts.optimiser:main",
    "pyar-tabu=pyar.scripts.tabu:main",
    "pyar-clustering=pyar.scripts.clustering:main",
    "pyar-benchmark-clustering=pyar.scripts.benchmark_clustering:main",
    "pyar-similarity=pyar.scripts.similarity:main",
    "pyar-descriptor=pyar.scripts.descriptor:main",
    "pyar-mlopt=pyar.interface.mlopt:main",
    "pyar-aimnet2-ase-opt=pyar.AIMNet2.calculators.aimnet2_ase_opt:main",
]

ROOT = Path(__file__).resolve().parent
DIST_INFO = f"{NAME}-{VERSION}.dist-info"
WHEEL_NAME = f"{NAME}-{VERSION}-py3-none-any.whl"


def _metadata_text() -> str:
    lines = [
        "Metadata-Version: 2.1",
        f"Name: {NAME}",
        f"Version: {VERSION}",
        f"Summary: {SUMMARY}",
        f"Home-page: {HOME_PAGE}",
        f"Author: {AUTHOR}",
        f"Author-email: {AUTHOR_EMAIL}",
        f"License: {LICENSE}",
        f"Requires-Python: {REQUIRES_PYTHON}",
    ]
    for requirement in INSTALL_REQUIRES:
        lines.append(f"Requires-Dist: {requirement}")
    return "\n".join(lines) + "\n"


def _wheel_text() -> str:
    return "\n".join(
        [
            "Wheel-Version: 1.0",
            "Generator: build_backend",
            "Root-Is-Purelib: true",
            "Tag: py3-none-any",
        ]
    ) + "\n"


def _entry_points_text() -> str:
    return "[console_scripts]\n" + "\n".join(ENTRY_POINTS) + "\n"


def _top_level_text() -> str:
    return "pyar\n"


def _iter_package_files() -> Iterable[Path]:
    for path in ROOT.joinpath("pyar").rglob("*"):
        if path.is_dir():
            continue
        if "__pycache__" in path.parts:
            continue
        if path.suffix in {".pyc", ".pyo"}:
            continue
        yield path


def _record_hash(data: bytes) -> str:
    digest = hashlib.sha256(data).digest()
    encoded = base64.urlsafe_b64encode(digest).decode("ascii").rstrip("=")
    return f"sha256={encoded}"


def _normalize_archive_path(path: Path) -> str:
    return path.as_posix()


def _build_file_map(editable: bool) -> Dict[str, bytes]:
    files: Dict[str, bytes] = {}
    if editable:
        files["pyar.pth"] = (str(ROOT) + "\n").encode("utf-8")
        return files

    for source in _iter_package_files():
        archive_name = _normalize_archive_path(source.relative_to(ROOT))
        files[archive_name] = source.read_bytes()

    return files


def _dist_info_files() -> Dict[str, bytes]:
    return {
        f"{DIST_INFO}/METADATA": _metadata_text().encode("utf-8"),
        f"{DIST_INFO}/WHEEL": _wheel_text().encode("utf-8"),
        f"{DIST_INFO}/entry_points.txt": _entry_points_text().encode("utf-8"),
        f"{DIST_INFO}/top_level.txt": _top_level_text().encode("utf-8"),
    }


def _write_wheel(archive_path: Path, files: Dict[str, bytes]) -> None:
    records: List[Tuple[str, str, str]] = []
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for archive_name, data in sorted(files.items()):
            info = zipfile.ZipInfo(archive_name)
            source = ROOT / archive_name
            mode = source.stat().st_mode if source.exists() else 0o100644
            info.external_attr = (mode & 0xFFFF) << 16
            zf.writestr(info, data)
            records.append((archive_name, _record_hash(data), str(len(data))))

        record_name = f"{DIST_INFO}/RECORD"
        record_rows = records + [(record_name, "", "")]
        record_buffer = io.StringIO()
        writer = csv.writer(record_buffer, lineterminator="\n")
        for row in record_rows:
            writer.writerow(row)
        zf.writestr(record_name, record_buffer.getvalue().encode("utf-8"))


def _metadata_directory(metadata_directory: str) -> Path:
    path = Path(metadata_directory) / DIST_INFO
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_metadata_tree(target: Path) -> None:
    (target / "METADATA").write_text(_metadata_text(), encoding="utf-8")
    (target / "WHEEL").write_text(_wheel_text(), encoding="utf-8")
    (target / "entry_points.txt").write_text(_entry_points_text(), encoding="utf-8")
    (target / "top_level.txt").write_text(_top_level_text(), encoding="utf-8")


def get_requires_for_build_wheel(config_settings=None):
    return []


def get_requires_for_build_editable(config_settings=None):
    return []


def prepare_metadata_for_build_wheel(metadata_directory, config_settings=None):
    target = _metadata_directory(metadata_directory)
    _write_metadata_tree(target)
    return DIST_INFO


def prepare_metadata_for_build_editable(metadata_directory, config_settings=None):
    target = _metadata_directory(metadata_directory)
    _write_metadata_tree(target)
    return DIST_INFO


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    wheel_path = Path(wheel_directory) / WHEEL_NAME
    files = _build_file_map(editable=False)
    files.update(_dist_info_files())
    _write_wheel(wheel_path, files)
    return WHEEL_NAME


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    wheel_path = Path(wheel_directory) / WHEEL_NAME
    files = _build_file_map(editable=True)
    files.update(_dist_info_files())
    _write_wheel(wheel_path, files)
    return WHEEL_NAME
