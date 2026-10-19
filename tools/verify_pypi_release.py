"""Verify a published PyAR version and exercise a clean PyPI installation."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import urllib.error
import urllib.request
import venv


class RegistryError(ValueError):
    """Raised when PyPI metadata or files disagree with the approved release."""


def validate_pypi_json(document: dict, version: str, expected: dict[str, str]) -> dict:
    info = document.get("info", {})
    if info.get("version") != version:
        raise RegistryError(
            f"PyPI reports latest version {info.get('version')!r}, expected {version!r}"
        )
    files = document.get("releases", {}).get(version)
    if not isinstance(files, list):
        raise RegistryError(f"PyPI has no release file list for {version}")
    actual = {
        item.get("filename"): item.get("digests", {}).get("sha256") for item in files
    }
    if len(actual) != len(files):
        raise RegistryError(f"PyPI contains duplicate filenames for {version}")
    if actual != expected:
        raise RegistryError(f"PyPI file names or SHA-256 hashes differ: {actual!r}")
    return actual


def fetch_pypi(version: str, timeout: float = 15.0) -> dict:
    url = f"https://pypi.org/pypi/pyar-chem/{version}/json"
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def verify_registry(
    version: str,
    expected: dict[str, str],
    *,
    attempts: int = 40,
    interval: float = 15.0,
) -> dict:
    last_error = None
    for attempt in range(attempts):
        try:
            document = fetch_pypi(version)
            return validate_pypi_json(document, version, expected)
        except (urllib.error.HTTPError, urllib.error.URLError, RegistryError) as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(interval)
    raise RegistryError(
        f"PyPI verification failed after {attempts} attempts: {last_error}"
    )


def install_and_smoke(version: str, *, timeout: int = 600) -> dict:
    with tempfile.TemporaryDirectory(prefix="pyar-pypi-smoke-") as temporary:
        root = Path(temporary)
        environment = root / "venv"
        venv.EnvBuilder(with_pip=True).create(environment)
        python = environment / "bin" / "python"
        clean_env = os.environ.copy()
        clean_env.pop("PYTHONPATH", None)
        clean_env.pop("PYTHONHOME", None)
        clean_env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
        subprocess.run(
            [str(python), "-m", "pip", "install", f"pyar-chem=={version}"],
            cwd=root,
            env=clean_env,
            check=True,
            timeout=timeout,
        )
        subprocess.run(
            [str(python), "-m", "pip", "check"],
            cwd=root,
            env=clean_env,
            check=True,
            timeout=timeout,
        )
        code = (
            "import importlib.metadata as m, pyar; "
            f"assert m.version('pyar-chem') == {version!r}; "
            "print(pyar.__file__)"
        )
        imported = subprocess.run(
            [str(python), "-c", code],
            cwd=root,
            env=clean_env,
            text=True,
            capture_output=True,
            check=True,
            timeout=timeout,
        ).stdout.strip()
        if not imported.startswith(str(environment)):
            raise RegistryError(
                f"installed import resolved outside clean venv: {imported}"
            )
        for command in ("pyar-cli", "pyar"):
            subprocess.run(
                [str(environment / "bin" / command), "--help"],
                cwd=root,
                env=clean_env,
                check=True,
                timeout=timeout,
                stdout=subprocess.DEVNULL,
            )
        return {
            "installed_version": version,
            "import_path": imported,
            "cli_smokes": ["pyar-cli --help", "pyar --help"],
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("approval", type=Path)
    args = parser.parse_args()
    approval = json.loads(args.approval.read_text(encoding="utf-8"))
    version = approval["version"]
    hashes = verify_registry(version, approval["artifacts"])
    smoke = install_and_smoke(version)
    print(
        json.dumps(
            {"status": "PASS", "version": version, "pypi_sha256": hashes, **smoke},
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
