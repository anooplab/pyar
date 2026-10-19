# encoding: utf-8
"""
Optimiser Module

Functions
---------

optimise(molecule, qc_params)
write_csv(csv_file, energy_dit)
bulk_optimize(input_files, qc_params)
"""

import json
import hashlib
import importlib.metadata
import logging
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from pyar import file_manager
from pyar.core.molecule import Molecule
from pyar.backend_capabilities import normalize_backend_name, validate_backend_capability
from pyar.backend_errors import BackendExecutionError
from pyar.io.xyz import XYZParseError, parse_xyz

optimiser_logger = logging.getLogger('pyar.optimiser')

SUCCESS_STATUSES = {True, 'converged'}
CYCLE_EXCEEDED_STATUSES = {'CycleExceeded', 'cycle_exceeded', 'cycleexceeded'}
RELEASE_CANDIDATE_STATUSES = {'ReleaseCandidate', 'release_candidate'}
_OPTIMIZATION_CACHE_SCHEMA = 1


def _backend_cache_provenance(qc_params):
    """Return a stable backend identity, or ``None`` when it cannot be verified.

    Cache reuse is deliberately disabled for unknown backend versions. A changed
    executable or model must not silently inherit results from an older run.
    """
    software = normalize_backend_name(qc_params.get("software", ""))
    provenance = {"backend": software, "method": qc_params.get("method")}
    if software == "xtb":
        executable = shutil.which("xtb")
        if not executable:
            return None
        try:
            completed = subprocess.run(
                [executable, "--version"], capture_output=True, text=True,
                timeout=5, check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        version = (completed.stdout + "\n" + completed.stderr).strip()
        if completed.returncode != 0 or not version:
            return None
        provenance.update(executable=str(Path(executable).resolve()), version=version)
    elif software in {"aimnet2", "xtb-aimnet2"}:
        try:
            from pyar.backends.aimnet2_assets import MODEL_PATH
            model = Path(MODEL_PATH)
            if not model.is_file():
                return None
            digest = hashlib.sha256(model.read_bytes()).hexdigest()
            torch_version = importlib.metadata.version("torch")
        except (ImportError, OSError, importlib.metadata.PackageNotFoundError):
            return None
        provenance.update(model_sha256=digest, torch_version=torch_version)
    else:
        # ORCA and other external engines have vendor-specific version probes;
        # do not guess or reuse caches until their version can be established.
        return None
    if qc_params.get("geometry_optimizer") == "geometric":
        try:
            provenance["geometric_version"] = importlib.metadata.version("geometric")
        except importlib.metadata.PackageNotFoundError:
            return None
    return _json_value(provenance)


def _json_value(value):
    """Convert request settings to stable JSON values or reject cacheability."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("nonfinite optimization setting")
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    raise TypeError(f"unsupported cache setting type: {type(value).__name__}")


def _request_digest(molecule, qc_params, backend_provenance):
    """Fingerprint input structure, electronic state, backend, and all settings."""
    payload = {
        "schema_version": _OPTIMIZATION_CACHE_SCHEMA,
        "symbols": list(molecule.atoms_list),
        "coordinates_angstrom": np.asarray(molecule.coordinates, dtype=float).tolist(),
        "charge": int(molecule.charge),
        "multiplicity": int(molecule.multiplicity),
        "scftype": str(molecule.scftype),
        "settings": _json_value(qc_params),
        "backend_provenance": backend_provenance,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _result_digest(path, *, geometric=False):
    """Hash parsed result content; geomeTRIC state owns the physical energy."""
    symbols, coordinates, _name, _title, energy = parse_xyz(path)
    value = {
        "symbols": symbols,
        "coordinates_angstrom": np.asarray(coordinates, dtype=float).tolist(),
    }
    if not geometric:
        value["energy_hartree"] = energy
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _archive_unverified_cache(result_path, metadata_path):
    """Move an unverified result aside before the backend writes a fresh one."""
    if not result_path.exists() and not metadata_path.exists():
        return
    archive = Path("stale_cache")
    archive.mkdir(exist_ok=True)
    suffix = 0
    while (archive / f"attempt_{suffix:03d}").exists():
        suffix += 1
    target = archive / f"attempt_{suffix:03d}"
    target.mkdir()
    for path in (result_path, metadata_path):
        if path.exists():
            shutil.move(str(path), str(target / path.name))


def _verified_cache(molecule, request_digest, backend_provenance, *, geometric=False):
    """Return a verified cached status/result, or archive and reject it."""
    result_path = Path(f"result_{molecule.name}.xyz")
    metadata_path = Path("optimization_cache.json")
    if not result_path.exists():
        return None
    try:
        symbols, coordinates, _name, _title, energy = parse_xyz(result_path)
    except (OSError, XYZParseError) as exc:
        raise BackendExecutionError(
            f"Malformed cached optimization result {result_path}: {exc}"
        ) from exc
    coordinates = np.asarray(coordinates, dtype=float).reshape((-1, 3))
    valid_geometry = (
        symbols == list(molecule.atoms_list)
        and coordinates.shape == (molecule.number_of_atoms, 3)
        and np.isfinite(coordinates).all()
        and energy is not None
        and math.isfinite(float(energy))
    )
    if not valid_geometry:
        raise BackendExecutionError(
            f"Cached optimization result {result_path} has invalid atoms, coordinates, or energy"
        )
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        status = metadata["status"]
        cached_energy = float(energy)
        valid = (
            metadata.get("schema_version") == _OPTIMIZATION_CACHE_SCHEMA
            and metadata.get("request_sha256") == request_digest
            and metadata.get("backend_provenance") == backend_provenance
            and metadata.get("result_sha256") == _result_digest(result_path, geometric=geometric)
            and status in {"converged", "cycle_exceeded", "release_candidate"}
        )
        if geometric:
            state = json.loads(Path("pyar_geometric_state.json").read_text(encoding="utf-8"))
            state_status = state.get("optimization_status")
            state_coordinates = np.asarray(state["positions_angstrom"], dtype=float)
            state_energy = float(state["backend_energy_hartree"])
            expected_state = status
            valid = valid and (
                state_status == expected_state
                and state_coordinates.shape == coordinates.shape
                and np.isfinite(state_coordinates).all()
                and np.allclose(state_coordinates, coordinates, rtol=0.0, atol=5.1e-6)
                and math.isfinite(state_energy)
                and np.asarray(state["forces_ev_per_angstrom"], dtype=float).shape == coordinates.shape
                and np.isfinite(np.asarray(state["forces_ev_per_angstrom"], dtype=float)).all()
            )
            if valid:
                cached_energy = state_energy
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
        valid = False
        status = None
    if not valid:
        _archive_unverified_cache(result_path, metadata_path)
        return None
    return status, cached_energy, coordinates


def _write_optimization_cache(molecule, request_digest, status, backend_provenance, *, geometric=False):
    """Atomically record a successful or explicitly provisional optimizer result."""
    result_path = Path(f"result_{molecule.name}.xyz")
    try:
        symbols, coordinates, _name, _title, energy = parse_xyz(result_path)
        coordinates = np.asarray(coordinates, dtype=float).reshape((-1, 3))
        if (symbols != list(molecule.atoms_list)
                or coordinates.shape != (molecule.number_of_atoms, 3)
                or not np.isfinite(coordinates).all()
                or energy is None or not math.isfinite(float(energy))):
            raise ValueError("optimizer wrote an invalid result geometry or energy")
        record = {
            "schema_version": _OPTIMIZATION_CACHE_SCHEMA,
            "request_sha256": request_digest,
            "backend_provenance": backend_provenance,
            "result_sha256": _result_digest(result_path, geometric=geometric),
        "status": (
            "cycle_exceeded" if is_cycle_exceeded(status)
            else "release_candidate" if status in RELEASE_CANDIDATE_STATUSES
            else "converged"
        ),
        }
        temporary = Path("optimization_cache.json.tmp")
        temporary.write_text(json.dumps(record, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace("optimization_cache.json")
    except (OSError, ValueError, TypeError, XYZParseError) as exc:
        raise BackendExecutionError(f"Could not validate optimization result for cache: {exc}") from exc


def is_success(status):
    return status in SUCCESS_STATUSES


def is_cycle_exceeded(status):
    return status in CYCLE_EXCEEDED_STATUSES


def is_usable(status):
    return is_success(status) or is_cycle_exceeded(status) or status in RELEASE_CANDIDATE_STATUSES


def _status_label(status):
    """Return a normalized human-readable optimization status label."""
    if is_success(status):
        return "success"
    if is_cycle_exceeded(status):
        return "cycle_exceeded_usable"
    if status in RELEASE_CANDIDATE_STATUSES:
        return "release_candidate_usable"
    if status is None:
        return "failed_none_status"
    return f"failed:{status}"


def apply_geometry_result(molecule, geometry):
    molecule.energy = geometry.energy
    if hasattr(geometry, 'optimized_coordinates'):
        molecule.optimized_coordinates = geometry.optimized_coordinates
        molecule.coordinates = geometry.optimized_coordinates


def build_geometry(molecule, qc_params):
    """Create the interface wrapper for a configured software backend."""
    software = normalize_backend_name(qc_params['software'])
    geometry_optimizer = qc_params.get('geometry_optimizer', 'native')
    gamma = qc_params.get('gamma', None)

    if geometry_optimizer == 'geometric':
        validate_backend_capability(
            software,
            ("energy_gradient",),
            context="geomeTRIC AFIR optimisation",
        )
        from pyar.backends import geometric
        return geometric.Geometric(molecule, qc_params)

    if software == 'mlatom_aiqm1':
        from pyar.backends import mlatom_aiqm1
        return mlatom_aiqm1.MlatomAiqm1(molecule, qc_params)
    if software == "gaussian":
        from pyar.backends import gaussian
        return gaussian.Gaussian(molecule, qc_params)
    if software == "orca":
        from pyar.backends import orca
        return orca.Orca(molecule, qc_params)
    if software == "orca-aiqm1":
        from pyar.backends import orca_aiqm1
        geometry = orca_aiqm1.OrcaAIQM1(molecule, qc_params)
        if gamma is not None:
            geometry.set_gamma(gamma)
        return geometry
    if software == "psi4":
        from pyar.backends import psi4
        return psi4.Psi4(molecule, qc_params)
    if software == "xtb":
        from pyar.backends import xtb
        return xtb.Xtb(molecule, qc_params)
    if software == 'xtb_turbo':
        if gamma == 0.0:
            from pyar.backends import xtb
            return xtb.Xtb(molecule, qc_params)
        from pyar.backends import xtb_turbo
        return xtb_turbo.XtbTurbo(molecule, qc_params)
    if software == 'turbomole':
        from pyar.backends import turbomole
        return turbomole.Turbomole(molecule, qc_params)
    if software == "mopac":
        from pyar.backends import mopac
        return mopac.Mopac(molecule, qc_params)
    if software == "aimnet_2":
        from pyar.backends import aimnet_2
        return aimnet_2.Aimnet2(molecule, qc_params)
    if software == "aiqm1_mlatom":
        from pyar.backends import aiqm1_mlatom
        return aiqm1_mlatom.AIQM1(molecule, qc_params)
    if software == "xtb-aimnet2":
        from pyar.backends import xtb_aimnet2
        return xtb_aimnet2.XtbAimnet2(molecule, qc_params)
    if software == "xtb-aiqm1":
        from pyar.backends import xtb_aiqm1
        return xtb_aiqm1.XtbAIQM1(molecule, qc_params)
    if software == 'obabel':
        from pyar.backends import babel
        return babel.OBabel(molecule)

    raise ValueError(f"Unsupported software backend: {software}")


def optimise(molecule, qc_params):
    """Run one optimization job and return the backend status object."""
    original_coordinates = (
        None if molecule.coordinates is None else np.asarray(molecule.coordinates, dtype=float).copy()
    )
    original_energy = molecule.energy
    had_optimized_coordinates = hasattr(molecule, "optimized_coordinates")
    original_optimized_coordinates = getattr(molecule, "optimized_coordinates", None)
    if original_optimized_coordinates is not None:
        original_optimized_coordinates = np.asarray(original_optimized_coordinates, dtype=float).copy()
    cwd = os.getcwd()
    if molecule.name == '':
        molecule.name = 'Opt job'
    job_dir = f'job_{molecule.name}'
    if not os.path.exists(job_dir):
        file_manager.make_directories(job_dir)
    final_coordinate_path = os.path.relpath(
        os.path.join(job_dir, f"result_{molecule.name}.xyz"),
        cwd,
    )
    os.chdir(job_dir)
    software = qc_params.get('software', 'unknown')
    try:
        gamma = qc_params.get('gamma', None)
        optimiser_logger.debug(
            "Optimization job start: name=%s software=%s geometry_optimizer=%s gamma=%s dir=%s",
            molecule.name,
            software,
            qc_params.get('geometry_optimizer', 'native'),
            gamma,
            os.path.abspath(job_dir),
        )
        backend_provenance = _backend_cache_provenance(qc_params)
        try:
            request_digest = (
                _request_digest(molecule, qc_params, backend_provenance)
                if backend_provenance is not None else None
            )
        except (TypeError, ValueError) as exc:
            request_digest = None
            optimiser_logger.info("Optimization cache disabled for %s: %s", molecule.name, exc)
        force_geometric_restart = bool(
            qc_params.get('bias_controller_restart') or qc_params.get('release_retry_attempt')
        )
        if force_geometric_restart:
            request_digest = None
        cached = _verified_cache(
            molecule,
            request_digest if not force_geometric_restart else None,
            backend_provenance,
            geometric=qc_params.get('geometry_optimizer') == 'geometric',
        )
        if cached is not None:
            cached_status, molecule.energy, molecule.coordinates = cached
            molecule.optimized_coordinates = molecule.coordinates.copy()
            optimiser_logger.info(
                "Optimization reused verified cache: name=%s status=%s energy=%15.6f",
                molecule.name,
                cached_status,
                float(molecule.energy),
            )
            optimiser_logger.info(f'     {molecule.name:35s}: {molecule.energy:15.6f}')
            return {
                "converged": True,
                "cycle_exceeded": "CycleExceeded",
                "release_candidate": "ReleaseCandidate",
            }[cached_status]

        geometry = build_geometry(molecule, qc_params)
        optimize_status = geometry.optimize()
        normalized_status = _status_label(optimize_status)
        if is_usable(optimize_status):
            try:
                result_coordinates = np.asarray(geometry.optimized_coordinates, dtype=float).reshape((-1, 3))
                result_energy = float(geometry.energy)
                result_symbols = list(getattr(geometry, "atoms_list", molecule.atoms_list))
                valid_result = (
                    result_symbols == list(molecule.atoms_list)
                    and result_coordinates.shape == (molecule.number_of_atoms, 3)
                    and np.isfinite(result_coordinates).all()
                    and math.isfinite(result_energy)
                )
            except (AttributeError, TypeError, ValueError) as exc:
                raise BackendExecutionError(
                    f"{software} returned an invalid optimization result for {molecule.name}: {exc}"
                ) from exc
            if not valid_result:
                raise BackendExecutionError(
                    f"{software} returned invalid atoms, coordinates, or energy for {molecule.name}"
                )
            apply_geometry_result(molecule, geometry)
            if request_digest is not None:
                try:
                    _write_optimization_cache(
                        molecule,
                        request_digest,
                        optimize_status,
                        backend_provenance,
                        geometric=qc_params.get('geometry_optimizer') == 'geometric',
                    )
                except BackendExecutionError as exc:
                    # A valid calculator result remains valid if cache bookkeeping fails.
                    # It will simply be recomputed on a later invocation without a verified record.
                    optimiser_logger.warning("Could not record optimization cache: %s", exc)
            optimiser_logger.info(
                "Optimization completed! status=%s final_coordinate=%s",
                normalized_status,
                final_coordinate_path,
            )
            optimiser_logger.info(f'     {molecule.name:35s}: {float(molecule.energy):15.6f}')
        else:
            molecule.energy = None
            molecule.coordinates = None
            optimiser_logger.warning(
                "Optimization failed: name=%s status=%s software=%s",
                molecule.name,
                normalized_status,
                software,
            )
            optimizer = qc_params.get('geometry_optimizer')
            output_file = (
                'geometric.out' if optimizer == 'geometric'
                else 'xtb.out' if software == 'xtb' else None
            )
            output_hint = (
                f"; see {os.path.join(os.getcwd(), output_file)}" if output_file else ""
            )
            raise BackendExecutionError(
                f"{software} optimization failed for {molecule.name}: "
                f"{normalized_status}{output_hint}"
            )
        return optimize_status
    except BackendExecutionError:
        molecule.energy = original_energy
        molecule.coordinates = original_coordinates
        if had_optimized_coordinates:
            molecule.optimized_coordinates = original_optimized_coordinates
        elif hasattr(molecule, "optimized_coordinates"):
            del molecule.optimized_coordinates
        raise
    except Exception as exc:
        optimiser_logger.exception("Optimization crashed: name=%s", molecule.name)
        molecule.energy = original_energy
        molecule.coordinates = original_coordinates
        if had_optimized_coordinates:
            molecule.optimized_coordinates = original_optimized_coordinates
        elif hasattr(molecule, "optimized_coordinates"):
            del molecule.optimized_coordinates
        raise BackendExecutionError(
            f"{software} optimization crashed for {molecule.name}: {exc}"
        ) from exc
    finally:
        os.chdir(cwd)


def write_csv_file(csv_filename, energy_dict):
    import csv
    with open(csv_filename, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Name", "Energy"])
        writer.writerows(energy_dict.items())


def bulk_optimize(input_molecules, qc_params):
    """Optimize molecules and keep only usable results."""
    # Modern requests resolve spin on every molecule, including mixed-parity pools.
    # Keep the historical parameter path unchanged for existing callers.
    status_list = []
    backend_failures = []
    for each_mol in input_molecules:
        parameters = ({
            **qc_params, 'charge': each_mol.charge,
            'multiplicity': each_mol.multiplicity, 'scftype': each_mol.scftype,
            'xtb_unpaired_electrons': each_mol.multiplicity - 1,
        } if qc_params.get('_electronic_state_from_molecule') else qc_params)
        try:
            status_list.append(optimise(each_mol, parameters))
        except BackendExecutionError as exc:
            status_list.append(None)
            backend_failures.append(f"{each_mol.name}: {exc}")
            optimiser_logger.error("Backend failed for %s: %s", each_mol.name, exc)
    status_counts = {}
    for status in status_list:
        key = _status_label(status)
        status_counts[key] = status_counts.get(key, 0) + 1
    optimiser_logger.info(
        "Bulk optimization summary: total=%d usable=%d status_counts=%s",
        len(input_molecules),
        sum(1 for s in status_list if is_usable(s)),
        status_counts,
    )
    if backend_failures:
        raise BackendExecutionError(
            "all input jobs were attempted; one or more molecule optimizations failed: "
            + "; ".join(backend_failures)
        )
    return [
        n
        for n, s in zip(input_molecules, status_list)
        if is_usable(s)
    ]


def main():
    pass


if __name__ == '__main__':
    main()
