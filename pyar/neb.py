"""Unbiased geomeTRIC reaction-path stages with validated artifact handoff."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from numbers import Integral
import os
from pathlib import Path
import tempfile
import time
from importlib.metadata import PackageNotFoundError, version as package_version

import numpy as np

from pyar.data import defualt_parameters
from pyar.backends.xtb_utils import canonical_xtb_model


STAGES = ("relax", "neb", "ts", "frequency", "irc", "endpoints")
_STAGE_OPTIONS = {
    "relax": ("product_relaxation_fmax", "product_relaxation_max_steps"),
    "neb": ("images", "max_cycles", "max_gradient", "average_gradient", "spring", "climb", "align",
            "interpolation", "idpp_fmax", "idpp_steps", "geodesic_tol", "geodesic_max_iter"),
    "ts": ("ts_max_cycles", "ts_optimizer", "ts_fmax", "sella_fmax"),
    "frequency": ("imaginary_frequency_threshold",),
    "irc": ("irc_max_cycles",),
    "endpoints": ("endpoint_max_cycles", "imaginary_frequency_threshold", "irc_endpoint_rmsd_tolerance"),
}
_LEGACY_NEB_INITIALIZATION_DEFAULTS = {
    "interpolation": "linear",
    "idpp_fmax": 0.1,
    "idpp_steps": 100,
}
_NEB_COMMON_PARAMETERS = (
    "images", "max_cycles", "max_gradient", "average_gradient", "spring", "climb", "align",
)
_IDPP_PARAMETERS = ("idpp_fmax", "idpp_steps")
_GEODESIC_PARAMETERS = ("geodesic_tol", "geodesic_max_iter")
_TS_SHARED_PARAMETERS = ("ts_max_cycles",)
_TS_SELLA_PARAMETERS = ("sella_fmax",)
_NEB_PARAMETER_KEYS = frozenset(
    (*_NEB_COMMON_PARAMETERS, "interpolation", *_IDPP_PARAMETERS, *_GEODESIC_PARAMETERS)
)
_STAGE_GATES = {
    "relax": ("reactant_relaxation_converged", "product_relaxation_converged",
              "reactant_connectivity_survived_relaxation", "product_connectivity_survived_relaxation"),
    "neb": ("converged",),
    "ts": ("ts_optimization_converged",),
    "frequency": ("first_order_saddle_confirmed",),
    "irc": ("irc_converged",),
    "endpoints": ("reactant_product_connection_confirmed",),
}


def _geometric_climbing_state(chain):
    """Return explicit geomeTRIC climbing state without guessing from gradients."""
    activated = getattr(chain, "climbSet", None)
    if not isinstance(activated, (bool, np.bool_)):
        return False, None, []
    if not activated:
        return True, False, []
    climbers = getattr(chain, "climbers", None)
    if climbers is None:
        return True, True, []
    try:
        return True, True, list(climbers)
    except TypeError:
        return True, True, []


def _select_neb_ts_candidate(frames, energies, *, climbing_activated, climber_indices):
    """Choose a valid geomeTRIC climber, or the historical interior maximum."""
    if not frames or not energies:
        raise ValueError("NEB must return at least one image and energy")
    highest = int(np.argmax(energies))
    interior = 0 < highest < len(frames) - 1

    valid_climbers = []
    for value in climber_indices:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
            continue
        index = int(value)
        if not 0 < index < len(frames) - 1 or index >= len(energies):
            continue
        if not np.all(np.isfinite(np.asarray(frames[index], dtype=float))):
            continue
        if not np.isfinite(energies[index]):
            continue
        if index not in valid_climbers:
            valid_climbers.append(index)

    selected_climber = None
    if climbing_activated and valid_climbers:
        selected_climber = max(
            valid_climbers, key=lambda index: (energies[index], -index),
        )
        selected, source = selected_climber, "climbing_image"
    elif interior and highest < len(frames) and np.isfinite(energies[highest]) \
            and np.all(np.isfinite(np.asarray(frames[highest], dtype=float))):
        selected, source = highest, "highest_energy_interior_image"
    else:
        selected, source = None, None

    return {
        "highest_energy_image_index": highest,
        "interior_maximum": interior,
        "climbing_image_indices": valid_climbers,
        "selected_climbing_image_index": selected_climber,
        "ts_guess_image_index": selected,
        "ts_guess_source": source,
    }


def _stage_gate_passed(result, stage):
    """Apply gates while retaining the interior-image rule for old NEB summaries."""
    if stage == "neb":
        if result.get("converged") is not True:
            return False
        if "ts_guess_image_index" in result:
            return result["ts_guess_image_index"] is not None
        return result.get("interior_maximum") is True
    return all(result.get(key) is True for key in _STAGE_GATES[stage])


def read_xyz(path):
    """Read exactly one finite, nonempty XYZ frame in atom-mapped order."""
    from ase.data import atomic_numbers

    lines = Path(path).read_text().splitlines()
    try:
        natoms = int(lines[0].strip())
        if natoms < 1 or len(lines) < natoms + 2:
            raise ValueError("missing atom records")
        symbols, coordinates = [], []
        for line in lines[2:natoms + 2]:
            fields = line.split()
            symbol = fields[0].capitalize()
            if symbol not in atomic_numbers or symbol == "X" or len(fields) < 4:
                raise ValueError("invalid atom record")
            symbols.append(symbol)
            coordinates.append([float(value) for value in fields[1:4]])
        coordinates = np.asarray(coordinates)
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("non-finite coordinates")
        if any(line.strip() for line in lines[natoms + 2:]):
            raise ValueError("expected exactly one XYZ frame")
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid XYZ file {path}: {exc}") from exc
    return symbols, coordinates, lines[1]


def _linear_neb_images(start_coordinates, ts_coordinates, end_coordinates, image_count):
    """Build the existing piecewise-linear path through the TS waypoint."""
    midpoint = image_count // 2
    images = []
    for index in range(image_count):
        if index <= midpoint:
            fraction, left, right = index / midpoint, start_coordinates, ts_coordinates
        else:
            fraction = (index - midpoint) / midpoint
            left, right = ts_coordinates, end_coordinates
        images.append((1 - fraction) * left + fraction * right)
    return images


def _ase_idpp_api():
    """Load ASE's IDPP implementation across old and current module paths."""
    try:
        from ase import Atoms
    except ImportError as exc:
        raise RuntimeError("ASE is required for IDPP initialization") from exc

    try:
        idpp_module = importlib.import_module("ase.mep.neb")
    except ImportError:
        # ASE moved NEB helpers from ase.neb to ase.mep.neb. Keep IDPP usable
        # with installations that still expose the legacy module path.
        try:
            idpp_module = importlib.import_module("ase.neb")
        except ImportError as exc:
            raise RuntimeError(
                "Installed ASE does not provide IDPP interpolation; "
                "tried ase.mep.neb and ase.neb"
            ) from exc
    idpp_interpolate = getattr(idpp_module, "idpp_interpolate", None)
    if idpp_interpolate is None:
        raise RuntimeError(
            "Installed ASE does not provide idpp_interpolate in "
            "ase.mep.neb or ase.neb"
        )
    return Atoms, idpp_interpolate


def _idpp_refine_segment(symbols, coordinates, *, label, fmax, steps):
    """Refine one piecewise-linear segment with ASE IDPP, without file output."""
    if len(coordinates) <= 2:
        return [np.asarray(frame, dtype=float).copy() for frame in coordinates]

    Atoms, idpp_interpolate = _ase_idpp_api()

    images = [Atoms(symbols=symbols, positions=frame) for frame in coordinates]
    fixed_start = images[0].get_positions().copy()
    fixed_end = images[-1].get_positions().copy()
    try:
        idpp_interpolate(images, traj=None, log=None, fmax=fmax, steps=steps)
    except Exception as exc:
        raise RuntimeError(f"ASE IDPP initialization failed for {label}: {exc}") from exc

    refined = [image.get_positions().copy() for image in images]
    # ASE NEB holds endpoints fixed; restore their exact input values as an
    # explicit contract, including across ASE-version implementation details.
    refined[0] = fixed_start
    refined[-1] = fixed_end
    if any(not np.all(np.isfinite(frame)) for frame in refined):
        raise ValueError(f"ASE IDPP initialization produced non-finite coordinates for {label}")
    return refined


def _installed_geodesic_version():
    """Return installed geodesic distribution provenance without importing it."""
    try:
        return package_version("geodesic-interpolate")
    except PackageNotFoundError:
        return None


def _geodesic_api():
    """Load the optional geodesic package and its installed distribution version."""
    try:
        from geodesic_interpolate import Geodesic
    except ImportError as exc:
        raise RuntimeError(
            "Geodesic NEB initialization requires the optional "
            "'geodesic-interpolate' dependency. Install PyAR with "
            "`pip install 'pyar-chem[geodesic]'`."
        ) from exc
    installed_version = _installed_geodesic_version()
    if installed_version is None:
        raise RuntimeError(
            "The geodesic_interpolate module is importable, but its "
            "geodesic-interpolate distribution version is unavailable."
        )
    return Geodesic, installed_version


def _align_rigid_frame(coordinates, reference):
    """Rigidly align coordinates to reference using a proper Kabsch rotation."""
    moving = np.asarray(coordinates, dtype=float)
    target = np.asarray(reference, dtype=float)
    moving_center = moving.mean(axis=0)
    target_center = target.mean(axis=0)
    moving_centered = moving - moving_center
    target_centered = target - target_center
    left, _, right = np.linalg.svd(moving_centered.T @ target_centered)
    correction = np.eye(3)
    correction[-1, -1] = 1.0 if np.linalg.det(left @ right) >= 0 else -1.0
    rotation = left @ correction @ right
    return moving_centered @ rotation + target_center


def _geodesic_refine_segment(
    symbols,
    seed_coordinates,
    *,
    tol=0.002,
    max_iter=15,
    geodesic_class=None,
):
    """Smooth one fixed-endpoint path segment and restore its input frames."""
    seed = np.asarray(seed_coordinates, dtype=float).copy()
    if seed.ndim != 3 or seed.shape[1:] != (len(symbols), 3):
        raise ValueError("Geodesic seed coordinates have an unexpected shape")
    if len(seed) <= 2:
        return [frame.copy() for frame in seed]
    if geodesic_class is None:
        geodesic_class, _ = _geodesic_api()

    try:
        smoother = geodesic_class(
            list(symbols), seed.copy(), scaler=1.7, threshold=3.0,
            min_neighbors=4, friction=0.01,
        )
        smoothed = np.asarray(smoother.smooth(tol=tol, max_iter=max_iter), dtype=float)
    except Exception as exc:
        raise RuntimeError(f"Geodesic smoothing failed: {exc}") from exc
    if smoothed.shape != seed.shape:
        raise RuntimeError(
            "Geodesic smoothing returned an unexpected path shape "
            f"{smoothed.shape}; expected {seed.shape}"
        )
    if not np.all(np.isfinite(smoothed)):
        raise ValueError("Geodesic smoothing produced non-finite coordinates")

    restored = [
        _align_rigid_frame(frame, reference)
        for frame, reference in zip(smoothed, seed)
    ]
    restored[0] = seed[0].copy()
    restored[-1] = seed[-1].copy()
    if any(not np.all(np.isfinite(frame)) for frame in restored):
        raise ValueError("Geodesic frame restoration produced non-finite coordinates")
    return restored


def canonical_neb_parameters(parameters):
    """Return only NEB settings that affect the selected initializer."""
    values = dict(parameters)
    unknown = set(values) - _NEB_PARAMETER_KEYS
    if unknown:
        names = ", ".join(sorted(unknown))
        raise ValueError(f"Unrecognized NEB stage parameter(s): {names}")
    interpolation = values.get("interpolation", "linear")
    canonical = {key: values[key] for key in _NEB_COMMON_PARAMETERS if key in values}
    canonical["interpolation"] = interpolation
    if interpolation in {"idpp", "geodesic"}:
        canonical.update({
            "idpp_fmax": values.get("idpp_fmax", _LEGACY_NEB_INITIALIZATION_DEFAULTS["idpp_fmax"]),
            "idpp_steps": values.get("idpp_steps", _LEGACY_NEB_INITIALIZATION_DEFAULTS["idpp_steps"]),
        })
    if interpolation == "geodesic":
        canonical.update({
            "geodesic_tol": values.get("geodesic_tol", 0.002),
            "geodesic_max_iter": values.get("geodesic_max_iter", 15),
        })
    return canonical


def canonical_ts_parameters(parameters):
    """Return only TS settings active for the selected optimizer."""
    values = dict(parameters)
    unknown = set(values) - {*_TS_SHARED_PARAMETERS, "ts_optimizer", "ts_fmax", *_TS_SELLA_PARAMETERS}
    if unknown:
        names = ", ".join(sorted(unknown))
        raise ValueError(f"Unrecognized TS stage parameter(s): {names}")
    optimizer = values.get("ts_optimizer", "geometric")
    if optimizer not in {"geometric", "sella"}:
        raise ValueError("ts_optimizer must be 'geometric' or 'sella'")
    canonical = {key: values[key] for key in _TS_SHARED_PARAMETERS if key in values}
    canonical["ts_optimizer"] = optimizer
    ts_fmax = values.get("ts_fmax")
    if ts_fmax is not None:
        try:
            if isinstance(ts_fmax, (bool, np.bool_)):
                raise ValueError
            ts_fmax = float(ts_fmax)
        except (TypeError, ValueError, OverflowError):
            raise ValueError("ts_fmax must be positive and finite (eV/angstrom)") from None
        if not np.isfinite(ts_fmax) or ts_fmax <= 0:
            raise ValueError("ts_fmax must be positive and finite (eV/angstrom)")
    canonical["ts_fmax"] = ts_fmax
    if optimizer == "sella" and ts_fmax is None:
        sella_fmax = values.get("sella_fmax", 0.05)
        try:
            if isinstance(sella_fmax, (bool, np.bool_)):
                raise ValueError
            sella_fmax = float(sella_fmax)
        except (TypeError, ValueError, OverflowError):
            raise ValueError("sella_fmax must be positive and finite (eV/angstrom)")
        if not np.isfinite(sella_fmax) or sella_fmax <= 0:
            raise ValueError("sella_fmax must be positive and finite (eV/angstrom)")
        canonical["sella_fmax"] = sella_fmax
    return canonical


def build_neb_images(
    start, end, ts_guess, image_count, *, interpolation="linear",
    idpp_fmax=0.1, idpp_steps=100, geodesic_tol=0.002, geodesic_max_iter=15,
):
    """Build an odd-sized band through the supplied TS waypoint.

    ``linear`` retains the historical piecewise Cartesian interpolation.
    ``idpp`` refines the two halves independently while fixing both endpoints
    and the supplied TS waypoint. ``geodesic`` adds geodesic smoothing to the
    same deterministic IDPP seed on each side of that waypoint.
    """
    if image_count < 3 or image_count % 2 != 1:
        raise ValueError("image_count must be an odd integer of at least 3")
    interpolation = str(interpolation).lower()
    if interpolation not in {"linear", "idpp", "geodesic"}:
        raise ValueError(
            f"Unknown interpolation {interpolation!r}; expected 'linear', 'idpp', or 'geodesic'"
        )
    if interpolation in {"idpp", "geodesic"}:
        _validate_idpp_options(idpp_fmax, idpp_steps)
    if interpolation == "geodesic":
        _validate_geodesic_options(geodesic_tol, geodesic_max_iter)
        geodesic_class, _ = _geodesic_api()
    else:
        geodesic_class = None
    frames = [read_xyz(path) for path in (start, end, ts_guess)]
    symbols = frames[0][0]
    for path, (other_symbols, _, _) in zip((end, ts_guess), frames[1:]):
        if other_symbols != symbols:
            raise ValueError(f"Atom count and ordered elements in {path} must match {start}")
    start_coordinates, end_coordinates, ts_coordinates = (
        frames[0][1], frames[1][1], frames[2][1]
    )
    images = _linear_neb_images(
        start_coordinates, ts_coordinates, end_coordinates, image_count,
    )
    if interpolation in {"idpp", "geodesic"}:
        midpoint = image_count // 2
        first_half = _idpp_refine_segment(
            symbols, images[:midpoint + 1], label="reactant-to-TS segment",
            fmax=idpp_fmax, steps=idpp_steps,
        )
        second_half = _idpp_refine_segment(
            symbols, images[midpoint:], label="TS-to-product segment",
            fmax=idpp_fmax, steps=idpp_steps,
        )
        images = first_half + second_half[1:]
        images[0] = start_coordinates.copy()
        images[midpoint] = ts_coordinates.copy()
        images[-1] = end_coordinates.copy()
        if len(images) != image_count:
            raise RuntimeError("IDPP initialization returned an unexpected number of images")
        if any(not np.all(np.isfinite(frame)) for frame in images):
            raise ValueError("IDPP initialization produced non-finite coordinates")
    if interpolation == "geodesic":
        midpoint = image_count // 2
        first_half = _geodesic_refine_segment(
            symbols, images[:midpoint + 1], tol=geodesic_tol,
            max_iter=geodesic_max_iter, geodesic_class=geodesic_class,
        )
        second_half = _geodesic_refine_segment(
            symbols, images[midpoint:], tol=geodesic_tol,
            max_iter=geodesic_max_iter, geodesic_class=geodesic_class,
        )
        images = first_half + second_half[1:]
        images[0] = start_coordinates.copy()
        images[midpoint] = ts_coordinates.copy()
        images[-1] = end_coordinates.copy()
        if len(images) != image_count:
            raise RuntimeError("Geodesic initialization returned an unexpected number of images")
        if any(not np.all(np.isfinite(frame)) for frame in images):
            raise ValueError("Geodesic initialization produced non-finite coordinates")
    return symbols, images


def _validate_idpp_options(fmax, steps):
    if not np.isfinite(fmax) or fmax <= 0:
        raise ValueError("idpp_fmax must be positive and finite")
    if not isinstance(steps, int) or steps < 1:
        raise ValueError("idpp_steps must be a positive integer")


def _validate_geodesic_options(tol, max_iter):
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("geodesic_tol must be positive and finite")
    if not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("geodesic_max_iter must be a positive integer")


def _bond_set(symbols, coordinates, scale=1.3):
    from pyar.data.new_atomic_data import covalent_radius

    return {
        (left, right)
        for left in range(len(symbols)) for right in range(left + 1, len(symbols))
        if np.linalg.norm(coordinates[left] - coordinates[right]) < scale * (
            covalent_radius[symbols[left]] + covalent_radius[symbols[right]]
        )
    }


def _aligned_rmsd(first, second):
    """Atom-mapped RMSD after rigid alignment without reflections."""
    first, second = np.asarray(first), np.asarray(second)
    if first.shape != second.shape:
        return float("inf")
    first, second = first - first.mean(axis=0), second - second.mean(axis=0)
    left, _, right = np.linalg.svd(first.T @ second)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right)
    return float(np.sqrt(np.mean(np.sum((first @ left @ correction @ right - second)**2, axis=1))))


def _write_xyz_trajectory(path, symbols, frames, energies):
    if len(frames) != len(energies) or not len(frames):
        raise ValueError("Geometry and energy frame counts must agree and be nonempty")
    with Path(path).open("w") as stream:
        for index, (coordinates, energy) in enumerate(zip(frames, energies)):
            coordinates = np.asarray(coordinates)
            if (coordinates.shape != (len(symbols), 3)
                    or not np.all(np.isfinite(coordinates)) or not np.isfinite(energy)):
                raise ValueError("Cannot write non-finite or malformed trajectory data")
            stream.write(f"{len(symbols)}\nimage={index} energy_hartree={energy:.14f}\n")
            for symbol, xyz in zip(symbols, coordinates):
                stream.write(f"{symbol:<2s} {xyz[0]: .14f} {xyz[1]: .14f} {xyz[2]: .14f}\n")


def _write_json(path, data):
    temporary = Path(path).with_suffix(".json.tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _physical_settings(calculator):
    """Exclude execution resources, but bind artifacts to the physical method."""
    qc_params = calculator.qc_params
    settings = {key: value for key, value in qc_params.items()
                if key not in {"nprocs", "xtb_model"}}
    if str(qc_params.get("software", "")).lower() == "xtb":
        settings["xtb_model"] = canonical_xtb_model(qc_params.get("xtb_model"))
    return settings


def _save_stage(output, stage, result, calculator, inputs, artifacts, dependencies=(), parameters=None):
    result = dict(result, stage=stage, status="complete", schema_version=2,
                  qc_params=_physical_settings(calculator), dependencies=list(dependencies))
    result["parameters"] = dict(parameters or {})
    result["legacy_parameters_unverified"] = False
    result["dependency_hashes"] = {
        dependency: _hash(output / f"{dependency}_summary.json") for dependency in dependencies
    }
    if str(calculator.qc_params.get("software", "")).lower() == "xtb":
        model = canonical_xtb_model(calculator.qc_params.get("xtb_model"))
        result["backend_model"] = (
            "g-xTB (--gxtb)" if model == "gxtb" else "GFN2-xTB (--gfn 2)"
        )
    else:
        result["backend_model"] = calculator.qc_params.get("method")
    result["inputs"] = {str(Path(path).resolve()): _hash(path) for path in inputs}
    result["artifacts"] = {str(Path(path).resolve()): _hash(path) for path in artifacts}
    _write_json(output / f"{stage}_summary.json", result)
    return result


def _load_stage(output, stage, calculator, visited=None, *,
                expected_parameters=None, expected_by_stage=None,
                reuse_legacy_summaries=False, expected_geodesic_version=None):
    """Validate a stage, optionally migrating verified schema-1 summaries."""
    visited = set() if visited is None else visited
    if stage in visited:
        raise ValueError("Cyclic stage dependency")
    visited = visited | {stage}
    path = output / f"{stage}_summary.json"
    if not path.is_file():
        raise ValueError(f"Run --stage {stage} first; missing {path.name}")
    result = json.loads(path.read_text())
    schema_version = result.get("schema_version")
    if result.get("status") != "complete" or schema_version not in {1, 2}:
        raise ValueError(f"Stage {stage} is incomplete; rerun it")
    if schema_version == 1 and not reuse_legacy_summaries:
        raise ValueError(
            f"Stage {stage} uses a legacy summary without stage-parameter provenance. "
            "Pass --reuse-legacy-summaries to verify its artifacts and reuse it, "
            "or rerun the stage."
        )
    recorded_qc_params = result.get("qc_params", {})
    requested_qc_params = _physical_settings(calculator)
    # Prior schema-2 xTB stages always ran the hard-coded --gxtb command.
    # Preserve their reuse under the historical default while distinguishing
    # them from explicitly requested GFN2-xTB stages.
    if (str(recorded_qc_params.get("software", "")).lower() == "xtb"
            and "xtb_model" not in recorded_qc_params
            and result.get("backend_model") == "g-xTB (--gxtb)"):
        recorded_qc_params = dict(recorded_qc_params, xtb_model="gxtb")
    if recorded_qc_params != requested_qc_params:
        raise ValueError(f"Stage {stage} used different backend or electronic-structure settings")
    for group in ("inputs", "artifacts"):
        for filename, digest in result[group].items():
            if not Path(filename).is_file() or _hash(filename) != digest:
                raise ValueError(f"Stale {stage} artifact: {filename}; rerun that stage")
    for dependency in result["dependencies"]:
        dependency_path = output / f"{dependency}_summary.json"
        if not dependency_path.is_file():
            raise ValueError(f"Stale {stage} dependency: {dependency}; rerun {stage}")
        if (schema_version == 2 and result.get("dependency_hashes", {}).get(dependency)
                != _hash(dependency_path)):
            raise ValueError(f"Stale {stage} dependency: {dependency}; rerun {stage}")
        upstream = _load_stage(
            output, dependency, calculator, visited,
            expected_parameters=(expected_by_stage or {}).get(dependency),
            expected_by_stage=expected_by_stage,
            reuse_legacy_summaries=(reuse_legacy_summaries or schema_version == 1),
            expected_geodesic_version=expected_geodesic_version,
        )
        # A schema-1 dependency may have been migrated during recursive
        # loading. Do not accept a schema-2 parent whose recorded dependency
        # hash now points to the pre-migration file.
        if (schema_version == 2 and result.get("dependency_hashes", {}).get(dependency)
                != _hash(dependency_path)):
            raise ValueError(f"Stale {stage} dependency: {dependency}; rerun {stage}")
        if not _stage_gate_passed(upstream, dependency):
            raise ValueError(f"Stage {stage} depends on scientifically invalid {dependency}; rerun that stage")
    if schema_version == 1:
        # Legacy records have strong hashes for their numerical files and
        # methods, but no stage-option record. Preserve that uncertainty rather
        # than claiming the current CLI options produced these artifacts.
        result.update(
            schema_version=2,
            parameters=None,
            legacy_parameters_unverified=True,
            legacy_reuse_authorized=True,
            dependency_hashes={
                dependency: _hash(output / f"{dependency}_summary.json")
                for dependency in result["dependencies"]
            },
        )
        _write_json(path, result)
    elif expected_parameters is not None and not result.get("legacy_parameters_unverified"):
        recorded_parameters = result.get("parameters")
        if stage == "neb" and isinstance(recorded_parameters, dict):
            recorded_parameters = canonical_neb_parameters(recorded_parameters)
        elif stage == "ts" and isinstance(recorded_parameters, dict):
            recorded_parameters = canonical_ts_parameters(recorded_parameters)
        parameters_match = recorded_parameters == expected_parameters
        if not parameters_match:
            raise ValueError(f"Stage {stage} used different stage-specific parameters; rerun it")
        if (stage == "neb" and expected_parameters.get("interpolation") == "geodesic"):
            recorded_version = result.get("geodesic_interpolate_version")
            if (recorded_version is not None and expected_geodesic_version is not None
                    and recorded_version != expected_geodesic_version):
                raise ValueError(
                    "Stage neb used geodesic-interpolate version "
                    f"{recorded_version}; installed version is {expected_geodesic_version}; rerun it"
                )
    return result


def _engine(symbols, coordinates, calculator):
    from geometric.ase_engine import EngineASE
    from geometric.molecule import Molecule

    molecule = Molecule()
    molecule.elem = list(symbols)
    molecule.xyzs = [np.asarray(coordinates).copy()]
    molecule.build_topology()
    return molecule, EngineASE(molecule, calculator)


def _optimize_geometry(symbols, coordinates, calculator, output, label, max_cycles,
                       *, transition=False, direction=None, hessian=None, fmax=None):
    """Run one geomeTRIC optimization or IRC direction and retain failures."""
    from ase.units import Bohr, Hartree
    from geometric.errors import GeomOptNotConvergedError
    from geometric.internal import DelocalizedInternalCoordinates
    from geometric.optimize import OPT_STATE, Optimizer
    from geometric.params import OptParams

    molecule, engine = _engine(symbols, coordinates, calculator)
    internal = DelocalizedInternalCoordinates(molecule, build=True, connect=False, addcart=False)
    # Fresh scratch prevents Hessian reuse across distinct geometries/methods.
    scratch = tempfile.mkdtemp(prefix=f"geometric_{label}_", dir=output)
    options = dict(maxiter=max_cycles, xyzout=str(output / f"{label}_path.xyz"),
                   convergence_set="GAU_TIGHT", frequency=False)
    if fmax is not None:
        # geomeTRIC tolerances are in Hartree/Bohr; the public option is eV/angstrom.
        options["converge"] = ["gmax", str(fmax * Bohr / Hartree),
                               "grms", str(fmax * Bohr / Hartree / 1.5)]
    if transition:
        options.update(transition=True, hessian="first")
    if direction:
        # geomeTRIC's IRC step uses the mass-weighted modes computed here.
        options.update(irc=True, irc_direction=direction, hessian=f"file:{hessian}", frequency=True)
    optimizer = Optimizer(
        np.asarray(coordinates).flatten() / Bohr, molecule, internal, engine,
        scratch, OptParams(**options), print_info=False,
    )
    try:
        progress = optimizer.optimizeGeometry()
    except GeomOptNotConvergedError:
        progress = optimizer.progress
    steps = getattr(optimizer, "Iteration", None)
    calculator._pyar_optimizer_steps = None if steps is None else int(steps)
    frames = [np.asarray(frame).copy() for frame in progress.xyzs]
    energies = [float(energy) for energy in progress.qm_energies]
    _write_xyz_trajectory(output / f"{label}_path.xyz", symbols, frames, energies)
    return frames, energies, optimizer.state == OPT_STATE.CONVERGED


def _sella_api():
    """Load the optional Sella saddle optimizer and return its distribution version."""
    try:
        from sella import Sella
    except ImportError as exc:
        raise RuntimeError(
            "Sella TS optimization requires the optional 'sella' dependency. "
            "Install PyAR with `pip install 'pyar-chem[sella]'`."
        ) from exc
    try:
        installed_version = package_version("sella")
    except PackageNotFoundError as exc:
        raise RuntimeError(
            "The Sella module is importable, but its installed distribution version is unavailable."
        ) from exc
    return Sella, installed_version


def _optimize_sella(symbols, coordinates, calculator, max_steps, fmax):
    """Optimize a TS candidate with Sella using PyAR's unbiased ASE calculator."""
    from ase import Atoms
    from ase.units import Hartree

    Sella, sella_version = _sella_api()
    atoms = Atoms(symbols=list(symbols), positions=np.asarray(coordinates, dtype=float))
    atoms.calc = calculator
    frames, energies = [], []

    def record_frame():
        frames.append(np.asarray(atoms.get_positions(), dtype=float).copy())
        energies.append(float(atoms.get_potential_energy()) / Hartree)

    try:
        optimizer = Sella(atoms, logfile=None, order=1, internal=False)
        optimizer.attach(record_frame, interval=1)
        converged = bool(optimizer.run(fmax=fmax, steps=max_steps))
        steps = getattr(optimizer, "nsteps", None)
        calculator._pyar_optimizer_steps = None if steps is None else int(steps)
    except Exception as exc:
        raise RuntimeError(f"Sella TS optimization failed: {exc}") from exc

    if not frames:
        raise RuntimeError("Sella returned an empty TS optimization trajectory")
    try:
        final_coordinates = np.asarray(atoms.get_positions(), dtype=float).copy()
        final_energy = float(atoms.get_potential_energy()) / Hartree
    except Exception as exc:
        raise RuntimeError(f"Sella final-state evaluation failed: {exc}") from exc
    if (not np.array_equal(frames[-1], final_coordinates)
            or energies[-1] != final_energy):
        frames.append(final_coordinates)
        energies.append(final_energy)
    if any(frame.shape != (len(symbols), 3) or not np.all(np.isfinite(frame)) for frame in frames):
        raise ValueError("Sella TS optimization produced malformed or non-finite coordinates")
    if any(not np.isfinite(energy) for energy in energies):
        raise ValueError("Sella TS optimization produced a non-finite energy")
    return frames, energies, converged, sella_version


def _relax_endpoint(symbols, coordinates, calculator, output, label, fmax, max_steps):
    """Require an input reaction endpoint to survive unbiased optimization."""
    frames, energies, converged = _optimize_geometry(
        symbols, coordinates, calculator, output, label, max_steps, fmax=fmax,
    )
    if not converged:
        raise RuntimeError(f"{label.capitalize()} endpoint did not converge during unbiased relaxation")
    bonds = _bond_set(symbols, frames[-1])
    if bonds != _bond_set(symbols, coordinates):
        raise RuntimeError(f"{label.capitalize()} endpoint did not survive unbiased relaxation: connectivity changed")
    return frames[-1], energies[-1], bonds


def _frequency(symbols, coordinates, calculator, output, label, threshold):
    """Verify stationarity and Hessian index on the unbiased physical surface."""
    from ase.units import Bohr
    from geometric.normal_modes import calc_cartesian_hessian, frequency_analysis

    # A numerically bent linear minimum can lose a bending mode when geomeTRIC
    # projects six rigid modes instead of five. Resolve only sub-1e-5 A departures
    # from a line, then evaluate BOTH gradients and Hessian at that geometry.
    coordinates = np.asarray(coordinates, dtype=float).copy()
    center = coordinates.mean(axis=0)
    centered = coordinates - center
    _, _, axes = np.linalg.svd(centered, full_matrices=False)
    line = np.outer(centered @ axes[0], axes[0]) + center
    linear_correction = float(np.max(np.linalg.norm(line - coordinates, axis=1)))
    linearized = len(symbols) > 2 and 0 < linear_correction <= 1e-5
    if linearized:
        coordinates = line
    molecule, engine = _engine(symbols, coordinates, calculator)
    scratch = tempfile.mkdtemp(prefix=f"geometric_{label}_frequency_", dir=output)
    coords_bohr = np.asarray(coordinates).flatten() / Bohr
    evaluation = engine.calc(coords_bohr, scratch)
    gradient = np.asarray(evaluation["gradient"]).reshape(-1, 3)
    norms = np.linalg.norm(gradient, axis=1)
    max_gradient = float(np.max(norms))
    rms_gradient = float(np.sqrt(np.mean(norms**2)))
    hessian = calc_cartesian_hessian(coords_bohr.copy(), molecule, engine, scratch, read_data=False)
    if hessian.shape != (coords_bohr.size, coords_bohr.size) or not np.all(np.isfinite(hessian)):
        raise ValueError("Invalid Cartesian Hessian")
    # Central differences have small numerical asymmetry; use the symmetric Hessian.
    hessian = (hessian + hessian.T) / 2
    hessian_path = output / f"{label}_hessian.txt"
    np.savetxt(hessian_path, hessian)
    wavenumbers, _, _ = frequency_analysis(
        coords_bohr, hessian, elem=list(symbols), energy=float(evaluation["energy"]),
        outfnm=str(output / f"{label}_frequencies.vdata"),
    )
    wavenumbers = np.asarray(wavenumbers)
    if not np.all(np.isfinite(wavenumbers)):
        raise ValueError("Non-finite vibrational frequencies")
    imaginary = wavenumbers[wavenumbers < -threshold]
    stationary = max_gradient < 4.5e-4 and rms_gradient < 3.0e-4
    return {
        "evaluated_coordinates_angstrom": coordinates.tolist(),
        "linear_geometry_correction_angstrom": linear_correction if linearized else 0.0,
        "energy_hartree": float(evaluation["energy"]),
        "frequencies_cm-1": wavenumbers.tolist(),
        "imaginary_frequencies_cm-1": imaginary.tolist(),
        "imaginary_frequency_count": int(imaginary.size),
        "imaginary_frequency_threshold_cm-1": threshold,
        "maximum_gradient_hartree_per_bohr": max_gradient,
        "rms_gradient_hartree_per_bohr": rms_gradient,
        "stationary": stationary,
        "first_order_saddle_confirmed": stationary and imaginary.size == 1,
        "minimum_confirmed": stationary and imaginary.size == 0,
    }


def _match_endpoints(symbols, observed, expected, tolerance):
    """Require a bijection: two branches returning to one minimum must fail."""
    rmsds = [[_aligned_rmsd(obs, exp) for exp in expected] for obs in observed]
    bonds = [[_bond_set(symbols, obs) == _bond_set(symbols, exp)
              for exp in expected] for obs in observed]
    assignments = ((0, 1), (1, 0))
    observed_rmsd = _aligned_rmsd(observed[0], observed[1])
    distinct = (_bond_set(symbols, observed[0]) != _bond_set(symbols, observed[1])
                or observed_rmsd > tolerance)
    connectivity = any(all(bonds[i][j] for i, j in enumerate(order)) for order in assignments)
    matching = [order for order in assignments if distinct and all(
        bonds[i][j] and rmsds[i][j] <= tolerance for i, j in enumerate(order)
    )]
    return {
        "observed_endpoints_distinct": distinct,
        "observed_endpoint_rmsd_angstrom": observed_rmsd,
        "irc_endpoint_connectivities_match": connectivity,
        "irc_endpoint_geometries_match": bool(matching),
        "irc_endpoint_rmsd_matrix_angstrom": rmsds,
        "endpoint_assignment": list(matching[0]) if matching else None,
    }


def _execute_stage(stage, output, calculator, options):
    """Execute a single stage; dependent results are checked before use."""
    expected_by_stage = {
        name: {key: options[key] for key in keys}
        for name, keys in _STAGE_OPTIONS.items()
    }
    expected_by_stage["neb"] = canonical_neb_parameters(expected_by_stage["neb"])
    expected_by_stage["ts"] = canonical_ts_parameters(expected_by_stage["ts"])
    # Optimizer choices describe how a TS was produced. Once the TS summary
    # and its hashes are validated, downstream scientific stages consume the
    # artifact independently of the optimizer currently selected or installed.
    if options["stage"] in {"frequency", "irc", "endpoints"}:
        expected_by_stage["ts"] = None

    def load(stage_name):
        return _load_stage(
            output, stage_name, calculator,
            expected_parameters=expected_by_stage[stage_name],
            expected_by_stage=expected_by_stage,
            reuse_legacy_summaries=options["reuse_legacy_summaries"],
            expected_geodesic_version=options.get("geodesic_interpolate_version"),
        )

    def save(result, inputs=(), artifacts=(), dependencies=()):
        if stage == "neb" and options["interpolation"] == "geodesic":
            result["geodesic_interpolate_version"] = options["geodesic_interpolate_version"]
        parameters = (
            canonical_neb_parameters({key: options[key] for key in _STAGE_OPTIONS["neb"]})
            if stage == "neb"
            else canonical_ts_parameters({key: options[key] for key in _STAGE_OPTIONS["ts"]})
            if stage == "ts"
            else {key: options[key] for key in _STAGE_OPTIONS[stage]}
        )
        return _save_stage(output, stage, result, calculator, inputs,
                           [output / name for name in artifacts], dependencies, parameters)

    start_path, end_path = output / "reactant_relaxed.xyz", output / "product_relaxed.xyz"
    if stage == "relax":
        start, end = options["start"], options["end"]
        symbols, start_xyz, _ = read_xyz(start)
        end_symbols, end_xyz, _ = read_xyz(end)
        if symbols != end_symbols:
            raise ValueError("Start and end atom count and ordered elements must match")
        result = {}
        for label, geometry, path in (("reactant", start_xyz, start_path), ("product", end_xyz, end_path)):
            relaxed, energy, bonds = _relax_endpoint(
                symbols, geometry, calculator, output, label,
                options["product_relaxation_fmax"], options["product_relaxation_max_steps"],
            )
            _write_xyz_trajectory(path, symbols, [relaxed], [energy])
            result.update({f"{label}_relaxation_converged": True,
                           f"{label}_connectivity_survived_relaxation": True,
                           f"{label}_energy_hartree": energy,
                           f"{label}_bonds_zero_based": sorted(map(list, bonds))})
        return save(result, [start, end], [start_path.name, end_path.name])

    if stage == "neb":
        from geometric.ase_engine import EngineASE
        from geometric.molecule import Molecule
        from geometric.neb import ElasticBand, OptimizeChain
        from geometric.params import NEBParams

        load("relax")
        guess = options["ts_guess"]
        symbols, frames = build_neb_images(
            start_path, end_path, guess, options["images"],
            interpolation=options["interpolation"],
            idpp_fmax=options["idpp_fmax"],
            idpp_steps=options["idpp_steps"],
            geodesic_tol=options["geodesic_tol"],
            geodesic_max_iter=options["geodesic_max_iter"],
        )
        molecule = Molecule()
        molecule.elem, molecule.xyzs = symbols, frames
        engine = EngineASE(molecule, calculator)
        params = NEBParams(images=options["images"], neb_maxcyc=options["max_cycles"],
                           maxg=options["max_gradient"], avgg=options["average_gradient"],
                           nebk=options["spring"], climb=options["climb"], align=options["align"],
                           ncimg=1)
        scratch = tempfile.mkdtemp(prefix="geometric_neb_", dir=output)
        chain, cycles = OptimizeChain(ElasticBand(molecule, engine, scratch, params, plain=0), engine, params)
        frames = [structure.M.xyzs[0].copy() for structure in chain.Structures]
        energies = [float(structure.energy) for structure in chain.Structures]
        climbing_state_available, climbing_activated, climber_indices = _geometric_climbing_state(chain)
        candidate = _select_neb_ts_candidate(
            frames, energies, climbing_activated=climbing_activated is True,
            climber_indices=climber_indices,
        )
        highest = candidate["highest_energy_image_index"]
        _write_xyz_trajectory(output / "neb_path.xyz", symbols, frames, energies)
        artifacts = ["neb_path.xyz"]
        selected = candidate["ts_guess_image_index"]
        if selected is not None:
            _write_xyz_trajectory(output / "ts_guess.xyz", symbols, [frames[selected]], [energies[selected]])
            artifacts.append("ts_guess.xyz")
        return save({
            "converged": bool(chain.avgg <= options["average_gradient"] and chain.maxg <= options["max_gradient"]),
            "images": len(frames), "cycles": int(cycles),
            "climbing_image_state_available": climbing_state_available,
            "climbing_image_activated": climbing_activated,
            "climbing_image_indices": candidate["climbing_image_indices"],
            "selected_climbing_image_index": candidate["selected_climbing_image_index"],
            "requested_climbing_images": 1,
            "interior_maximum": candidate["interior_maximum"],
            "highest_energy_image_index": highest,
            "ts_guess_image_index": selected,
            "ts_guess_source": candidate["ts_guess_source"],
            "highest_energy_hartree": energies[highest],
            "average_rms_gradient_ev_per_angstrom": float(chain.avgg),
            "maximum_rms_gradient_ev_per_angstrom": float(chain.maxg),
        }, [start_path, end_path, guess], artifacts, ["relax"])

    if stage == "ts":
        guess = options["ts_geometry"] or options["ts_guess"]
        dependencies = []
        if guess is None:
            neb = load("neb")
            has_candidate = (
                neb["ts_guess_image_index"] is not None
                if "ts_guess_image_index" in neb else neb["interior_maximum"]
            )
            if not neb["converged"] or not has_candidate:
                raise ValueError(
                    "TS optimization requires a converged NEB with an interior maximum "
                    "or a usable climbing image"
                )
            guess, dependencies = output / "ts_guess.xyz", ["neb"]
        symbols, geometry, _ = read_xyz(guess)
        evaluation_count_start = getattr(calculator, "backend_energy_gradient_evaluations", 0)
        optimization_started = time.perf_counter()
        effective_fmax = (
            options["ts_fmax"] if options["ts_fmax"] is not None
            else options["sella_fmax"] if options["ts_optimizer"] == "sella"
            else None
        )
        if options["ts_optimizer"] == "sella":
            frames, energies, converged, sella_version = _optimize_sella(
                symbols, geometry, calculator, options["ts_max_cycles"], effective_fmax,
            )
            _write_xyz_trajectory(output / "ts_path.xyz", symbols, frames, energies)
        else:
            frames, energies, converged = _optimize_geometry(
                symbols, geometry, calculator, output, "ts", options["ts_max_cycles"],
                transition=True, fmax=effective_fmax,
            )
        optimization_wall_seconds = time.perf_counter() - optimization_started
        optimizer_steps = getattr(calculator, "_pyar_optimizer_steps", None)
        _write_xyz_trajectory(output / "ts_optimized.xyz", symbols, [frames[-1]], [energies[-1]])
        result = {
            "ts_optimization_converged": converged,
            "ts_energy_hartree": energies[-1],
            "backend_energy_gradient_evaluations": (
                getattr(calculator, "backend_energy_gradient_evaluations", 0) - evaluation_count_start
            ),
            "optimizer_steps": optimizer_steps,
            "wall_seconds": float(optimization_wall_seconds),
            "effective_ts_convergence": {
                "max_steps": options["ts_max_cycles"],
                "fmax_ev_per_angstrom": effective_fmax,
                "geometric_convergence_set": (
                    "GAU_TIGHT" if options["ts_optimizer"] == "geometric"
                    and effective_fmax is None else None
                ),
            },
        }
        if options["ts_optimizer"] == "sella":
            result["sella_version"] = sella_version
        return save(result,
                    [guess], ["ts_optimized.xyz", "ts_path.xyz"], dependencies)

    if stage == "frequency":
        geometry_path = options["ts_geometry"]
        dependencies = []
        if geometry_path is None:
            ts = load("ts")
            if not ts["ts_optimization_converged"]:
                raise ValueError("TS optimization did not converge; cannot validate a saddle")
            geometry_path, dependencies = output / "ts_optimized.xyz", ["ts"]
        symbols, geometry, _ = read_xyz(geometry_path)
        result = _frequency(symbols, geometry, calculator, output, "ts", options["imaginary_frequency_threshold"])
        _write_xyz_trajectory(output / "frequency_geometry.xyz", symbols,
                              [result["evaluated_coordinates_angstrom"]], [result["energy_hartree"]])
        return save(result, [geometry_path], ["frequency_geometry.xyz", "ts_hessian.txt", "ts_frequencies.vdata"], dependencies)

    if stage == "irc":
        frequency = load("frequency")
        if not frequency["first_order_saddle_confirmed"]:
            raise ValueError("IRC requires a stationary TS with exactly one significant imaginary frequency")
        symbols, geometry, _ = read_xyz(output / "frequency_geometry.xyz")
        result, paths, artifacts = {"irc_run": True}, [], []
        for direction in ("forward", "backward"):
            frames, energies, converged = _optimize_geometry(
                symbols, geometry, calculator, output, f"irc_{direction}", options["irc_max_cycles"],
                direction=direction, hessian=output / "ts_hessian.txt",
            )
            result[f"irc_{direction}_converged"] = converged
            path = f"irc_{direction}_endpoint.xyz"
            _write_xyz_trajectory(output / path, symbols, [frames[-1]], [energies[-1]])
            artifacts.extend([path, f"irc_{direction}_path.xyz"])
            paths.append((frames, energies))
        _write_xyz_trajectory(output / "irc_path.xyz", symbols,
                              paths[0][0][::-1] + paths[1][0][1:],
                              paths[0][1][::-1] + paths[1][1][1:])
        result["irc_converged"] = result["irc_forward_converged"] and result["irc_backward_converged"]
        result["reactant_product_connection_confirmed"] = False
        return save(result, [output / "frequency_geometry.xyz", output / "ts_hessian.txt"],
                    artifacts + ["irc_path.xyz"], ["frequency"])

    if stage == "endpoints":
        irc = load("irc")
        load("relax")
        symbols, reactant, _ = read_xyz(start_path)
        end_symbols, product, _ = read_xyz(end_path)
        if symbols != end_symbols:
            raise ValueError("Relaxed endpoint atom order differs")
        result, observed, artifacts, inputs = {}, [], [], [start_path, end_path]
        for direction in ("forward", "backward"):
            path = output / f"irc_{direction}_endpoint.xyz"
            branch_symbols, geometry, _ = read_xyz(path)
            if symbols != branch_symbols:
                raise ValueError("IRC atom order differs from relaxed endpoints")
            label = f"irc_{direction}_relaxed"
            frames, energies, converged = _optimize_geometry(
                symbols, geometry, calculator, output, label, options["endpoint_max_cycles"],
            )
            # Topology is free to change here: classify the actual final minimum.
            observed.append(frames[-1])
            _write_xyz_trajectory(output / f"{label}.xyz", symbols, [frames[-1]], [energies[-1]])
            result[f"{direction}_optimization_converged"] = converged
            result[f"{direction}_topology_changed_on_relaxation"] = _bond_set(symbols, geometry) != _bond_set(symbols, frames[-1])
            verification = _frequency(symbols, frames[-1], calculator, output, label,
                                      options["imaginary_frequency_threshold"])
            observed[-1] = np.asarray(verification["evaluated_coordinates_angstrom"])
            _write_xyz_trajectory(output / f"{label}.xyz", symbols, [observed[-1]], [verification["energy_hartree"]])
            result[f"{direction}_frequency"] = verification
            artifacts.extend([f"{label}.xyz", f"{label}_path.xyz", f"{label}_hessian.txt", f"{label}_frequencies.vdata"])
            inputs.append(path)
        result.update(_match_endpoints(symbols, observed, [reactant, product], options["irc_endpoint_rmsd_tolerance"]))
        result["endpoints_are_minima"] = all(
            result[f"{direction}_optimization_converged"] and result[f"{direction}_frequency"]["minimum_confirmed"]
            for direction in ("forward", "backward")
        )
        distinct = (_bond_set(symbols, reactant) != _bond_set(symbols, product)
                    or _aligned_rmsd(reactant, product) > options["irc_endpoint_rmsd_tolerance"])
        result["reference_endpoints_distinct"] = distinct
        result["reactant_product_connection_confirmed"] = bool(
            irc["irc_converged"] and result["endpoints_are_minima"] and distinct
            and result["observed_endpoints_distinct"]
            and result["irc_endpoint_connectivities_match"] and result["irc_endpoint_geometries_match"]
        )
        return save(result, inputs, artifacts, ["irc", "relax"])
    raise ValueError(f"Unknown stage: {stage}")


def _run_neb_in_directory(start, end, ts_guess, *, software, output="neb_run", images=11,
                          max_cycles=100, method=None, basis=None, charge=0, multiplicity=1,
                          nprocs=1, xtb_model="gxtb", max_gradient=0.05, average_gradient=0.025, spring=1.0,
                          climb=0.5, align=False, product_relaxation_fmax=0.05,
                          product_relaxation_max_steps=200, ts_max_cycles=200,
                          irc_max_cycles=200, imaginary_frequency_threshold=20.0,
                          irc_endpoint_rmsd_tolerance=0.5, stage="all", ts_geometry=None,
                          endpoint_max_cycles=300, reuse_legacy_summaries=False,
                          interpolation="linear", idpp_fmax=0.1, idpp_steps=100,
                          geodesic_tol=0.002, geodesic_max_iter=15,
                          ts_optimizer="geometric", ts_fmax=None, sella_fmax=0.05):
    from pyar.backends.geometric import PyarGeometricCalculator

    options = dict(locals())
    options["geodesic_interpolate_version"] = None
    output = Path(output).resolve()
    if stage not in ("all",) + STAGES:
        raise ValueError(f"Unknown stage: {stage}")
    if stage in {"all", "relax"} and (start is None or end is None):
        raise ValueError("This stage requires --start and --end")
    if stage in {"all", "neb"}:
        if ts_guess is None:
            raise ValueError("The NEB stage requires --ts-guess")
        if images < 3 or images % 2 != 1:
            raise ValueError("images must be an odd integer of at least 3")
    if ts_geometry is not None and stage not in {"ts", "frequency"}:
        raise ValueError("--ts-geometry is supported only for ts and frequency stages")
    overwritten = {
        "relax": {"reactant_relaxed.xyz", "product_relaxed.xyz"},
        "neb": {"neb_path.xyz", "ts_guess.xyz"},
        "ts": {"ts_optimized.xyz", "ts_path.xyz"},
        "frequency": {"frequency_geometry.xyz"},
        "irc": set(), "endpoints": set(),
    }
    reserved = set().union(*overwritten.values()) if stage == "all" else overwritten[stage]
    for path in (start, end, ts_guess, ts_geometry):
        if path is not None and Path(path).resolve() in {output / name for name in reserved}:
            raise ValueError("An input would be overwritten by this stage; use a new output directory")
    for key in ("max_cycles", "product_relaxation_max_steps", "ts_max_cycles", "irc_max_cycles", "endpoint_max_cycles", "nprocs", "multiplicity"):
        if not isinstance(options[key], int) or options[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    for key in ("max_gradient", "average_gradient", "spring", "climb", "product_relaxation_fmax", "imaginary_frequency_threshold", "irc_endpoint_rmsd_tolerance"):
        if not np.isfinite(options[key]) or options[key] <= 0:
            raise ValueError(f"{key} must be positive and finite")
    if options["interpolation"] not in {"linear", "idpp", "geodesic"}:
        raise ValueError("interpolation must be 'linear', 'idpp', or 'geodesic'")
    if interpolation in {"idpp", "geodesic"}:
        _validate_idpp_options(idpp_fmax, idpp_steps)
    if interpolation == "geodesic":
        _validate_geodesic_options(geodesic_tol, geodesic_max_iter)
        if stage in {"all", "neb"}:
            _, options["geodesic_interpolate_version"] = _geodesic_api()
        else:
            # Downstream stages consume the validated NEB artifacts and their
            # recorded provenance; the local package version is irrelevant.
            options["geodesic_interpolate_version"] = None
    if stage in {"all", "ts"}:
        canonical_ts_parameters({
            "ts_optimizer": ts_optimizer,
            "ts_max_cycles": ts_max_cycles,
            "ts_fmax": ts_fmax,
            "sella_fmax": sella_fmax,
        })
    qc_params = dict(software=software, method=method or defualt_parameters.values["method"],
                     basis=basis or defualt_parameters.values["basis"], charge=charge,
                     multiplicity=multiplicity, nprocs=nprocs, gamma=0.0)
    if str(software).lower() == "xtb":
        qc_params["xtb_model"] = canonical_xtb_model(xtb_model)
    calculator = PyarGeometricCalculator(qc_params=qc_params)
    stages = STAGES if stage == "all" else (stage,)
    results = {"reactant_product_connection_confirmed": False}
    for current in stages:
        _write_json(output / f"{current}_summary.json", {"stage": current, "status": "running"})
        try:
            # During all, the TS must consume the optimized NEB maximum, not the input waypoint.
            current_options = dict(options)
            if stage == "all" and current == "ts":
                current_options["ts_guess"] = None
            result = _execute_stage(current, output, calculator, current_options)
        except Exception as exc:
            _write_json(output / f"{current}_summary.json", {"stage": current, "status": "failed", "error": str(exc)})
            if stage == "all":
                _write_json(output / "workflow_summary.json", dict(results, status="failed", failed_stage=current, error=str(exc)))
            raise
        results[current] = result
        if stage != "all":
            return result
        results["reactant_product_connection_confirmed"] = bool(result.get("reactant_product_connection_confirmed", False))
        _write_json(output / "workflow_summary.json", dict(results, status="running", completed_stage=current))
    results["status"] = "complete"
    _write_json(output / "workflow_summary.json", results)
    return results


def run_neb(start=None, end=None, ts_guess=None, **kwargs):
    """Isolate backend artifacts in output and restore cwd even on failure."""
    output = Path(kwargs.pop("output", "neb_run")).resolve()
    paths = [Path(path).resolve() if path is not None else None for path in (start, end, ts_guess)]
    if kwargs.get("ts_geometry") is not None:
        kwargs["ts_geometry"] = Path(kwargs["ts_geometry"]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    previous = Path.cwd()
    try:
        os.chdir(output)
        return _run_neb_in_directory(*paths, output=output, **kwargs)
    finally:
        os.chdir(previous)


def _build_parser():
    parser = argparse.ArgumentParser(prog="pyar-neb", description="Run unbiased geomeTRIC reaction-path stages.")
    parser.add_argument("--stage", choices=("all",) + STAGES, default="all")
    parser.add_argument(
        "--reuse-legacy-summaries", action="store_true",
        help="verify and reuse schema-1 stage files; their original stage options were not recorded",
    )
    parser.add_argument("--start", help="Reactant XYZ for all/relax")
    parser.add_argument("--end", help="Product XYZ for all/relax")
    parser.add_argument("--ts-guess", help="NEB waypoint for all/neb, or explicit guess for ts")
    parser.add_argument("--ts-geometry", help="Explicit geometry for ts or frequency")
    parser.add_argument("--software", required=True, help="PyAR energy-gradient backend")
    parser.add_argument(
        "--xtb-model", choices=("gxtb", "gfn2"), default="gxtb",
        help="xTB Hamiltonian for the xTB energy-gradient provider (default: gxtb)",
    )
    parser.add_argument("--output", default="neb_run")
    parser.add_argument("--images", type=int, default=11)
    parser.add_argument(
        "--interpolation", choices=("linear", "idpp", "geodesic"), default="linear",
        help="initial NEB path interpolation (linear is the backward-compatible default)",
    )
    parser.add_argument("--idpp-fmax", type=float, default=0.1,
                        help="ASE IDPP force convergence for idpp initialization")
    parser.add_argument("--idpp-steps", type=int, default=100,
                        help="maximum ASE IDPP optimization steps per path half")
    parser.add_argument("--geodesic-tol", type=float, default=0.002,
                        help="geodesic smoothing tolerance for geodesic initialization")
    parser.add_argument("--geodesic-max-iter", type=int, default=15,
                        help="maximum geodesic smoothing iterations per path half")
    parser.add_argument("--max-cycles", type=int, default=100)
    parser.add_argument("--method", default=defualt_parameters.values["method"])
    parser.add_argument("--basis", default=defualt_parameters.values["basis"])
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--nprocs", type=int, default=1)
    parser.add_argument("--max-gradient", type=float, default=0.05)
    parser.add_argument("--average-gradient", type=float, default=0.025)
    parser.add_argument("--spring", type=float, default=1.0)
    parser.add_argument("--climb", type=float, default=0.5)
    parser.add_argument("--align", action="store_true")
    parser.add_argument("--product-relaxation-fmax", type=float, default=0.05)
    parser.add_argument("--product-relaxation-max-steps", type=int, default=200)
    parser.add_argument("--ts-max-cycles", type=int, default=200)
    parser.add_argument("--ts-optimizer", choices=("geometric", "sella"), default="geometric",
                        help="TS optimizer (geometric is the default; sella is optional)")
    parser.add_argument("--ts-fmax", type=float, default=None,
                        help="optional shared TS force convergence in eV/angstrom for geomeTRIC and Sella")
    parser.add_argument("--sella-fmax", type=float, default=0.05,
                        help="Sella force convergence in eV/angstrom (used only with --ts-optimizer sella)")
    parser.add_argument("--irc-max-cycles", type=int, default=200)
    parser.add_argument("--endpoint-max-cycles", type=int, default=300)
    parser.add_argument("--imaginary-frequency-threshold", type=float, default=20.0)
    parser.add_argument("--irc-endpoint-rmsd-tolerance", type=float, default=0.5)
    return parser


def main(argv=None):
    parser = _build_parser()
    try:
        result = run_neb(**vars(parser.parse_args(argv)))
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"pyar-neb: {exc}\n")
    print(json.dumps(result, indent=2))
    # A completed numerical calculation can still fail scientific validation.
    failed = any(result.get(key) is False for key in (
        "converged", "interior_maximum", "ts_optimization_converged", "first_order_saddle_confirmed",
        "irc_converged", "endpoints_are_minima",
    ))
    if result.get("stage") == "endpoints" or "endpoints" in result:
        failed |= not result.get("reactant_product_connection_confirmed", False)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
