"""Unbiased geomeTRIC reaction-path stages with validated artifact handoff."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from pyar.data import defualt_parameters


STAGES = ("relax", "neb", "ts", "frequency", "irc", "endpoints")
_STAGE_OPTIONS = {
    "relax": ("product_relaxation_fmax", "product_relaxation_max_steps"),
    "neb": ("images", "max_cycles", "max_gradient", "average_gradient", "spring", "climb", "align"),
    "ts": ("ts_max_cycles",),
    "frequency": ("imaginary_frequency_threshold",),
    "irc": ("irc_max_cycles",),
    "endpoints": ("endpoint_max_cycles", "imaginary_frequency_threshold", "irc_endpoint_rmsd_tolerance"),
}
_STAGE_GATES = {
    "relax": ("reactant_relaxation_converged", "product_relaxation_converged",
              "reactant_connectivity_survived_relaxation", "product_connectivity_survived_relaxation"),
    "neb": ("converged", "interior_maximum"),
    "ts": ("ts_optimization_converged",),
    "frequency": ("first_order_saddle_confirmed",),
    "irc": ("irc_converged",),
    "endpoints": ("reactant_product_connection_confirmed",),
}


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


def build_neb_images(start, end, ts_guess, image_count):
    """Build an odd-sized band through the supplied TS waypoint."""
    if image_count < 3 or image_count % 2 != 1:
        raise ValueError("image_count must be an odd integer of at least 3")
    frames = [read_xyz(path) for path in (start, end, ts_guess)]
    symbols = frames[0][0]
    for path, (other_symbols, _, _) in zip((end, ts_guess), frames[1:]):
        if other_symbols != symbols:
            raise ValueError(f"Atom count and ordered elements in {path} must match {start}")
    midpoint = image_count // 2
    images = []
    for index in range(image_count):
        if index <= midpoint:
            fraction, left, right = index / midpoint, frames[0][1], frames[2][1]
        else:
            fraction = (index - midpoint) / midpoint
            left, right = frames[2][1], frames[1][1]
        images.append((1 - fraction) * left + fraction * right)
    return symbols, images


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
    return {key: value for key, value in calculator.qc_params.items() if key != "nprocs"}


def _save_stage(output, stage, result, calculator, inputs, artifacts, dependencies=(), parameters=None):
    result = dict(result, stage=stage, status="complete", schema_version=2,
                  qc_params=_physical_settings(calculator), dependencies=list(dependencies))
    result["parameters"] = dict(parameters or {})
    result["legacy_parameters_unverified"] = False
    result["dependency_hashes"] = {
        dependency: _hash(output / f"{dependency}_summary.json") for dependency in dependencies
    }
    result["backend_model"] = (
        "g-xTB (--gxtb)" if str(calculator.qc_params.get("software")).lower() == "xtb"
        else calculator.qc_params.get("method")
    )
    result["inputs"] = {str(Path(path).resolve()): _hash(path) for path in inputs}
    result["artifacts"] = {str(Path(path).resolve()): _hash(path) for path in artifacts}
    _write_json(output / f"{stage}_summary.json", result)
    return result


def _load_stage(output, stage, calculator, visited=None, *,
                expected_parameters=None, expected_by_stage=None,
                reuse_legacy_summaries=False):
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
    if result.get("qc_params") != _physical_settings(calculator):
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
        )
        # A schema-1 dependency may have been migrated during recursive
        # loading. Do not accept a schema-2 parent whose recorded dependency
        # hash now points to the pre-migration file.
        if (schema_version == 2 and result.get("dependency_hashes", {}).get(dependency)
                != _hash(dependency_path)):
            raise ValueError(f"Stale {stage} dependency: {dependency}; rerun {stage}")
        if not all(upstream.get(key) is True for key in _STAGE_GATES[dependency]):
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
    elif (expected_parameters is not None
          and not result.get("legacy_parameters_unverified")
          and result.get("parameters") != expected_parameters):
        raise ValueError(f"Stage {stage} used different stage-specific parameters; rerun it")
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
    frames = [np.asarray(frame).copy() for frame in progress.xyzs]
    energies = [float(energy) for energy in progress.qm_energies]
    _write_xyz_trajectory(output / f"{label}_path.xyz", symbols, frames, energies)
    return frames, energies, optimizer.state == OPT_STATE.CONVERGED


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

    def load(stage_name):
        return _load_stage(
            output, stage_name, calculator,
            expected_parameters=expected_by_stage[stage_name],
            expected_by_stage=expected_by_stage,
            reuse_legacy_summaries=options["reuse_legacy_summaries"],
        )

    def save(result, inputs=(), artifacts=(), dependencies=()):
        return _save_stage(output, stage, result, calculator, inputs,
                           [output / name for name in artifacts], dependencies,
                           {key: options[key] for key in _STAGE_OPTIONS[stage]})

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
        symbols, frames = build_neb_images(start_path, end_path, guess, options["images"])
        molecule = Molecule()
        molecule.elem, molecule.xyzs = symbols, frames
        engine = EngineASE(molecule, calculator)
        params = NEBParams(images=options["images"], neb_maxcyc=options["max_cycles"],
                           maxg=options["max_gradient"], avgg=options["average_gradient"],
                           nebk=options["spring"], climb=options["climb"], align=options["align"])
        scratch = tempfile.mkdtemp(prefix="geometric_neb_", dir=output)
        chain, cycles = OptimizeChain(ElasticBand(molecule, engine, scratch, params, plain=0), engine, params)
        frames = [structure.M.xyzs[0].copy() for structure in chain.Structures]
        energies = [float(structure.energy) for structure in chain.Structures]
        highest = int(np.argmax(energies))
        _write_xyz_trajectory(output / "neb_path.xyz", symbols, frames, energies)
        interior = 0 < highest < len(frames) - 1
        artifacts = ["neb_path.xyz"]
        if interior:
            _write_xyz_trajectory(output / "ts_guess.xyz", symbols, [frames[highest]], [energies[highest]])
            artifacts.append("ts_guess.xyz")
        return save({
            "converged": bool(chain.avgg <= options["average_gradient"] and chain.maxg <= options["max_gradient"]),
            "images": len(frames), "cycles": int(cycles),
            "interior_maximum": interior, "highest_energy_image_index": highest,
            "highest_energy_hartree": energies[highest],
            "average_rms_gradient_ev_per_angstrom": float(chain.avgg),
            "maximum_rms_gradient_ev_per_angstrom": float(chain.maxg),
        }, [start_path, end_path, guess], artifacts, ["relax"])

    if stage == "ts":
        guess = options["ts_geometry"] or options["ts_guess"]
        dependencies = []
        if guess is None:
            neb = load("neb")
            if not neb["converged"] or not neb["interior_maximum"]:
                raise ValueError("TS optimization requires a converged NEB with an interior maximum")
            guess, dependencies = output / "ts_guess.xyz", ["neb"]
        symbols, geometry, _ = read_xyz(guess)
        frames, energies, converged = _optimize_geometry(
            symbols, geometry, calculator, output, "ts", options["ts_max_cycles"], transition=True,
        )
        _write_xyz_trajectory(output / "ts_optimized.xyz", symbols, [frames[-1]], [energies[-1]])
        return save({"ts_optimization_converged": converged, "ts_energy_hartree": energies[-1]},
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
                          nprocs=1, max_gradient=0.05, average_gradient=0.025, spring=1.0,
                          climb=0.5, align=False, product_relaxation_fmax=0.05,
                          product_relaxation_max_steps=200, ts_max_cycles=200,
                          irc_max_cycles=200, imaginary_frequency_threshold=20.0,
                          irc_endpoint_rmsd_tolerance=0.5, stage="all", ts_geometry=None,
                          endpoint_max_cycles=300, reuse_legacy_summaries=False):
    from pyar.backends.geometric import PyarGeometricCalculator

    options = dict(locals())
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
    qc_params = dict(software=software, method=method or defualt_parameters.values["method"],
                     basis=basis or defualt_parameters.values["basis"], charge=charge,
                     multiplicity=multiplicity, nprocs=nprocs, gamma=0.0)
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
    parser.add_argument("--output", default="neb_run")
    parser.add_argument("--images", type=int, default=11)
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
