"""Batch execution and scientific qualification of fixed geometries."""

import csv
import hashlib
from pathlib import Path
import time

import numpy as np

from pyar import __version__
from pyar.singlepoint.io import (
    BOHR_ANGSTROM,
    geometry_hash,
    validate_result,
    write_json,
    write_xyz,
)
from pyar.singlepoint.native import EVALUATORS
from pyar.singlepoint.vibrations import (
    analyze_hessian,
    gradient_metrics,
    numerical_hessian,
    qualification,
)


def make_evaluator(settings, molecule, directory):
    def evaluate(coordinates, label="reference", gradient=True):
        target = directory / "backend" / label
        target.mkdir(parents=True, exist_ok=False)
        energy, values, metadata = EVALUATORS[settings["software"]](
            settings, molecule, coordinates, target, gradient
        )
        energy, values = validate_result(
            energy, values, molecule.number_of_atoms, gradient
        )
        write_json(
            target / "evaluation.json",
            {
                "energy_hartree": energy,
                "gradient_hartree_per_bohr": None
                if values is None
                else values.tolist(),
                "coordinates_bohr": np.asarray(coordinates).tolist(),
                "metadata": metadata,
            },
        )
        return energy, values

    return evaluate


def frequency_result(evaluate, molecule, coordinates, gradient, directory, controls):
    metrics = gradient_metrics(
        gradient, controls["max_gradient"], controls["rms_gradient"]
    )
    hessian, asymmetry = numerical_hessian(evaluate, coordinates, controls["step"])
    np.savetxt(
        directory / "hessian.txt",
        hessian,
        header="Cartesian Hessian; Hartree/bohr^2; input atom order",
    )
    analysis = analyze_hessian(
        hessian,
        molecule.coordinates,
        molecule.atomic_mass,
        controls["imaginary_threshold"],
    )
    write_json(
        directory / "normal_modes.json",
        {
            "normal_modes_cartesian": analysis.pop("normal_modes_cartesian"),
            "masses_amu": analysis["masses_amu"],
        },
    )
    state = qualification(metrics["stationary"], analysis["imaginary_frequency_count"])
    return {
        **metrics,
        **analysis,
        "qualification": state,
        "minimum_confirmed": state == "minimum_confirmed",
        "hessian_source": "central_finite_difference_of_gradients",
        "hessian_step_bohr": controls["step"],
        "hessian_asymmetry_max": asymmetry,
        "hessian_gradient_evaluations": 6 * molecule.number_of_atoms,
        "warning": None
        if state == "minimum_confirmed"
        else "Minimum not confirmed. Do not use as a validated minimum or for minimum-based thermochemistry. "
        "Review the gradient and imaginary modes, then optimize at the intended level.",
    }


def input_record(source, molecule, settings, workflow):
    return {
        "schema_version": 1,
        "workflow": workflow,
        "pyar_version": __version__,
        "source_file": str(Path(source).resolve()),
        "source_sha256": hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        "geometry_sha256": geometry_hash(molecule.atoms_list, molecule.coordinates),
        "symbols": molecule.atoms_list,
        "coordinates_angstrom": molecule.coordinates.tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
        "settings": settings,
        "geometry_level": "unknown; not inferred from XYZ title",
    }


def calculate_one(workflow, source, molecule, settings, root, controls):
    directory = root / f"job_{molecule.name}"
    directory.mkdir()
    record = input_record(source, molecule, settings, workflow)
    write_json(directory / "request.json", record)
    write_xyz(directory / "input.xyz", molecule.atoms_list, molecule.coordinates)
    start = time.perf_counter()
    result = {
        "schema_version": 1,
        "name": molecule.name,
        "workflow": workflow,
        "geometry_sha256": record["geometry_sha256"],
        "units": {
            "energy": "Hartree",
            "gradient": "Hartree/bohr",
            "coordinates": "Angstrom",
        },
        "backend_evidence": "backend/reference/evaluation.json",
        "directory": str(directory),
        "calculation_status": "failed",
        "qualification": "not_assessed",
        "energy_hartree": None,
    }
    try:
        result.update(evaluate_one(workflow, molecule, settings, directory, controls))
        result["calculation_status"] = "completed"
    except Exception as exc:
        result.update(error=str(exc), exception_type=type(exc).__name__)
    result["wall_seconds"] = time.perf_counter() - start
    write_json(directory / "result.json", result)
    return result


def evaluate_one(workflow, molecule, settings, directory, controls):
    coordinates = molecule.coordinates / BOHR_ANGSTROM
    evaluate = make_evaluator(settings, molecule, directory)
    energy, gradient = evaluate(coordinates, gradient=workflow != "energy")
    result = {
        "energy_hartree": energy,
        "gradient_hartree_per_bohr": None,
        "qualification": "minimum_not_assessed",
    }
    if gradient is not None:
        result["gradient_hartree_per_bohr"] = gradient.tolist()
        np.savetxt(
            directory / "gradient.txt",
            gradient,
            header="Energy gradient; Hartree/bohr; input atom order; force = -gradient",
        )
        result.update(gradient_metrics(gradient))
    if workflow == "frequency":
        result.update(
            frequency_result(
                evaluate, molecule, coordinates, gradient, directory, controls
            )
        )
    return result


def count_unqualified(results):
    completed = [item for item in results if item["calculation_status"] == "completed"]
    return sum(item["qualification"] != "minimum_confirmed" for item in completed)


def batch_status(failed, unqualified, total):
    if failed:
        return "failed" if failed == total else "partial_failed"
    return "scientific_gate_failed" if unqualified else "completed"


def summarize(workflow, root, results):
    failed = sum(item["calculation_status"] == "failed" for item in results)
    unqualified = count_unqualified(results) if workflow == "frequency" else 0
    status = batch_status(failed, unqualified, len(results))
    return {
        "schema_version": 1,
        "workflow": workflow,
        "status": status,
        "run_directory": str(root),
        "input_count": len(results),
        "completed_count": len(results) - failed,
        "failed_count": failed,
        "unqualified_count": unqualified,
        "results": results,
    }


def write_summary(root, report):
    write_json(root / "summary.json", report)
    columns = (
        "name",
        "calculation_status",
        "qualification",
        "energy_hartree",
        "directory",
        "error",
    )
    with (root / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report["results"])


def run_batch(workflow, sources, molecules, settings, root, controls):
    root.mkdir(parents=True, exist_ok=False)
    write_json(
        root / "request.json",
        {
            "schema_version": 1,
            "workflow": workflow,
            "pyar_version": __version__,
            "settings": settings,
            "controls": controls,
            "inputs": [
                input_record(s, m, settings, workflow)
                for s, m in zip(sources, molecules)
            ],
        },
    )
    results = []
    for source, molecule in zip(sources, molecules):
        result = calculate_one(workflow, source, molecule, settings, root, controls)
        results.append(result)
        print(
            f"{molecule.name}: {result['calculation_status']}; {result['qualification']}; E={result['energy_hartree']}"
        )
        report = summarize(workflow, root, results)
        write_summary(root, report)
    return report
