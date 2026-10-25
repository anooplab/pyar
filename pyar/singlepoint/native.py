"""Native fixed-geometry inputs; no inherited optimization keywords or bias."""

import os
import re
import subprocess

import numpy as np

from pyar.backends import require_executable
from pyar.energy_gradient_providers import (
    _read_gaussian_energy,
    _read_gaussian_forces,
    _read_orca_energy,
    _read_orca_engrad,
    _read_xtb_energy,
    _read_xtb_gradient,
)
from pyar.singlepoint.io import BOHR_ANGSTROM, write_json, write_xyz


def run_program(command, directory, settings, *, stdin=None, label="backend"):
    """Retain the exact invocation and both streams, including failed jobs."""
    environment = dict(
        os.environ,
        OMP_NUM_THREADS=str(settings["nprocs"]),
        MKL_NUM_THREADS=str(settings["nprocs"]),
        OPENBLAS_NUM_THREADS=str(settings["nprocs"]),
    )
    write_json(directory / f"{label}_command.json", {"argv": command})
    with (
        (directory / f"{label}.out").open("w") as stdout,
        (directory / f"{label}.err").open("w") as stderr,
    ):
        process = subprocess.run(
            command,
            cwd=directory,
            env=environment,
            input=stdin,
            stdout=stdout,
            stderr=stderr,
            text=True,
            check=False,
        )
    write_json(directory / f"{label}_exit.json", {"returncode": process.returncode})
    output = (
        (directory / f"{label}.out").read_text()
        + "\n"
        + (directory / f"{label}.err").read_text()
    )
    if process.returncode:
        raise RuntimeError(
            f"Backend exited with status {process.returncode}; see {directory}. "
            + output[-1200:]
        )
    return output


def program_metadata(output, backend):
    patterns = {
        "xtb": r"xtb version\s+(\S+)",
        "orca": r"Program Version\s+(\S+)",
        "gaussian": r"(Gaussian\s+\d+[^\n]*)",
    }
    match = re.search(patterns[backend], output, re.I)
    return {
        "backend": backend,
        "version": match.group(1) if match else "unknown; see raw output",
    }


def require_completion(output, marker):
    if marker.lower() not in output.lower():
        raise RuntimeError(f"Backend completion not verified: missing {marker!r}")
    bad = re.search(
        r"(SCF NOT CONVERGED|SCC (?:failed|not converged)|convergence failure|unknown option|error termination)",
        output,
        re.I,
    )
    if bad:
        raise RuntimeError(f"Backend failed: {bad.group(0)}")


def xtb(settings, molecule, coordinates, directory, gradient):
    write_xyz(
        directory / "input.xyz",
        molecule.atoms_list,
        (np.asarray(coordinates) * BOHR_ANGSTROM),
    )
    model = settings["xtb_model"]
    selector = {
        "gfn0": ["--gfn", "0"],
        "gfn1": ["--gfn", "1"],
        "gfn2": ["--gfn", "2"],
        "gfnff": ["--gfnff"],
        "gxtb": ["--gxtb"],
    }[model]
    command = [
        require_executable("xtb"),
        "input.xyz",
        *selector,
        "--chrg",
        str(molecule.charge),
        "--uhf",
        str(molecule.multiplicity - 1),
        "--parallel",
        str(settings["nprocs"]),
        "--acc",
        str(settings["accuracy"]),
        "--iterations",
        str(settings["scf_cycles"]),
    ]
    if gradient:
        command.append("--grad")
    if settings["solvent"]:
        command.extend(["--" + settings["solvent_model"], settings["solvent"]])
    output = run_program(command, directory, settings)
    require_completion(output, "normal termination of xtb")
    expected = {
        "gfn0": "GFN0-xTB",
        "gfn1": "GFN1-xTB",
        "gfn2": "GFN2-xTB",
        "gfnff": "GFN-FF",
        "gxtb": "g-xTB",
    }[model]
    if expected.lower() not in output.lower():
        raise RuntimeError(
            f"Requested xTB Hamiltonian {expected} was not verified in output"
        )
    value = _read_xtb_gradient(directory / "gradient") if gradient else None
    return (
        _read_xtb_energy(output.splitlines()),
        value,
        {**program_metadata(output, "xtb"), "model_verified": expected},
    )


def orca_input(settings, molecule, coordinates, gradient):
    keywords = [settings["method"]]
    if settings["basis"]:
        keywords.append(settings["basis"])
    keywords.extend(["TightSCF", "NoAutoStart"])
    if gradient:
        keywords.append("EnGrad")
    if settings["dispersion"] != "none":
        keywords.append(settings["dispersion"])
    if settings["solvent"]:
        keywords.append(f"{settings['solvent_model']}({settings['solvent']})")
    head = "! " + " ".join(keywords) + "\n"
    head += f"%pal nprocs {settings['nprocs']} end\n%maxcore {settings['memory_mb']}\n"
    head += f"%scf MaxIter {settings['scf_cycles']} end\n"
    head += f"* xyz {molecule.charge} {molecule.multiplicity}\n"
    body = "\n".join(
        f"{s} {x:.15f} {y:.15f} {z:.15f}"
        for s, (x, y, z) in zip(
            molecule.atoms_list, (np.asarray(coordinates) * BOHR_ANGSTROM)
        )
    )
    return head + body + "\n*\n"


def orca(settings, molecule, coordinates, directory, gradient):
    (directory / "input.inp").write_text(
        orca_input(settings, molecule, coordinates, gradient)
    )
    output = run_program([require_executable("orca"), "input.inp"], directory, settings)
    require_completion(output, "ORCA TERMINATED NORMALLY")
    energy = _read_orca_energy(output)
    if not gradient:
        return energy, None, program_metadata(output, "orca")
    gradient_energy, value = _read_orca_engrad(
        directory / "input.engrad", len(molecule.atoms_list)
    )
    if abs(energy - gradient_energy) > 1e-7:
        raise ValueError("ORCA output and engrad energies disagree")
    return energy, value, program_metadata(output, "orca")


def gaussian_input(settings, molecule, coordinates, gradient):
    basis = re.sub(r"^def2-", "def2", settings["basis"], flags=re.I)
    route = f"{settings['method']}/{basis} NoSymm SCF=(Tight,MaxCycle={settings['scf_cycles']})"
    if gradient:
        route += " Force"
    if settings["dispersion"] != "none":
        route += (
            " EmpiricalDispersion="
            + {"D3BJ": "GD3BJ", "D3ZERO": "GD3"}[settings["dispersion"]]
        )
    if settings["solvent"]:
        route += f" SCRF=({settings['solvent_model']},Solvent={settings['solvent']})"
    head = f"%nprocshared={settings['nprocs']}\n%mem={settings['memory_mb']}MB\n#p {route}\n\nPyAR single point\n\n"
    head += f"{molecule.charge} {molecule.multiplicity}\n"
    body = "\n".join(
        f"{s} {x:.15f} {y:.15f} {z:.15f}"
        for s, (x, y, z) in zip(
            molecule.atoms_list, (np.asarray(coordinates) * BOHR_ANGSTROM)
        )
    )
    return head + body + "\n\n"


def gaussian(settings, molecule, coordinates, directory, gradient):
    content = gaussian_input(settings, molecule, coordinates, gradient)
    (directory / "input.com").write_text(content)
    output = run_program(
        [require_executable("g16")], directory, settings, stdin=content
    )
    require_completion(output, "Normal termination of Gaussian")
    lines = output.splitlines()
    value = _read_gaussian_forces(lines, len(molecule.atoms_list)) if gradient else None
    return _read_gaussian_energy(lines), value, program_metadata(output, "gaussian")


def python_backend(settings, molecule, coordinates, directory, gradient):
    """Isolate optional runtimes and process-global calculator state."""
    import sys

    payload = {
        "settings": settings,
        "symbols": molecule.atoms_list,
        "coordinates_angstrom": np.asarray(
            (np.asarray(coordinates) * BOHR_ANGSTROM)
        ).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
        "gradient": gradient,
    }
    write_json(directory / "worker_request.json", payload)
    run_program([sys.executable, "-m", "pyar.singlepoint.worker"], directory, settings)
    import json

    result = json.loads((directory / "worker_result.json").read_text())
    return (
        result["energy_hartree"],
        result.get("gradient_hartree_per_bohr"),
        result["metadata"],
    )


EVALUATORS = {
    "xtb": xtb,
    "orca": orca,
    "gaussian": gaussian,
    "aimnet_2": python_backend,
    "ani": python_backend,
    "turbomole": python_backend,
}
