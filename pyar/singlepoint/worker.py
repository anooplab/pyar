"""Isolated optional Python calculators. Invoked inside each evaluation folder."""

import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np

from pyar.singlepoint.io import write_json


def aimnet_calculator(request):
    from pyar.AIMNet2.calculators.aimnet2ase import AIMNet2Calculator
    from pyar.backends.aimnet_2 import load_default_aimnet2_model

    return AIMNet2Calculator(load_default_aimnet2_model(), charge=request["charge"])


def ani_calculator(request):
    import torchani

    name = request["settings"]["model"]
    model = getattr(torchani.models, name)().double()
    return model.ase()


def turbomole_calculator(request):
    from ase.calculators.turbomole import Turbomole

    settings = request["settings"]
    functional = {"bp86": "b-p", "b3lyp": "b3-lyp"}.get(
        settings["method"].lower(), settings["method"].lower()
    )
    return Turbomole(
        multiplicity=request["multiplicity"],
        **{
            "total charge": request["charge"],
            "basis set name": settings["basis"],
            "use dft": functional != "hf",
            "density functional": functional,
            "use resolution of identity": False,
            "scf iterations": settings["scf_cycles"],
            "scf energy convergence": 1e-8 * 27.211386245988,
        },
    )


CALCULATORS = {
    "aimnet_2": aimnet_calculator,
    "ani": ani_calculator,
    "turbomole": turbomole_calculator,
}


def calculator_metadata(request):
    backend = request["settings"]["software"]
    metadata = {"backend": backend, "ase_version": version("ase")}
    if backend in {"ani", "aimnet_2"}:
        metadata["torch_version"] = version("torch")
    if backend == "ani":
        metadata.update(
            torchani_version=version("torchani"), model=request["settings"]["model"]
        )
    if backend == "aimnet_2":
        from pyar.backends.aimnet2_assets import MODEL_PATH

        metadata.update(
            model=Path(MODEL_PATH).name,
            model_sha256=hashlib.sha256(Path(MODEL_PATH).read_bytes()).hexdigest(),
        )
    return metadata


def main():
    from ase import Atoms
    from ase.units import Hartree, Bohr

    request = json.loads(Path("worker_request.json").read_text())
    atoms = Atoms(request["symbols"], positions=request["coordinates_angstrom"])
    atoms.calc = CALCULATORS[request["settings"]["software"]](request)
    energy = atoms.get_potential_energy() / Hartree
    gradient = None
    if request["gradient"]:
        gradient = (-np.asarray(atoms.get_forces()) * Bohr / Hartree).tolist()
    metadata = calculator_metadata(request)
    write_json(
        "worker_result.json",
        {
            "energy_hartree": energy,
            "gradient_hartree_per_bohr": gradient,
            "metadata": metadata,
        },
    )


if __name__ == "__main__":
    main()
