#!/usr/bin/env python3
"""TorchANI-based ANI backend and CLI wrapper."""

import argparse
import logging

from pyar.optional_dependencies import optional_dependency_error

try:
    import ase.io
    import ase.optimize
    from ase.calculators.calculator import Calculator, all_changes
    import torch
except ImportError as exc:
    module_name = getattr(exc, "name", None) or "torchani"
    raise optional_dependency_error(module_name, feature="ANI backend", extra="ml") from exc

logger = logging.getLogger("pyar.backends.ani")

_MODEL_ALIASES = {
    "ANI-1x": "ANI1x",
    "ANI-1ccx": "ANI1ccx",
    "ANI-2x": "ANI2x",
}

try:
    import torchani
except ImportError:  # pragma: no cover - optional runtime dependency
    torchani = None


class ANICalculationFailed(Exception):
    """Raised when an ANI model lookup, evaluation, or workflow step fails."""


def _require_torchani():
    """Return the optional TorchANI import or raise a clear import error."""
    if torchani is None:
        raise optional_dependency_error("torchani", feature="ANI backend", extra="ml")
    return torchani


def _load_model(model_name):
    """Load a TorchANI model by name, honoring PyAR compatibility aliases."""
    torchani_mod = _require_torchani()
    normalized_name = _MODEL_ALIASES.get(model_name, model_name)
    try:
        return torchani_mod.models.__dict__[normalized_name]()
    except KeyError as exc:
        msg = f"Model '{model_name}' not found in torchani models."
        logger.error(msg)
        raise ANICalculationFailed(msg) from exc


class ANI(Calculator):
    """ASE calculator backed by a TorchANI potential."""

    implemented_properties = ["energy"]

    def __init__(self, species, model="ANI-1x"):
        """Create a calculator for a given species list and ANI model."""
        super().__init__()
        self.species = species
        self.model_name = model
        self.model = _load_model(model)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        """Compute the potential energy for an ASE atoms object."""
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ANICalculationFailed("Atoms are required for an ANI calculation")

        try:
            species = torch.as_tensor(
                [atoms.get_atomic_numbers()],
                dtype=torch.long,
            )
            positions = torch.as_tensor(
                [atoms.get_positions()],
                dtype=torch.float32,
            )
            energy = self.model((species, positions)).energies.squeeze(0).item()
        except Exception as exc:
            msg = f"ANI model calculation failed: {exc}"
            logger.error(msg)
            raise ANICalculationFailed(msg) from exc

        self.results = {"energy": energy}


class ANIInterface:
    """Small convenience wrapper for single-point and optimization jobs."""

    def __init__(self, xyzfile):
        """Load an XYZ file into an ASE atoms object."""
        try:
            self.atoms = ase.io.read(xyzfile)
        except OSError as exc:
            msg = f"Could not read XYZ file {xyzfile}"
            logger.error(msg)
            raise ANICalculationFailed(msg) from exc

    def _attach_calculator(self, model):
        """Attach the requested ANI calculator to the loaded structure."""
        self.atoms.calc = ANI(self.atoms.get_chemical_symbols(), model=model)

    def optimize(self, model="ANI-1x", fmax=1e-3, trajectory="optimization.traj"):
        """Run an ASE geometry optimization and return the final energy."""
        try:
            self._attach_calculator(model)
            dyn = ase.optimize.BFGS(self.atoms, trajectory=trajectory)
            dyn.run(fmax=fmax)
            return self.atoms.get_potential_energy()
        except ANICalculationFailed:
            raise
        except Exception as exc:
            msg = f"Geometry optimization failed: {exc}"
            logger.error(msg)
            raise ANICalculationFailed(msg) from exc

    def singlepoint(self, model="ANI-1x"):
        """Return a single-point ANI energy for the loaded structure."""
        try:
            self._attach_calculator(model)
            return self.atoms.get_potential_energy()
        except ANICalculationFailed:
            raise
        except Exception as exc:
            msg = f"Single-point calculation failed: {exc}"
            logger.error(msg)
            raise ANICalculationFailed(msg) from exc


def main():
    """Run the ANI helper CLI."""
    parser = argparse.ArgumentParser(description="TorchANI helper for PyAR")
    parser.add_argument("command", choices=["optimize", "singlepoint"])
    parser.add_argument("xyzfile")
    parser.add_argument("--model", default="ANI-1x")
    parser.add_argument("--fmax", type=float, default=1e-3)
    parser.add_argument("--trajectory", default="optimization.traj")
    args = parser.parse_args()

    interface = ANIInterface(args.xyzfile)
    if args.command == "optimize":
        energy = interface.optimize(
            model=args.model,
            fmax=args.fmax,
            trajectory=args.trajectory,
        )
    else:
        energy = interface.singlepoint(model=args.model)

    print(energy)


if __name__ == "__main__":
    main()
