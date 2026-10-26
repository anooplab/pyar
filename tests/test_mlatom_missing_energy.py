import ast
from pathlib import Path
import tempfile
import types

import ase
import numpy as np
import pytest


def _load_calculator_class():
    """Compile the exact adapter class with narrow MLatom/ASE boundary stubs."""
    source = Path("pyar/backends/ASE.py")
    tree = ast.parse(source.read_text(encoding="utf-8"))
    calculator_node = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "MLatomCalculator"
    )

    class Calculator:
        def __init__(self, *args, **kwargs):
            self.results = {}

        def calculate(self, atoms=None, *args, **kwargs):
            self.atoms = atoms

    class Molecule:
        def read_from_xyz_file(self, **kwargs):
            self.xyz_coordinates = [[0.0, 0.0, 0.0]]

        def copy(self):
            return Molecule()

        def get_energy_gradients(self):
            return np.array([[0.1, -0.2, 0.3]])

    namespace = {
        "Calculator": Calculator,
        "all_changes": (),
        "ase": ase,
        "tempfile": tempfile,
        "io": types.SimpleNamespace(write=lambda *args, **kwargs: None),
        "initial_molecule": Molecule(),
        "data": types.SimpleNamespace(molecule=Molecule),
    }
    exec(
        compile(
            ast.Module(body=[calculator_node], type_ignores=[]), str(source), "exec"
        ),
        namespace,
    )
    return namespace["MLatomCalculator"], Molecule


def test_missing_mlatom_energy_raises_actionable_value_error():
    calculator_type, _molecule_type = _load_calculator_class()
    calculator = calculator_type(types.SimpleNamespace(predict=lambda **kwargs: None))

    with pytest.raises(ValueError, match="model did not return any energy"):
        calculator.calculate(atoms=None)


def test_valid_mlatom_energy_and_gradient_keep_existing_units_and_sign():
    calculator_type, _molecule_type = _load_calculator_class()

    def predict(*, molecule, **kwargs):
        molecule.energy = -0.5

    calculator = calculator_type(types.SimpleNamespace(predict=predict))

    calculator.calculate(atoms=None, properties=["energy", "forces"])

    assert calculator.results["energy"] == -0.5 * ase.units.Hartree
    np.testing.assert_array_equal(
        calculator.results["forces"],
        np.array([[-0.1, 0.2, -0.3]]) * ase.units.Hartree,
    )
