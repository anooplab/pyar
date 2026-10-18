import pytest

from pyar.backend_errors import BackendExecutionError
from pyar import optimiser


class DummyMolecule:
    def __init__(self, name):
        self.name = name
        self.charge = 0
        self.multiplicity = 1
        self.scftype = "rhf"
        self.energy = None
        self.coordinates = None


def test_bulk_optimization_attempts_remaining_jobs_before_reporting_failure(
    monkeypatch,
):
    molecules = [DummyMolecule("bad"), DummyMolecule("good")]
    attempted = []

    def optimize(molecule, _params):
        attempted.append(molecule.name)
        if molecule.name == "bad":
            raise BackendExecutionError("backend rejected this geometry")
        molecule.energy = -1.0
        molecule.coordinates = [[0.0, 0.0, 0.0]]
        return True

    monkeypatch.setattr(optimiser, "optimise", optimize)

    with pytest.raises(BackendExecutionError, match="all input jobs were attempted"):
        optimiser.bulk_optimize(molecules, {"software": "xtb"})

    assert attempted == ["bad", "good"]
    assert molecules[1].energy == -1.0
