from pathlib import Path

import numpy as np
import pytest

from pyar.backend_errors import BackendExecutionError
from pyar.workflows import reaction


class DummyMolecule:
    def __init__(self, name):
        self.name = name
        self.energy = 0.0
        self.coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        self.atoms_list = ["H", "H"]
        self.fragments = [[0], [1]]

    def mol_to_xyz(self, filename):
        Path(filename).write_text("2\ninput\nH 0 0 0\nH 1 0 0\n")

    def copy(self):
        copy = DummyMolecule(self.name)
        copy.energy = self.energy
        copy.coordinates = self.coordinates.copy()
        return copy

    def is_bonded(self):
        return False


class DummyRunState:
    data = {"reactant_identity": None}

    def __init__(self):
        self.jobs = []

    def current_survivor_molecules(self):
        return []

    def record_job(self, *args, **kwargs):
        self.jobs.append((args, kwargs))


def test_failed_orientation_is_recorded_and_remaining_jobs_continue(
    tmp_path, monkeypatch
):
    state = DummyRunState()
    product_dir = tmp_path / "products"
    product_dir.mkdir()
    outcomes = iter([BackendExecutionError("backend rejected this geometry"), True])

    def optimize(*_args):
        outcome = next(outcomes)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(reaction, "write_disconnected_reference", lambda *_args: None)
    monkeypatch.setattr(
        reaction,
        "molecule_identity_from_xyz",
        lambda *_args: {"inchi": "reactants", "smiles": "reactants"},
    )
    monkeypatch.setattr(
        reaction,
        "_optimize_with_diagnostics",
        optimize,
    )

    survivors = reaction.optimize_all(
        gamma_id="0001",
        orientations=[
            DummyMolecule("orientation_00000001"),
            DummyMolecule("orientation_00000002"),
        ],
        run_state=state,
        product_dir=str(product_dir),
        qc_param={"gamma": 1.0},
    )

    assert len(survivors) == 1
    assert [job[0][2] for job in state.jobs] == ["failed_backend", True]
    assert "backend rejected this geometry" in state.jobs[0][1]["failure"]
    assert state.jobs[0][0][3] == []


def test_all_failed_orientations_still_raise_backend_failure(tmp_path, monkeypatch):
    state = DummyRunState()
    product_dir = tmp_path / "products"
    product_dir.mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(reaction, "write_disconnected_reference", lambda *_args: None)
    monkeypatch.setattr(
        reaction,
        "molecule_identity_from_xyz",
        lambda *_args: {"inchi": "reactants", "smiles": "reactants"},
    )

    def fail_job(*_args):
        raise BackendExecutionError("backend is unavailable")

    monkeypatch.setattr(reaction, "_optimize_with_diagnostics", fail_job)

    with pytest.raises(
        BackendExecutionError, match="All 1 reaction orientation backend jobs failed"
    ):
        reaction.optimize_all(
            gamma_id="0001",
            orientations=[DummyMolecule("orientation_00000001")],
            run_state=state,
            product_dir=str(product_dir),
            qc_param={"gamma": 1.0},
        )
    assert state.jobs[0][0][2] == "failed_backend"
    assert "backend is unavailable" in state.jobs[0][1]["failure"]
