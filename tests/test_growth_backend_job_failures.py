import json

from pyar.backend_errors import BackendExecutionError
from pyar.growth import service


class DummyMolecule:
    def __init__(self, name):
        self.name = name


def test_growth_batch_records_failure_and_attempts_remaining_jobs(
    tmp_path, monkeypatch
):
    molecules = [DummyMolecule("failed_candidate"), DummyMolecule("good_candidate")]
    attempted = []

    def optimize(molecule, _params, **_kwargs):
        attempted.append(molecule.name)
        if molecule.name == "failed_candidate":
            raise BackendExecutionError("backend rejected candidate")
        return True, None, None

    monkeypatch.setattr(service, "optimize_workflow_once", optimize)
    failure_file = tmp_path / "backend_failures.json"

    statuses = service._optimize_growth_batch(
        molecules, {"software": "xtb"}, "growth_preselection", str(failure_file)
    )

    assert attempted == ["failed_candidate", "good_candidate"]
    assert statuses == [False, True]
    failure = json.loads(failure_file.read_text())[0]
    assert failure["molecule"] == "failed_candidate"
    assert "backend rejected candidate" in failure["error"]
