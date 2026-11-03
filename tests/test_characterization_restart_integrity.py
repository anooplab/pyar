"""Tampered restart evidence must be rejected read-only, before backend work."""

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from pyar.core.molecule import Molecule
from pyar.persistence import ArtifactError, ArtifactRef
from pyar.state.reaction_v2 import ReactionRunStateV2
from pyar.state.reaction import ReactionStateError
from pyar.workflows import reaction_characterization as workflow


def write_xyz(path, symbols=("H", "He")):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"2\nsynthetic boundary fixture\n{symbols[0]} 0 0 0\n{symbols[1]} 2 0 0\n"
    )
    return path


def snapshot(root):
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file()
    }


def edit_json(path, edit):
    data = json.loads(path.read_text())
    edit(data)
    path.write_text(json.dumps(data))


@pytest.fixture
def prepared_route(tmp_path, monkeypatch):
    reaction = tmp_path / "reaction"
    job = reaction / "gamma_0100/orientation_x/job_gamma_0100_x"
    candidates = job / "candidate_ts"
    product_path = write_xyz(reaction / "products/accepted.xyz")
    write_xyz(job.parent / "trial_gamma_0100_x.xyz")
    write_xyz(candidates / "highest_backend_energy.xyz")
    (candidates / "metadata.json").write_text("{}\n")
    (job / "reaction_trace").mkdir()
    (job / "reaction_trace/trace.jsonl").write_text("")
    request_data = {
        "backend_parameters": {"software": "xtb"},
        "gamma_schedule": [100.0],
    }
    molecule = Molecule.from_xyz(str(product_path))
    state = ReactionRunStateV2.create(tmp_path, request_data, [], [molecule])
    state.record_product(
        "gamma_0100_x",
        100.0,
        "fixture-identity",
        "fixture-smiles",
        product_path,
        trace_summary={"candidate_ts_directory": str(candidates)},
    )
    state.finish("completed_products_found")
    evidence_candidate = (
        reaction / state.data["products"][0]["trace_summary"]["candidate_ts_directory"]
    )
    request = SimpleNamespace(
        restart_request=request_data,
        qc_params={"software": "xtb", "charge": 1, "multiplicity": 1},
    )
    options = {"images": 5, "pathway_ts_source": "highest-backend-energy"}
    run = Mock(
        return_value={"status": "complete", "completed_stages": ["relax", "neb"]}
    )
    monkeypatch.setattr(workflow, "run_path_stages", run)
    result = workflow.characterize_reaction(tmp_path, request, "neb", options)
    assert result["status"] == "complete"
    directory = Path(result["selected_routes"][0]["directory"])
    artifact = directory / "stage-result.xyz"
    write_xyz(artifact)
    reference = ArtifactRef.capture(directory, artifact.name, "test-stage-result")
    commit = directory / "stages/neb/attempts/attempt_000/commit.json"
    commit.parent.mkdir(parents=True)
    commit.write_text(
        json.dumps({"status": "completed", "artifacts": [reference.to_dict()]})
    )
    run.reset_mock()
    return SimpleNamespace(
        root=tmp_path,
        reaction=reaction,
        job=job,
        evidence_job=evidence_candidate.parent,
        directory=directory,
        request=request,
        options=options,
        backend=run,
        commit=commit,
        artifact=artifact,
        state=reaction / "characterization/pathways.json",
    )


def validate_read_only(route, error=None, match=None):
    before = snapshot(route.root)
    if error is None:
        result = workflow.validate_characterization_restart(
            route.root, route.request, "neb", route.options
        )
        assert result == json.loads(route.state.read_text())
    else:
        with pytest.raises(error, match=match):
            workflow.validate_characterization_restart(
                route.root, route.request, "neb", route.options
            )
    assert snapshot(route.root) == before
    route.backend.assert_not_called()


def test_valid_restart_verifies_artifacts_without_backend_or_writes(prepared_route):
    validate_read_only(prepared_route)


@pytest.mark.parametrize(
    "field,message",
    [
        ("reaction_request", "does not match completed reaction discovery"),
        ("discovery_commit_sha256", "Completed reaction state changed"),
    ],
)
def test_changed_discovery_binding_is_rejected(prepared_route, field, message):
    edit_json(prepared_route.state, lambda data: data.update({field: "changed"}))
    validate_read_only(prepared_route, ValueError, message)


@pytest.mark.parametrize("owner", ["summary", "pathway", "invocation"])
def test_changed_stage_options_are_rejected_at_their_owner(prepared_route, owner):
    route = prepared_route
    if owner == "summary":
        edit_json(
            route.state,
            lambda data: data["routes"]["product_001"]["path_options"].update(images=7),
        )
    elif owner == "pathway":
        edit_json(
            route.directory / "pathway.json",
            lambda data: data["path_options"].update(images=7),
        )
    else:
        route.options["images"] = 7
    validate_read_only(route, ValueError, "settings differ")


@pytest.mark.parametrize("role", ["reactant", "product", "ts_guess"])
@pytest.mark.parametrize("mutation", ["missing", "corrupt"])
def test_prepared_input_geometry_cannot_be_replaced_or_removed(
    prepared_route, role, mutation
):
    artifact = prepared_route.directory / f"inputs/{role}.xyz"
    if mutation == "missing":
        artifact.unlink()
    else:
        artifact.write_text("corrupted geometry\n")
    validate_read_only(
        prepared_route, ValueError, "input geometry was changed or removed"
    )


def test_missing_route_identity_is_rejected(prepared_route):
    (prepared_route.directory / "pathway.json").unlink()
    validate_read_only(prepared_route, ValueError, "restart metadata is missing")


def test_unknown_prepared_input_role_is_rejected(prepared_route):
    edit_json(
        prepared_route.directory / "pathway.json",
        lambda data: data["pathway_input_files"].update(unknown="inputs/reactant.xyz"),
    )
    validate_read_only(
        prepared_route, ValueError, "Unknown pathway input geometry role"
    )


def test_incomplete_stage_commit_cannot_authorize_reuse(prepared_route):
    edit_json(prepared_route.commit, lambda data: data.update(status="failed"))
    validate_read_only(
        prepared_route, ValueError, "Invalid characterization stage commit"
    )


@pytest.mark.parametrize("mutation", ["missing", "corrupt"])
def test_stage_artifact_hash_is_restart_authority(prepared_route, mutation):
    if mutation == "missing":
        prepared_route.artifact.unlink()
    else:
        prepared_route.artifact.write_text("tampered\n")
    validate_read_only(prepared_route, ArtifactError)


@pytest.mark.parametrize(
    "relative",
    [
        "candidate_ts/metadata.json",
        "reaction_trace/trace.jsonl",
        "candidate_ts/highest_backend_energy.xyz",
    ],
)
def test_mutable_job_files_do_not_replace_committed_discovery_evidence(
    prepared_route, relative
):
    path = prepared_route.job / relative
    path.write_text(path.read_text() + "\n")
    validate_read_only(prepared_route)


@pytest.mark.parametrize(
    "relative",
    [
        "candidate_ts/metadata.json",
        "reaction_trace/trace.jsonl",
        "candidate_ts/highest_backend_energy.xyz",
    ],
)
def test_committed_discovery_evidence_is_verified_before_route_reuse(
    prepared_route, relative
):
    path = prepared_route.evidence_job / relative
    path.write_text(path.read_text() + "\n")
    validate_read_only(prepared_route, ReactionStateError, "checksum mismatch")


@pytest.mark.parametrize(
    "key",
    [
        "reactant_sha256",
        "product_sha256",
        "ts_guess_sha256",
        "reaction_trace_sha256",
        "candidate_metadata_sha256",
    ],
)
def test_route_input_hash_must_match_verified_discovery_sources(prepared_route, key):
    edit_json(
        prepared_route.state,
        lambda data: data["routes"]["product_001"]["inputs"].update({key: "0" * 64}),
    )
    validate_read_only(prepared_route, ValueError, "Reaction evidence.*changed")


@pytest.mark.parametrize("role", ["reactant", "waypoint"])
def test_committed_atom_order_cannot_be_changed_before_a_new_route(
    prepared_route, role
):
    edit_json(prepared_route.state, lambda data: data.update(routes={}))
    path = (
        prepared_route.evidence_job.parent / "trial_gamma_0100_x.xyz"
        if role == "reactant"
        else prepared_route.evidence_job / "candidate_ts/highest_backend_energy.xyz"
    )
    write_xyz(path, symbols=("He", "H"))
    validate_read_only(prepared_route, ReactionStateError, "checksum mismatch")


def test_pathway_directory_without_restart_state_is_rejected(prepared_route):
    prepared_route.state.unlink()
    validate_read_only(prepared_route, ValueError, "has no restart state")


def test_continuation_without_discovery_creates_no_files(tmp_path):
    request = SimpleNamespace(restart_request={})
    with pytest.raises(ValueError, match="requires completed reaction discovery"):
        workflow.validate_characterization_restart(tmp_path, request, "neb", {})
    assert snapshot(tmp_path) == {}
