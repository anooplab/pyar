import hashlib
import urllib.error
from contextlib import nullcontext

import pytest

from tools.release_promotion import PromotionError, assert_unpublished, validate_release
from tools.verify_pypi_release import RegistryError, validate_pypi_json

pytestmark = [pytest.mark.unit, pytest.mark.regression]


SHA = "a" * 40
WHEEL = "pyar_chem-1.5.2-py3-none-any.whl"
SDIST = "pyar_chem-1.5.2.tar.gz"


def run(run_id, workflow):
    name = (
        "R1 required checks"
        if workflow.endswith("/ci.yml")
        else "Real xTB GFN2 scan smoke"
    )
    return {
        "id": run_id,
        "path": workflow,
        "head_sha": SHA,
        "status": "completed",
        "conclusion": "success",
        "jobs": [{"name": name, "conclusion": "success"}],
    }


def test_release_promotion_requires_exact_approval_runs_and_bytes(tmp_path):
    wheel = tmp_path / WHEEL
    sdist = tmp_path / SDIST
    wheel.write_bytes(b"approved wheel")
    sdist.write_bytes(b"approved source")
    approval = {
        "schema_version": 1,
        "status": "APPROVED",
        "version": "1.5.2",
        "candidate_sha": SHA,
        "target_registry": "https://upload.pypi.org/legacy/",
        "approved_by": "project owner",
        "approved_at_utc": "2026-10-19T00:00:00Z",
        "scientific_scope_status": "PASS_SCOPED",
        "scientific_exclusions": ["reaction pathways"],
        "workflow_runs": {"ci": "1", "tier2": "2"},
        "required_jobs": {
            "ci": ["R1 required checks"],
            "tier2": ["Real xTB GFN2 scan smoke"],
        },
        "artifacts": {
            WHEEL: hashlib.sha256(wheel.read_bytes()).hexdigest(),
            SDIST: hashlib.sha256(sdist.read_bytes()).hexdigest(),
        },
    }
    hashes = validate_release(
        approval,
        tmp_path,
        version="1.5.2",
        candidate_sha=SHA,
        ci_run=run("1", ".github/workflows/ci.yml"),
        tier2_run=run("2", ".github/workflows/scientific-tier2.yml"),
    )
    assert hashes == approval["artifacts"]


@pytest.mark.parametrize("field", ["candidate_sha", "target_registry", "status"])
def test_release_promotion_rejects_approval_mismatch(tmp_path, field):
    approval = {
        "schema_version": 1,
        "status": "APPROVED",
        "version": "1.5.2",
        "candidate_sha": SHA,
        "target_registry": "https://upload.pypi.org/legacy/",
        "approved_by": "project owner",
        "approved_at_utc": "2026-10-19T00:00:00Z",
        "scientific_scope_status": "PASS_SCOPED",
        "scientific_exclusions": ["reaction pathways"],
        "workflow_runs": {"ci": "1", "tier2": "2"},
        "required_jobs": {
            "ci": ["R1 required checks"],
            "tier2": ["Real xTB GFN2 scan smoke"],
        },
        "artifacts": {WHEEL: "0" * 64, SDIST: "1" * 64},
    }
    approval[field] = "wrong"
    with pytest.raises(PromotionError):
        validate_release(
            approval,
            tmp_path,
            version="1.5.2",
            candidate_sha=SHA,
            ci_run=run("1", ".github/workflows/ci.yml"),
            tier2_run=run("2", ".github/workflows/scientific-tier2.yml"),
        )


def test_release_promotion_rejects_run_or_artifact_drift(tmp_path):
    (tmp_path / WHEEL).write_bytes(b"wheel")
    (tmp_path / SDIST).write_bytes(b"sdist")
    approval = {
        "schema_version": 1,
        "status": "APPROVED",
        "version": "1.5.2",
        "candidate_sha": SHA,
        "target_registry": "https://upload.pypi.org/legacy/",
        "approved_by": "project owner",
        "approved_at_utc": "2026-10-19T00:00:00Z",
        "scientific_scope_status": "PASS_SCOPED",
        "scientific_exclusions": ["reaction pathways"],
        "workflow_runs": {"ci": "1", "tier2": "2"},
        "required_jobs": {
            "ci": ["R1 required checks"],
            "tier2": ["Real xTB GFN2 scan smoke"],
        },
        "artifacts": {WHEEL: "0" * 64, SDIST: "1" * 64},
    }
    wrong_run = run("999", ".github/workflows/scientific-tier2.yml")
    with pytest.raises(PromotionError, match="run id"):
        validate_release(
            approval,
            tmp_path,
            version="1.5.2",
            candidate_sha=SHA,
            ci_run=run("1", ".github/workflows/ci.yml"),
            tier2_run=wrong_run,
        )
    with pytest.raises(PromotionError, match="hashes"):
        validate_release(
            approval,
            tmp_path,
            version="1.5.2",
            candidate_sha=SHA,
            ci_run=run("1", ".github/workflows/ci.yml"),
            tier2_run=run("2", ".github/workflows/scientific-tier2.yml"),
        )


def test_pypi_verification_requires_exact_file_hashes():
    expected = {WHEEL: "a" * 64, SDIST: "b" * 64}
    document = {
        "info": {"version": "1.5.2"},
        "releases": {
            "1.5.2": [
                {"filename": name, "digests": {"sha256": digest}}
                for name, digest in expected.items()
            ]
        },
    }
    assert validate_pypi_json(document, "1.5.2", expected) == expected
    document["releases"]["1.5.2"][0]["digests"]["sha256"] = "c" * 64
    with pytest.raises(RegistryError, match="hashes differ"):
        validate_pypi_json(document, "1.5.2", expected)


def test_pypi_version_endpoint_urls_are_verified():
    expected = {WHEEL: "a" * 64, SDIST: "b" * 64}
    document = {
        "info": {"version": "1.5.2"},
        "urls": [
            {"filename": name, "digests": {"sha256": digest}}
            for name, digest in expected.items()
        ],
    }
    assert validate_pypi_json(document, "1.5.2", expected) == expected


def test_pypi_verification_rejects_duplicate_release_filenames():
    expected = {WHEEL: "a" * 64, SDIST: "b" * 64}
    item = {"filename": WHEEL, "digests": {"sha256": expected[WHEEL]}}
    document = {"info": {"version": "1.5.2"}, "releases": {"1.5.2": [item, item]}}
    with pytest.raises(RegistryError, match="duplicate filenames"):
        validate_pypi_json(document, "1.5.2", expected)


def test_assert_unpublished_only_accepts_a_pypi_404(monkeypatch):
    import tools.release_promotion as promotion

    monkeypatch.setattr(
        promotion.urllib.request, "urlopen", lambda *a, **k: nullcontext()
    )
    with pytest.raises(PromotionError, match="already exists"):
        assert_unpublished("1.5.2")

    def missing(*args, **kwargs):
        raise urllib.error.HTTPError("url", 404, "not found", {}, None)

    monkeypatch.setattr(promotion.urllib.request, "urlopen", missing)
    assert_unpublished("1.5.2")

    def unavailable(*args, **kwargs):
        raise urllib.error.HTTPError("url", 503, "unavailable", {}, None)

    monkeypatch.setattr(promotion.urllib.request, "urlopen", unavailable)
    with pytest.raises(PromotionError, match="HTTP 503"):
        assert_unpublished("1.5.2")
