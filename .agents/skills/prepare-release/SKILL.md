---
name: prepare-release
description: Prepare PyAR release-candidate evidence and evaluate mandatory engineering gates when release preparation or readiness review is requested. Preparation alone does not authorize tags or package publication.
---

# Prepare a PyAR release

Read repository `AGENTS.md` and the canonical
[release gates](../../../docs/release_engineering/gates.md). Consult the
[audit](../../../docs/release_engineering/audit.md) for existing infrastructure
and the [staged plan](../../../docs/release_engineering/plan.md) for missing
automation. These are release-engineering stages, separate from the ORE series.
Do not implement planned infrastructure merely because this skill was invoked.

## Establish the candidate

1. Identify the requested version, target registry, intended commit and approved
   scientific/backend support matrix. Inspect `git status -sb`, `git rev-parse
   HEAD`, `pyproject.toml`, `CHANGELOG.md`, `MANIFEST.in` and both GitHub workflows.
   Preserve unrelated edits; use a separate clean candidate checkout when needed.
   Preparation does not silently bump the version or update scientific methods.
2. Choose a fresh evidence directory outside the checkout. Record the interpreter,
   dependency/build-tool versions, candidate SHA and exact commands/exits. A
   provisioned test environment may lack pip or release tools: report that fact
   or create an authorized isolated environment; never mutate a user's shared
   environment implicitly. Use pip with dependency-group support and the current
   metadata/CI install commands, not nonexistent `[test,docs]` extras.
3. Inspect the available candidate-manifest and gate validators. Where the plan
   still lists a required validator/corpus as missing, mark that gate BLOCKED.
   Collect available evidence without replacing the missing check with an AI
   assessment. Never execute a guessed script path or record an unrun command as
   successful. Require a frozen case manifest, tolerances, assets and budget
   before running scientific campaigns.

## Execute deterministic checks

Use the existing CI jobs and repository tools, retaining their structured/raw
outputs. `EVIDENCE` below must be a new absolute directory outside the checkout,
with the candidate checkout as the starting working directory. Run commands with
normal Python; repository archive/smoke assertions must not be disabled with `-O`.

```bash
python -m pytest --cov=pyar --cov-report=term-missing --junitxml="$EVIDENCE/pytest.xml"
python -m sphinx -E -b html -W --keep-going docs "$EVIDENCE/docs"
python -m build --outdir "$EVIDENCE/dist"
python tools/check_distributions.py "$EVIDENCE/dist"
python -m twine check "$EVIDENCE"/dist/*
sha256sum "$EVIDENCE"/dist/*
```

- Use the exact candidate environment/extras from `ci.yml`; provision optional
  Sella/geodesic environments where required. Inspect coverage and every
  skip/xfail against the approved inventory. Log command failures as failures;
  no success from the final command of a failed multi-command sequence.
- Sphinx rewrites `docs/generated_api.rst`. Review the diff and restore only
  generator changes introduced by this check, preserving prior user edits.
- Use build isolation for qualification. A diagnosed offline fallback may use
  recorded installed build tools with `--no-isolation`, but it does not by
  itself satisfy isolated-build/dependency-resolution evidence.
- Install the exact wheel and sdist into separate fresh virtual environments
  without inherited system packages; execute outside checkout with checkout
  `PYTHONPATH` entries removed. Invoke each installed interpreter's `-m pip
  check` and the absolute candidate path to `tools/installation_smoke.py core`.
  In separate appropriate environments invoke its supported workflow modes
  required by the support matrix. Follow the existing installation/archive CI
  jobs for commands; do not use source/editable imports as distribution proof.
  Confirm import origins and installed metadata match the approved artifact.
- Scientific commands must come from the approved case manifest and existing
  benchmark interfaces. `benchmarks/scan_bond/validate_backend_scans.py` is an
  available small scan smoke, not a replacement for the scientific release
  corpus. See the existing TS benchmark README and conformer/microsolvation
  runners for their actual interfaces. Require the numerical comparator's exit
  and structured verdict, raw backend outputs, versions, reference hashes,
  objective/identity/pathway checks and complete cost accounting. An absent
  comparator, unavailable backend/model or missing case is BLOCKED/NOT_RUN.

Set `CANDIDATE_SHA` to the frozen checkout's `git rev-parse HEAD` output.
Select `CI_RUN_ID` from the matching `ci.yml` run list, then verify its
`headSha` before accepting any jobs. Inspect exact-SHA remote evidence
with existing read-only GitHub interfaces:

```bash
gh run list --workflow ci.yml --commit "$CANDIDATE_SHA"
gh run view "$CI_RUN_ID" --json headSha,status,conclusion,jobs
gh api "repos/anooplab/pyar/commits/$CANDIDATE_SHA/check-runs" --paginate
```

Confirm head SHA, required job identities, completion and conclusions against
the release matrix; do not hard-code a historical check count. A previous green
run or a local subset is not a complete matrix. Inspect protected-environment
and repository settings read-only; an inaccessible setting is unknown. Manual
dispatch of `publish.yml` is publication, not a readiness check.

## Produce evidence and stop at approval

Report each gate using the states defined in `gates.md`, with exact candidate,
artifact hashes, command exits, CI run IDs, case results, exclusions and evidence
locations. Missing required checks or failed scientific assertions yield NO-GO.
Only complete G0–G5 evidence yields READY FOR OWNER APPROVAL. Record software,
scientific and installation evidence separately; no general chemistry claim
from mocks, process exit, an import or a small smoke calculation.

Preparation alone stops before creating/pushing a release tag, dispatching the
publish workflow or uploading packages. Honor existing explicit authorization
only when it covers the exact version, candidate, artifacts and registry; a
generic commit/merge instruction does not. A candidate or artifact change
requires renewed applicable evidence and approval. Current publishing gaps in
the audit are blockers, not reasons to bypass the gates manually.

For a separately authorized publication task, require G6 and the implemented
protected promotion procedure, then G7's deterministic version-specific registry
hash checks and fresh registry installs. Do not treat this paragraph as a
publish command or authority to implement R5. Until its verifier exists, report
post-release verification as missing rather than claiming a successful release.
