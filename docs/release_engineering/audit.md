# Release engineering R0 audit

Date: 2026-10-18. Baseline: `master`,
`9afd8b23ac8189a16fa54fe2bbe907a79bf76b71`.
The checkout was clean and aligned with `origin/master` before this audit.
Deliverables are policies and a preparation skill only. R1–R5 automation,
scientific algorithms, public APIs, CI configuration and packaging metadata are
unchanged; no release, tag, version bump or publishing dispatch was initiated.

## Instructions and skill reconciliation

- Preserve AGENTS.md's source layout, style, testing, scoped commit/PR and
  external-backend instructions. Two stale details need alignment with current
  metadata: test/docs tools are dependency **groups**, not the nonexistent
  `[test,docs]` extras; `setup.py` is absent and must not become a second
  metadata source. The updated development install includes the selection/xTB
  capabilities used by CI plus the `dev` group.
- Add six concise permanent engineering policies and link the canonical
  [release gates](gates.md), rather than duplicating the gate matrix in AGENTS.
- `.agents/` and `.codex/` had no repository skill files at baseline. The new
  [prepare-release skill](../../.agents/skills/prepare-release/SKILL.md) is the
  portable repository entry point. Existing external `pyar-pr-delivery` and
  historical local `pyar-release-checklist` guidance was reviewed: preserve
  their isolation, exact-head evidence, archive inspection, original-work
  preservation and explicit publication-authorization rules. Do not copy their
  incidental past counts or assume their installed availability for collaborators.
- The skill uses existing scripts and CI jobs. Missing deterministic release
  validators block readiness; it does not invent results or implement scripts
  as a side effect. Release gates apply to future publication, separately from
  the ORE feature scope and its frozen scientific amendment.

## Infrastructure and gaps

| Area | Existing evidence / infrastructure | Missing component and priority |
| --- | --- | --- |
| CI/regressions | `ci.yml`: Python 3.10–3.13 test matrix, OpenBabel, `pip check`, pytest-cov and legacy CLI help; coverage configuration has `fail_under = 50`. Dedicated Python 3.11 Sella and geodesic jobs run `tests/test_neb.py`. | **P0/P1, R1:** candidate-bound aggregate required-check/evidence validation; controlled skip/xfail inventory and artifacts/retention. The optional jobs install editable source despite step wording referring to a released package. No dedicated lint/type job is configured; consider scoped static checks without wholesale cleanup. |
| Known scientific failures | Recovery audit tests explicitly retain strict expected failures for native convergence, energy/geometry validation, stale markers and cache provenance. Negative/replay/path-gate tests and fixed fixtures exist. | **P0, R3:** issue-owned release disposition for safety xfails; green pytest including xfails is not proof those hazards are fixed. No release-wide scientific baseline/comparator or frozen quantitative acceptance manifest. |
| Scientific runners | `benchmarks/scan_bond/validate_backend_scans.py` produces structured smoke results and a failing exit for unsuccessful scans. Conformer and TS benchmark CLIs, microsolvation runners, fixtures, provenance and dated analyses exist. | **P0, R3:** connect the appropriate runners to a licensed, versioned representative release corpus and deterministic verdicts; record candidate/reference hashes, tolerances, seeds, costs and backend/model versions. The scan script explicitly is not a TS benchmark. Historical reports do not qualify a new candidate. |
| Optional backends | Capability/dependency registries; missing-extra/executable tests; installation matrix covers core, aggregate, conformer, reaction, AIMNet2 composition, NEB, ML, analysis and microsolvation on Python 3.10 and 3.13. Sella/geodesic have separate jobs. | **P1, R1/R3/R4:** explicit supported/excluded backend/version/platform matrix, final-artifact tests for remaining extras/aliases/compositions, real licensed and model-enabled runs. Import/discovery and failure-path evidence must remain distinct from numerical chemistry validation. |
| Packaging | `pyproject.toml` is authoritative: setuptools, `pyar-chem` 1.5.1, `pyar` namespace, Python >=3.10, runtime extras and test/docs/release/dev groups. `MANIFEST.in` excludes model weights and selected benchmark reruns; license/provenance files exist. Metadata and dependency-contract tests are present. | **P1, R1/R4:** unified version/changelog/import/archive consistency check, recorded constraints, fuller license/required-content checks and protection from stale local archives. Python >=3.10 is broader than the 3.10–3.13 CI/classifier evidence; do not silently claim a newer Python is release-validated. |
| Installation/artifacts | `tools/check_distributions.py` inspects wheel and sdist contents. `tools/installation_smoke.py` checks installed imports/CLI/doctor/dependency behavior outside checkout. Archives job builds and installs fresh wheel/sdist on Python 3.11; wheel checks core, aggregate and reaction. | **P1, R1/R4:** expand exact final-artifact coverage across the declared support matrix; dependency snapshots, candidate hashes and a single evidence record. Source-install matrix is not equivalent to final-wheel matrix. Existing archive checks use assertions and must run with normal Python, not `python -O`. |
| Documentation | Strict Sphinx CI and ReadTheDocs configuration; generated API hook; installation/backend guides, migration notes, dated packaging/layout audits and changelog. | **P1, R2:** executable example/CLI contracts and systematic synchronization. Changelog `Unreleased` does not describe all current recovery/intramolecular additions. Sphinx writes tracked `docs/generated_api.rst`; establish generated-output ownership. Keep historical records dated, not treated as current guarantees. |
| Publishing | `publish.yml` supports `v*` tag pushes and manual dispatch, builds distributions, checks version only for tag refs, uploads a `dist` artifact and uses OIDC in the `pypi` environment. | **P0, R1/R5:** publish depends only on `build`, not CI/science/install/docs gates. Manual dispatch skips tag/version verification. No approved-candidate/hash check or promotion of the already tested CI artifacts. Current environment has no protection rules or deployment restriction. |
| Approval/post-release | Historical release practice includes explicit authorization and version-specific registry verification. | **P0/P1, R4/R5:** enforce approval bound to version/SHA/artifact hashes, verify repository/environment settings, retain signed-off evidence, automate registry hash/install checks and document incidents. Preparation must not trigger existing publishing while these gaps remain. |

## Execution evidence collected in this audit

These are observations at the audited baseline, not a release-readiness verdict
for the new uncommitted documentation or a future package:

| Check actually performed | Result / limit |
| --- | --- |
| `git status -sb`, `git rev-parse HEAD`, repository file inspection | Clean baseline `9afd8b2`; no existing repository-local skill files. |
| `gh api repos/anooplab/pyar/commits/9afd8b23ac8189a16fa54fe2bbe907a79bf76b71/check-runs` | 29 check runs, all completed/success at query time: 26 configured CI jobs and three Pages-related checks. [CI run](https://github.com/anooplab/pyar/actions/runs/39136962382). This is observed success, not a hard-coded future check count or a scientific benchmark verdict. |
| `gh api repos/anooplab/pyar/environments/pypi` | Returned `protection_rules: []`, `deployment_branch_policy: null`. Naming an environment does not supply a required reviewer. |
| `gh api repos/anooplab/pyar/branches/master/protection` | HTTP 404. Classic branch protection could not be established through this endpoint; absence, token visibility and separate rulesets are not distinguished by this observation. Repository-level enforcement needs administrator verification. |
| Provisioned interpreter: `python -m pytest -q tests/test_packaging_metadata.py tests/test_installation_contract.py -k 'not built_wheel_excludes_large_model_assets'` | **29 passed, 1 deselected**, exit 0. Deselect the in-tree wheel check to avoid crediting an existing `dist/` artifact to this candidate. No new distribution was built or installed during this audit. |
| Provisioned interpreter: `python -m pip --version` | Exit 1: `No module named pip`. The provisioned test environment is not itself a clean release-tooling environment. No dependencies were installed or environment repaired here. |

The full test suite, scientific benchmark campaign, clean-install matrix,
distribution build, Twine checks and post-release registry checks were **not
executed locally as part of this audit**. CI evidence above is explicitly remote
and candidate-bound. Instruction-file validation is recorded in the delivery
report after the skill is written; it does not satisfy G0–G7 by itself.

R0 deliverable verification: the bundled skill-creator `quick_validate.py`
passed against the repository skill using system Python 3. The provisioned test
interpreter lacked PyYAML for this validator, so no dependencies were installed
and the already available system interpreter was used. A one-off local check
resolved all **11 local Markdown links** across the five deliverables and
confirmed the changed-path set contains only those five files. `git diff
--check` passed. No release orchestrator or gate script was added to obtain
these instruction-file checks.

## Implementation priorities and disposition

Start with candidate evidence and publishing safety design, then executable
documentation checks and scientific/backend qualification, then exact-artifact
RC testing and protected promotion/post-release verification. See [plan.md](plan.md)
for R0–R5 ownership, scope and exit criteria. Missing gate automation and
scientific evidence mean **release readiness is not established**.

R0 intentionally leaves CI, tools, algorithms, packaging, existing scientific
plans and public APIs untouched. These gaps are implementation work for future
authorized stages, not findings silently repaired by this instruction change.
