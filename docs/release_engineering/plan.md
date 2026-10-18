# Staged release engineering plan

Status: R0 audit/instruction deliverables only, 2026-10-18.
R1–R5 below are **proposed, not implemented or authorized by this document**.
This is a release-engineering series, independent of the ORE milestone names.
No scientific algorithm, public API, package version, tag or publication changes
are part of R0. Use [gates.md](gates.md) as the sole release gate definition.

| Stage | Work when separately authorized | Deterministic exit evidence |
| --- | --- | --- |
| R0 Audit | Inventory instructions, skills, CI, packaging, docs, tests and benchmarks. Establish permanent policies, release preparation skill, gates and prioritized gaps without editing infrastructure. | Reviewable documents, valid skill metadata/local links, scoped diff and a dated audit distinguishing executed checks from historical/configured evidence. |
| R1 CI and packaging quality | Define candidate/evidence manifest schemas and required-check validation. Consolidate existing build/archive/install checks without losing matrix coverage; add consistent version, changelog, dependency, archive/license and skip/xfail checks. Record tooling/dependency constraints and candidate-bound artifacts. Specify publish enforcement interface without enabling publication yet. | Validator rejects wrong SHA/version, missing or stale jobs, omitted/failed checks and wrong artifact hashes. Fresh archives pass content/Twine/clean-install checks. Supported Python/extra matrix passes with retained logs. No gate can be satisfied by old local `dist/` files. |
| R2 Documentation cleanup and validation | Reconcile AGENTS, installation/developer commands, README, changelog, modern/legacy CLI and supported-backend/version claims. Inventory runnable examples; add deterministic help/preflight/smoke checks with declared fixtures. Decide ownership/output handling of generated API docs. Preserve historical audit records as dated evidence. | Strict Sphinx and example/CLI contract checks pass; deliberate option/version drift fails a check; documentation build side effects are isolated. No unsupported scientific claim or release/version claim remains in current user guidance. |
| R3 Scientific and backend validation | With scientific-owner approval, freeze the release corpus, reference hashes, numerical tolerances, seeds, hardware/thread policy, budgets and backend support/exclusion matrix. Connect existing benchmark runners to structured numerical assertions and a release verdict. Define authorized licensed/model-enabled runner execution and artifact retention. No change to methods or public APIs to fit results. | Repeated controlled runs produce acceptable results within declared tolerances. Comparator rejects wrong units/method/spin, nonconvergence, stale caches, bad constraints, identity changes and invalid saddle/IRC connections. Missing prerequisites produce BLOCKED/NOT_RUN. Required candidate cases and positive/negative backend cases have raw and structured evidence. |
| R4 Release candidate testing | Select an explicit version/commit and build the candidate once. Exercise exact wheel and sdist artifacts in clean environments outside checkout; validate supported optional compositions, real entry points, metadata/import origins, documentation and the release scientific corpus. Test the approval/evidence boundary and artifact tampering using non-publishing fixtures. | G0–G5 pass on the frozen candidate with an aggregate manifest and artifact hashes. Unauthorized/missing approval, different commit/bytes, failed or missing required checks produce NO-GO. Rebuilding or changing the candidate invalidates approval. No production tag/upload is needed to test this stage. |
| R5 Publishing and post-release verification | Only after explicit owner authorization: configure protected publishing/environment restrictions, enforce successful candidate gates and tag/version policy for every trigger, promote the approved bytes using existing OIDC publishing, then verify the version-specific registry and fresh installs. Document incident handling and safe deprecation/yank decisions with owner approval. | G6 records exact approval; rejected-approval and wrong-commit/hash/dispatch tests block. Published registry hashes equal approved hashes; G7 clean installs and smoke checks pass. Retain release evidence and release notes. No forced tag movement or silent replacement of released artifacts. |

## Priority and sequencing

1. **P0 — R1/R5 design:** close the publish bypass and artifact-identity gaps.
   Environment protection and repository-required checks must be verified by an
   authorized administrator; current publishing must not be treated as gated
   merely because a `pypi` environment is named.
2. **P0 — R3:** agree the scientific corpus and quantitative acceptance criteria;
   CI success alone cannot qualify scientific behavior. Resolve known safety
   xfails that conflict with release claims without weakening the tests.
3. **P1 — R1/R4:** make final-artifact installation and optional backend coverage
   explicit; avoid validating source-only environments in place of distributions.
4. **P1 — R2:** synchronize current behavior and version boundaries with examples
   and release notes, then enforce this with executable checks.
5. **P1 — R4/R5:** assemble the candidate record, protected approval, exact-byte
   promotion, registry verification and incident procedure.

Each stage is a separate scoped change with a reviewed diff and its own
executed exit evidence. Plan approval does not authorize starting the next
stage, tagging, dispatching publishing or contacting other people. Preserve
existing development/PR instructions and frozen scientific contracts throughout.

## Responsibility and decisions still needed

The project owner approves release scope, scientific thresholds, explicit
backend exclusions and the final publication. Implementers supply deterministic
checks/evidence; an authorized administrator controls protected environments and
repository rules. A named reviewer must verify the candidate evidence before
owner approval. Exact people, retention periods, backend-runner access, numeric
benchmark thresholds and any new Python/platform claims remain decisions for
their corresponding stage, not assumptions made by R0.
