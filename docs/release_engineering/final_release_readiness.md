# PyAR 1.5.2 release readiness and outcome

**Decision: released.** The owner-authorized PyAR 1.5.2 maintenance release is published on PyPI and GitHub. This is a narrow maintenance release and does not establish readiness for a major release.

## Release identity

- Version: `1.5.2`
- Package: `pyar-chem`
- Candidate/tag commit: `f621275c1ef19950bc562264f1319112eb3bc6ad`
- Tag: `v1.5.2`
- Post-release verifier fix commit: `7383f9140d82fef3b4b6d27667e3410498e7efca`
- PyPI: <https://pypi.org/project/pyar-chem/1.5.2/>
- GitHub release: <https://github.com/anooplab/pyar/releases/tag/v1.5.2>
- Exact candidate CI: <https://github.com/anooplab/pyar/actions/runs/39211560723>
- Exact candidate Scientific Tier 2: <https://github.com/anooplab/pyar/actions/runs/39211573722>
- Post-publication verification: <https://github.com/anooplab/pyar/actions/runs/39213822582>

## Release gates

| Gate | Status | Evidence |
|---|---|---|
| G0 — Approved maintenance scope | PASS | Owner approved the narrow 1.5.2 maintenance scope and explicitly excluded broader scientific claims. The exact release authorization is recorded in `releases/1.5.2.json`. |
| G1 — Code and regression checks | PASS | Exact candidate CI completed successfully. Python 3.10–3.13 each reported 1,948 tests, zero failures and zero errors. Fast regressions: 312 tests, zero failures/errors, one documented packaging-fixture skip. Scientific Tier 1: 79 tests, zero failures/errors/skips. |
| G2 — Scientific validation | PASS, scoped | Owner-approved H₂ GFN2 energy/gradient consistency evidence and the final candidate's real xTB GFN2 bond-scan smoke are recorded in `scientific_validation_final.md`. Broader scientific categories remain excluded and unvalidated. |
| G3 — Documentation | PASS | Strict Sphinx documentation job passed in exact candidate CI. |
| G4 — Installation and compatibility | PASS | Wheel/sdist clean installations passed in CI; optional-extra installation matrix passed on Python 3.10 and 3.13. Post-release Python 3.12 and local Python 3.14 clean PyPI installations passed, including `pip check`, import-origin and CLI checks. |
| G5 — Distribution | PASS | Exact candidate CI built and inspected the wheel and sdist. Published files match the approved SHA-256 hashes in the approval record. |
| G6 — Owner authorization and protected promotion | PASS | The versioned approval record was merged to `master`; the PyPI environment permits deployments only from `master` and administrator bypass is disabled. The guarded workflow verified run identity, job inventory and artifact hashes before upload. |
| G7 — Registry and installed-package verification | PASS | The verification-only workflow confirmed version-specific PyPI hashes and a clean installed package. The initial combined publisher's post-job was canceled after an API-response parsing defect was found; the defect was fixed, covered by CI and the verification-only workflow passed. |

## Scientific limits

This release's scoped evidence does not approve a broad benchmark corpus or universal numerical thresholds. It does not validate reaction exploration, NEB, transition-state identity, frequencies, IRC, aggregation/conformer ranking, microsolvation, or scientific compatibility with backends other than the recorded xTB cases. Do not cite 1.5.2 as major-release scientific qualification.

The next major release remains blocked on the broader scientific corpus, scientific-owner review of references and tolerances, and real-backend workflow coverage specified by the release gates.
