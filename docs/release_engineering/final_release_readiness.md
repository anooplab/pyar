# Final release readiness

**Decision: NO-GO — PyAR is not release ready.** State recorded at master
`4c94f87217cf9a85af29077884fe2efefcf0fc20` (19 October 2026). Package
metadata is 1.5.1. No release tag, GitHub release or package publication was
created. PyAR 2.0 must not be published from this evidence.

## Gate status

| Gate | Status | Evidence / remaining requirement |
|---|---|---|
| G0 Candidate | BLOCKED | Exact SHA recorded; no 2.0 version candidate or frozen release manifest/changelog approval. |
| G1 Code and regressions | PASS (exact-head CI) | PR #78 head `a30cf238cea504c384dfe6dcc178cc5aeb02dae7` passed required CI before merge to this master. Python 3.10–3.13, lint/fast regressions, optional installs, docs, clean wheel/sdist installs and R1 aggregate checks passed. The former 17 optimization safety cases are ordinary regression tests. |
| G2 Scientific/backend | BLOCKED | Tier 2 real xTB scan passed on exact master SHA; it is only a constrained scan smoke. No approved benchmark corpus/tolerances; real TS/IRC endpoint validation incomplete; AIMNet2 optimization is rejected by result validation. See [scientific report](scientific_validation_final.md). |
| G3 Documentation | PASS | Exact-head documentation CI and a fresh strict local Sphinx build both passed after correcting the release-engineering toctree links. |
| G4 Installation/optional compatibility | PASS (CI evidence) | Exact-head clean wheel/sdist installation and core, aggregate, conformer, reaction, AIMNet2, NEB, ML, analysis and microsolvation combinations passed on Python 3.10 and 3.13. |
| G5 Distribution | BLOCKED for final artifacts | Local wheel and sdist build plus `twine check` passed for metadata version 1.5.1; these were exploratory, uncommitted-tree artifacts and are not frozen release artifacts. Exact-head CI also passed clean artifact installation. |
| G6 Approval/publish safety | NOT RUN | No approval of version, SHA, artifacts or publication; no tag/publish action. |
| G7 Post-release | NOT APPLICABLE | Nothing was published. |

## Integration and scientific Tier 2 evidence

PR #75 merged as `255ce0e50242d94e233152547b848dc89a18c4aa`; PR #76 as
`2b2c9c99dfefb8d27f8ff649ec1570857fd87322`; PR #78 cache backend-version
provenance as `4c94f87217cf9a85af29077884fe2efefcf0fc20`. PR #78 exact-head
checks: [PR run 39192224563](https://github.com/anooplab/pyar/actions/runs/39192224563)
and [push run 39192219869](https://github.com/anooplab/pyar/actions/runs/39192219869).

Scientific Tier 2 passed on exact master SHA:
[run 39192735168](https://github.com/anooplab/pyar/actions/runs/39192735168).
The 32-file artifact is `scientific-tier2-4c94f87217cf9a85af29077884fe2efefcf0fc20`.
Recorded versions: xTB 6.7.1 (`edcfbbe`), Open Babel 3.2.1, Python 3.11.17,
PyAR 1.5.1. Three finite H₂ scan energies were −0.971086678722,
−0.980467429012 and −0.981983694723 Ha; final constrained distance was
0.7400000005 Å. The first point is the scan maximum, the result explicitly says
it is not a confirmed TS, and no product was registered. This is not reaction,
barrier or TS validation.

## Remaining blockers and recommendation

1. Scientific owner approval of initial reference corpus, provenance and
   per-case thresholds (decision package is in the scientific report).
2. Real-backend validation for supported workflow scope, especially optimization
   convergence, reaction exploration, TS Hessian/mode and both IRC endpoint
   connections. Resolve AIMNet2 optimization rejection or explicitly exclude it.
3. Freeze a version-specific release candidate, changelog/compatibility notes,
   final artifact hashes and clean installation evidence for those exact bytes.
4. A fresh strict local Sphinx build, wheel/sdist build, and `twine check` passed. A local clean install did not complete because dependency downloads were too slow in the session; exact-head CI clean installation checks passed.

R3 scientific completion and publication readiness are not claimed. Do not tag or
publish PyAR 2.0. The next safe action is scientific-owner review of the corpus,
tolerances and supported backend/workflow scope, followed by the missing bounded
calculations on a frozen candidate.
