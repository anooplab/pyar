# PyAR 1.5.2 scientific validation record

**Status: PASS for the explicitly approved maintenance scope; broader scientific validation remains incomplete.** This report distinguishes backend execution evidence from evidence that a scientific method or workflow is generally correct.

## Exact candidate and evidence

- Release candidate: `f621275c1ef19950bc562264f1319112eb3bc6ad`
- Candidate CI: <https://github.com/anooplab/pyar/actions/runs/39211560723>
- Scientific Tier 2: <https://github.com/anooplab/pyar/actions/runs/39211573722>
- Raw Tier 2 artifact: `scientific-tier2-f621275c1ef19950bc562264f1319112eb3bc6ad`; GitHub artifact ID `12167911527`; archive SHA-256 `f1392bd3320cac4708fdd54f77daba70b5ae68afcf58ad5dd1a32adf348246ea`.
- The raw artifact contains 32 files. A SHA-256 manifest is retained with the release evidence. The candidate SHA in the artifact is the release candidate SHA above.

## Executed cases

| Case | Purpose and method | Observed evidence | Status and interpretation |
|---|---|---|---|
| H₂ energy/gradient consistency | Owner-approved neutral-singlet H₂, GFN2-xTB; compare backend gradient with independent finite differences of energy. | Executed against the immediately preceding candidate `3c3e92af07d702d4b3b4e44e39b238fed99d0b61`; xTB `6.7.1-6.fc44.fc44` (Fedora build, 2026-01-18). Maximum absolute component difference `1.113897222342064e-07 Hartree/Bohr` versus approved `1e-6`; relative L2 difference `4.854621070550338e-06` versus approved `1e-4`. Input manifest SHA-256 `b07cca1bba079445b72b9f317b80d1c19e4952419810433ae4db56ed97c0a33f`; raw-artifact hashes verified (109 files). | PASS for this one approved derivative-consistency case. This does not establish general energy or gradient accuracy. The final candidate adds release version/changelog and release-promotion integration; it does not change this scientific calculation implementation. |
| H₂ three-point bond scan | Real backend smoke of the final candidate using xTB GFN2; neutral singlet, one process, constrained scan from 0.96 to 0.74 Å. | xTB `6.7.1` (`edcfbbe`, conda build, 2026-06-10). Input fragments each SHA-256 `7dc85789455588c7e5be41c6d16fa7e63a01967ccb0258edbb54a756de022e96`. Energies at 0.96, 0.85 and 0.74 Å: `-0.971086678722`, `-0.980467429012`, `-0.981983694723` Hartree; all finite. Run and raw output report success. | PASS as a backend and workflow-execution smoke only. It is not a reference-energy benchmark, reaction-path validation or confirmed transition state. The reported scan maximum is explicitly a candidate only. |

## Coverage and exclusions

The release suite also executed deterministic scientific Tier 1 tests (79 passed, no skips) and the full Python test suite (1,948 per Python version; 13 tests skipped per version because optional `irmsd` or `geodesic-interpolate` dependencies were unavailable). The fast-regression report had one skip because a built wheel was not available in that isolated test step; wheel/sdist installation and archive checks ran in their dedicated CI job.

The following were not established by this release: independently sourced molecular-energy references beyond the approved derivative check, general optimization accuracy, reaction exploration correctness, NEB path continuity, a converged first-order saddle and Hessian signature, IRC endpoint connection, conformer ranking, aggregation, microsolvation, or real numerical validation of ORCA, AIMNet2, MLatom, Sella, and other non-xTB integrations. These are explicit maintenance-release exclusions and remain requirements for the broader R3/major-release validation effort.

Mocked adapter tests and successful process exits are not counted as scientific validation. No claim of R3 scientific completion is made.
