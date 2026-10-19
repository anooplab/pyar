# Scientific validation status for the release candidate

**Decision: INCOMPLETE — scientific release gate G2 is not passed.** This report
separates software execution evidence from scientific reference approval. It
records evidence on master `4c94f87217cf9a85af29077884fe2efefcf0fc20` (19 October
2026). A successful workflow or a finite output is not, by itself, validation
of a chemical prediction.

## Executed real-backend evidence

Scientific Tier 2 ran on the default branch at the exact candidate SHA above:
[workflow run 39192735168](https://github.com/anooplab/pyar/actions/runs/39192735168).
The job completed successfully and its artifact is named
`scientific-tier2-4c94f87217cf9a85af29077884fe2efefcf0fc20`. The artifact was
downloaded and inspected (32 files). Provenance recorded in it: Python 3.11.17,
pip 26.2.1, PyAR 1.5.1, Open Babel 3.2.1, xTB 6.7.1 (`edcfbbe`, conda build
2026-06-10), GFN2-xTB, one core. The request records neutral H atoms, singlet
total multiplicity, and reactant identity `InChI=1S/2H`, `[H].[H]`; the identity
error field is null.

The real calculation is a bounded three-point constrained H–H scan from 0.96 Å
to 0.74 Å. The reported energies are −0.971086678722, −0.980467429012 and
−0.981983694723 Ha. The final scan geometry has H–H distance 0.7400000005 Å;
all reported energies are finite, optimizer and physical-energy diagnostic gates
passed, and atom mapping and constraints passed. The reported maximum is at the
first point (0.96 Å), with a reported barrier of 0.0 kcal/mol relative to that
point. The structured result explicitly labels this a “scan maximum candidate;
not a confirmed transition state.” The product registry has no products. This
result validates this execution path only; it is not evidence for bond formation,
a reaction mechanism, a barrier, a TS, or an IRC connection.

Input SHA-256 values: H.xyz and both fragment XYZ files
`7dc85789455588c7e5be41c6d16fa7e63a01967ccb0258edbb54a756de022e96`.
Request SHA-256 `573e8b611a60b74f026301195cf7070c3d00a97f1a15506a66865c3ec60a7a19`;
scan-request SHA-256
`b813c4be394d0a28f797cdf03eac309c30785fd49df676e2efe3073188ceb46b`;
trajectory SHA-256
`bff41e30c6318b758d8fe85d14b532d6fd098adaadcced84aa504beb183f5614`;
final-scan XYZ SHA-256
`7185c791bac795e4a55a2672d58c09f62f367a0a83e91dd9372c1d7c6be151b2`.
The complete raw and structured files remain in the GitHub Actions artifact.
The workflow does not preserve backend stdout/stderr as separate raw files.

Additional local real-backend explorations are listed below. They are informative
integration evidence, not frozen reference benchmarks. Exact source SHA,
commands, inputs, hashes and raw outputs are recorded in the pre-existing
R3 report and local evidence directories named here; these local files are not
part of the repository or durable CI artifact, so rerun before relying on them
for a release decision:

| Case | Observed result | Status / limit |
|---|---|---|
| xTB H₂ analytic vs central-difference gradient | GFN2 energy −0.981983697485 Ha; component −0.021144456237862 Ha/Bohr vs finite difference −0.021144460000122 Ha/Bohr (step 1e−4 Bohr; difference 3.76e−9) | Exploratory local real calculation; no owner-approved tolerance. |
| xTB H₂ optimization | 0.96 Å start; converged at 0.7766217643 Å; energy −0.982686172097 Ha; gradient norm 5.46e−5 Ha/Bohr | Real backend path smoke; no independent equilibrium reference/tolerance approval. |
| xTB/RDKit butane conformer workflow | Seed 41; 3 requested, 2 generated, 1 selected; selected energy −13.665127571354 Ha | Exploratory; not an independent conformer ranking benchmark. |
| xTB water dimer microsolvation | One solvent/orientation; completed; selected energy −10.149004732546 Ha | Workflow smoke only; not solvation thermodynamics or a frozen candidate comparison. |
| AIMNet2 distorted methane gradient | Energy −40.5328364485656 Ha; analytic component 0.02677793882605 vs finite difference 0.02676932341608 Ha/Bohr (step 0.02 Bohr; difference 8.62e−6) | Real model evaluation; exploratory, no approved tolerance. Model SHA `a6d366b05c5997668f24a834584b4033d0cab67d13d2f8caf143b965018d3735`; Torch 2.14.0. |
| NH₃ inversion TS/frequency/IRC/NEB pilot | Sella 2.6.0, xTB GFN2. TS energy −4.416506797945 Ha; one imaginary frequency −971.19 cm⁻¹. Both IRC branches reported converged; PyAR did not confirm reactant/product connection. Three-image NEB converged in one cycle. | Exploratory only. Endpoints were optimizer-generated and no independent reference was used; IRC connection confirmation is false. Not a validated mechanism. |
| AIMNet2 methane optimization | Backend wrote a geometry, but PyAR classified it `INVALID_RESULT` / `BACKEND_UNSUPPORTED`; no result was accepted as converged. | Correctly not counted as a pass. This exposes a remaining real-backend optimization validation gap. |

## Initial corpus and tolerance decisions for the scientific owner

No entry in `benchmarks/scientific_validation/catalogue.json` is approved for
release. This table is a review package, not an automatic approval. No new
reference energies or thresholds are proposed from PyAR's own outputs.

| Candidate / objective | Independent reference and method | Proposed observable and tolerance | Rationale and status |
|---|---|---|---|
| Analytical energy-gradient consistency (small analytic potential) | Closed-form potential and derivative, independently evaluated by finite differences | Energy and gradient compared to analytic values; tolerance to be set from arithmetic precision and finite-difference step study before freezing | Best deterministic unit-level case; implementation exists. Owner must approve the exact fixture and numerical threshold. |
| H₂ GFN2 energy, gradient and optimization | Independent xTB release/version and input settings must be frozen; reference calculation and raw output retained separately | Energy/gradient/optimized distance thresholds remain unset pending independent calculation and repeatability study | Tier 2 is execution evidence, not an independent reference. **Pending owner approval.** |
| Structure and atom mapping invariants | Explicit atom correspondence, rigid transformations, connectivity and pair-distance invariants | Exact identity/mapping rules; numerical geometric thresholds must be fixed per fixture | Current structural corpus has constructed labels and ambiguous cases; eligible only as structural software validation, not energetic evidence. **Pending scope approval.** |
| Bias potential energy and gradient consistency | Closed-form AFIR/bias test cases with independent derivative evaluation | Compare component energies and gradients; set scale-aware tolerance after step-size/convergence analysis | No release threshold is frozen in the repository. **Pending independent fixture review.** |
| Backend parser and cache/convergence integrity | Deterministic captured backend outputs plus adversarial stale, malformed and nonconverged fixtures | Exact schema/state assertions; reject unverified convergence and incompatible provenance | Software safety tests, not scientific accuracy benchmarks. Code regressions are present as ordinary tests; review corpus inclusion. |
| Reaction-path benchmark (RGD1 or another licensed corpus) | RGD1 DOI 10.5281/zenodo.20489312 version 1.0; published dataset includes reference structures, but current provenance uses tblite 0.6.0 while xTB executable differs | No numeric thresholds approved; TS requires gradient/Hessian, reaction-coordinate mode, both IRC branches and endpoint basin identity | Dataset redistribution permission is not established and calculator mismatch matters. Obtain permission and independently reproduce references before use. **Blocked.** |
| Microsolvation candidate comparison | Existing MS0 protocol and independent method/settings | No numeric threshold approved; compare fixed candidate pool with common-method rescoring and clash/identity checks | Pilot is finite microsolvation, not bulk thermodynamics; external backend availability varies. **Pending.** |

The scientific owner must select the initial corpus, approve input/method/version
metadata, provenance, units, and per-case tolerances before these candidates can
be marked release eligible. Tolerances must be set before candidate comparison;
observed PyAR values above do not set them. Reference updates require old/new
hashes, rationale, independent evidence and explicit scientific review.

## Coverage and release blockers

- **Real xTB:** Tier 2 scan passed at the exact master SHA; local real energy/gradient,
  optimization, conformer and microsolvation explorations exist. None is an
  owner-approved reference corpus. Full reaction exploration, validated TS,
  frequency, IRC endpoint-connection and failure/recovery campaigns are not
  established by the Tier 2 scan.
- **AIMNet2:** Real energy/gradient exploratory calculation exists. Optimization
  result was rejected by existing safeguards; convergence semantics need a
  separately reviewed technical/scientific resolution before a release claim.
- **ORCA:** Executable was present locally but no real calculation was run.
  Gaussian, Psi4, MOPAC and Turbomole were not available for real validation.
- **Mocks and interface tests:** They verify software contracts only, not
  backend numerical accuracy or chemistry.
- **Release gate:** G2 remains **BLOCKED** until the owner approves the corpus and
  tolerances, required real-backend cases are executed with raw provenance, and
  outstanding convergence/path coverage is resolved or explicitly excluded
  from the supported release scope.
