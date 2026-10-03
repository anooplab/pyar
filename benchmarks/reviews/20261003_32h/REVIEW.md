# Review of the preceding 32 hours — 2026-10-03

## Scope

Reviewed recent production-code changes and the existing working tree at
`57fb99a6d086eba3a091c5be1a368c4964ed7f15`. The cumulative committed baseline was
`f48a6b662f052f5bdf1518d26ccc21cc19647a2c` (parent of the first commit in the window).
The scope includes the auto/HDBSCAN policy and basin-memory changes, clustering
benchmark and analysis tools, xTB model selection, climbing-image handoff,
TS benchmark instrumentation/adapters, and uncommitted ORCA comparisons and
backend-independent scan continuation. Large retained dataset/run files were
checked through source identities and report consistency rather than treated
as source-code edits. Existing work was preserved; no commit or push was made.

## Findings and fixes

| Area | Defect | Repair and regression evidence |
|---|---|---|
| Generic scans | Calculator state was written in the caller's working directory. | Execute within the scan directory, restore the caller directory even on exceptions; analytic constrained-optimizer test protects a caller sentinel. |
| Scan reuse | Completed scan-plus-relaxation runs trusted artifact existence alone. | Hash start, relaxed, scan trajectory, profile and candidate artifacts; rebuild when altered. Older results without hashes undergo validation/reconstruction. Generic scan loading also rejects different charge/spin settings. |
| Benchmark stage handoff | Custom frequency, endpoint and relaxation settings were dropped in later stages, causing dependency checks to use defaults. | One common validation-settings payload reaches every stage in both PyAR and ORCA runners; mocked end-to-end tests use nondefault settings. |
| Model provenance | Uppercase `XTB` passed one branch but acquired irrelevant DFT defaults in another. | Canonicalize software before model dispatch; reject unused method/basis consistently. |
| ORCA comparison | An explicitly selected xTB executable could differ from the one PyAR later used for validation. Validation settings and input hashes were incompletely recorded. | Reject executable mismatch before calculation; record all shared validation settings, model, coordinate system and copied input hashes. |
| Clustering benchmark cache | Changed inputs could reuse a stale `xtbopt.xyz` or restart file. A modified log was not hash checked; energy parsing selected the first TOTAL ENERGY occurrence. | Use a fresh retained attempt directory, hash output and log, take the final energy and reject nonfinite results. Tests exercise missing fresh output and log tampering. |
| Conformer analysis | Reference labels could be paired with modified source geometries or duplicate frame names. | Verify source ensemble SHA-256 against its run manifest and require an exact unique frame-label mapping. Existing CAMVES_I, FGG55 and WG01 source hashes all match. |
| ORCA report | Report generation inserted pilot-specific control results and conclusions without deriving them from supplied records. It checked starts but not recorded physical/validation compatibility. Exact sign-flip enumeration was unbounded. | Reject contradictory recorded provenance and duplicate/misfiled results; list missing evidence; remove hard-coded numerical claims and sample size; bound larger sign-flip analyses with seeded Monte Carlo. |
| Failed optimizer artifacts | Malformed final XYZ files could prevent a failed job from producing a benchmark result. | ORCA rejects invalid geometry before validation; both runners preserve analyzable failure records instead of crashing during final diagnostics. |

## Review and fix cycle

1. Inspected committed and uncommitted code paths, including scientific gates,
   cache identity, physical settings, stage dependencies and command construction.
2. Applied the first fixes and ran focused tests: **204 passed**. Initial full
   validation after additional integrity fixes: **864 passed**, **192 subtests**;
   strict Sphinx and package builds passed.
3. Reviewed the fixes again. Extended coverage for malformed optimizer geometry,
   ORCA execution/validation consistency, exception-safe scan directory restoration,
   changed spin, changed energy logs, and initial-geometry artifact integrity.
   Corrected the new ORCA mock to include the real convergence/termination banners.
4. Re-ran focused tests: **124 passed**. Final full validation is listed below.

The installed geomeTRIC **1.1.1** implementation was inspected: `climbSet`,
`climbers`, and `NEBParams.ncimg` correspond to the integration. Climbing-image
selection still supplies a candidate only; convergence, frequency and IRC/
endpoint confirmation remain separate gates. No optimizer, clustering,
deduplication or Hamiltonian default was changed by these review fixes.

## Final validation

Commands use the repository environment (`.venv/bin/python`):

- `python -m pytest -q`: **868 passed, 192 subtests passed**.
- `python -m sphinx -E -b html -W --keep-going docs docs/_build/html`: passed.
- `python -m build --no-isolation`: passed; wheel and source distribution built.
- `pyar-neb --help` and `pyar-cli scan-bond --help`: passed.
- `git diff --check`: passed.
- `python benchmarks/scan_bond/validate_backend_scans.py --output <fresh-directory> --backend xtb-gfn2`:
  passed with the real local xTB executable. Three finite H2 scan energies were
  -0.971086678722, -0.980467429012 and -0.981983694723 hartree.

Dependency deprecation warnings and the sandbox's unavailable Sella compilation
cache remain warnings; they did not fail tests. Packaging used installed build
dependencies rather than creating a new isolated environment. Sphinx regenerates
`docs/generated_api.rst`; its unrelated generated changes were restored after
inspection because that file was initially clean.

## Historical benchmark reanalysis

[Reanalysis](matched_internal_reanalysis.md) was regenerated from the retained
45-start matched-internal ORCA, geomeTRIC and Sella result collections. Success
counts remain **9/45, 10/45 and 11/45**, respectively. This was an analysis of
saved results, not a new optimizer benchmark run.

The old ORCA result schema lacks reference/endpoints input hashes and some
relaxation/endpoint validation settings. Those absent fields are reported as
unverified; existing records were not rewritten to invent provenance. The
new runner records these fields for future comparisons. The old numerical
outcomes alone do not establish a generally superior optimizer.

Reproduce the analysis with the updated script:

```sh
python benchmarks/ts_optimizers/analysis/analyze_orca_optts_comparison.py \
  <orca_matched_internal_run> <paired_geometric_sella_run> \
  --sella-run <sella_internal_run> --output <report.md>
```

## Limits

This review includes mocked external-program boundaries, real analytic
geomeTRIC/Sella tests and a real small GFN2-xTB scan. It does not constitute a
new real full reaction-cycle benchmark across ORCA, Gaussian and AIMNet2, or a
multi-Python CI run. The local xTB installation cannot run g-xTB; unsupported
requests fail explicitly. Generic scans and the reaction workflow change the
process working directory and should run concurrently in separate processes,
not threads. Previously unhashed cache records are conservatively invalidated.
