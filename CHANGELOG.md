# Changelog

## 1.2.0 - 2026-06-12

This release adds a new RDKit-based conformer workflow and tightens release compatibility across supported Python versions.

### Added

- A new conformer workflow that builds, ranks, and optionally refines conformers from SMILES, SDF, or MOL input.
- A dedicated `pyar-conformer` console script for the new workflow.
- A `conformer` optional dependency group for RDKit-based conformer generation.
- A conformer benchmark and failure-diagnosis framework to classify whether missed minima are due to generation, selection, backend refinement, ranking, or input chemistry.

### Changed

- The package version was bumped to `1.2.0`.
- Release publishing now follows the tagged `v1.2.0` workflow path.

### Fixed

- Python 3.10 now has a TOML fallback for runtime version detection.
- Reaction workflow summaries now prefer coordinate paths that actually exist.
- CI test isolation was improved so OpenBabel and AIMNet2-dependent tests do not rely on local machine state.

## 1.1.1 - 2026-05-31

This release focuses on cleanup, documentation, and workflow stabilization.

### Added

- A documented module layout with canonical imports and compatibility shims.
- Dedicated documentation pages for aggregate, react, solvate, backends, external programs, migration, installation, quickstart, and API reference material.
- A more explicit GitHub Actions CI workflow with test, docs, and build jobs.
- A project changelog for release tracking.

### Changed

- Canonical implementations now live in the planned modules under `pyar.core`, `pyar.io`, `pyar.sampling`, `pyar.selection`, `pyar.workflows`, `pyar.backends`, and `pyar.biases`.
- Legacy import paths remain as thin compatibility shims where needed.
- Workflow state objects now record structured metadata for restart and provenance.
- Sampling, selection, backend capability, and workflow logic were reorganized without changing the intended scientific behavior.
- Public docs were shortened and reorganized around user tasks and installation guidance.

### Fixed

- Version metadata is now consistent across packaging and runtime checks.
- CI now exercises the package across multiple Python versions and validates docs and source distributions.
- Test isolation issues that could create stray files in the repository root were removed.
- Connectivity filtering is now policy-controlled. It is off by default for molecular aggregates and solvation, preferred for atomic/formula growth, and user-overridable through `--connectivity-policy`.

## 1.1.0 - 2026-05-30

The 1.1.0 release concentrated on packaging modernization and dependency cleanup.

### Changed

- Packaging moved to `pyproject.toml` with setuptools as the build backend.
- Optional dependency groups were split out for test, docs, selection, ML, xTB, AIMNet2, and OpenBabel-related features.
- Heavy backend dependencies were externalized where possible.
- The package kept the `pyar` import namespace while distributing as `pyar-chem`.
