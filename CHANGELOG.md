# Changelog

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
