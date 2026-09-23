# Changelog

## Unreleased

### Added

- Structured request models for aggregation, conformer search, and solvation.
- gXTB energy-and-gradient provider support.
- A conformer benchmark and failure-diagnosis workflow that separates
  generation, selection, backend-refinement, ranking, and input-chemistry
  failures.
- A soft-minimum interfragment reaction bias, selectable through
  `--bias-potential softmin`.
- Fixed, scheduled, and adaptive reaction-bias controllers in `pyar-cli react`
  and `pyar-react`. Adaptive runs use the accepted-step geomeTRIC driver and
  expose alpha bounds, safety margin, and smoothing; fixed bias remains the
  default.

### Changed

- Conformer generation now uses a broader candidate pool and removes duplicate
  structures more aggressively before backend refinement. Basin selection
  preserves low-energy, geometrically diverse folded, open, and outlier
  candidates.
- Aggregation, conformer, and solvation workflows now use dedicated request
  models to validate and preserve workflow inputs.
- Reaction searches now use the generic `--bias-min` and `--bias-max` strength
  options. The historical `-gmin`/`--gmin` and `-gmax`/`--gmax` options remain
  accepted as compatibility aliases.
- geomeTRIC reaction optimization now selects either the AFIR or soft-min bias
  through `--bias-potential`; AFIR remains the default.

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
- Conformer backend-pool selection now protects low-energy, RMSD-diverse,
  contact-rich folded, open, and outlier basins instead of applying a global
  compactness bias.
- Conformer torsion search now uses a single stratified random torsion-kick
  implementation by default; experimental Bonobo, adaptive, basin-hopping,
  grid, Monte Carlo, and evolution probes were removed from the public workflow.
- Conformer generation defaults were widened to use more RDKit conformers and
  seeds, and the post-generation pool is now collapsed more aggressively before
  backend refinement so similarity duplicates do not crowd out distinct basins.

### Fixed

- Python 3.10 now has a TOML fallback for runtime version detection.
- Reaction workflow summaries now prefer coordinate paths that actually exist.
- Conformer benchmark diagnosis now reports RDKit and backend energy fields with
  explicit units and keeps RDKit energy-window classification in native
  force-field units.
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

## Historical development - 2017-12-01 to 2024-10-06

The repository predates its consistently versioned releases. The summaries
below group the major, dated changes from the Git history without assigning
retrospective version numbers.

### 2024 - 2024-06-11 to 2024-10-06

#### Added

- AIMNet2, AIQM1, xTB-AIMNet2, and xTB-AIQM1 integration and associated model
  packaging support.
- New clustering descriptors and an AFIR test path.
- NetworkX-based tabu functionality and a wall-potential update.

#### Changed

- Updated command-line scripts, clustering, exploration, reaction, and
  optimization workflows.
- Reworked OpenBabel and MOPAC integration, runtime requirements, and model
  executable paths.

### 2022-2023 - 2022-02-28 to 2023-11-10

#### Changed

- Corrected scan and tabu atom-index handling, distance scaling, and xTB
  optimization controls.
- Replaced the slow `mendeleev` dependency with local atomic data and updated
  molecular representations.

#### Fixed

- Resolved aggregation traversal and representation-import errors, and restored
  compatibility with Python versions later than 3.10.

### 2020-2021 - 2020-04-02 to 2021-09-24

#### Added

- The `pyar-cli`, `pyar-tabu`, `pyar-optimiser`, and `pyar-clustering` command
  entry points.
- Generalized aggregation for arbitrary numbers and types of components,
  together with per-molecule charge, multiplicity, and SCF-type settings.
- Bond-scan support, molecular representations, XYZ crawling, clustering
  descriptors, and reaction restart/checkpoint support.

#### Changed

- Consolidated binary and ternary aggregation into a common aggregation
  workflow and moved the package into the `pyar/` source directory.
- Reworked AFIR restraint gradients with Autograd and improved tabu placement,
  random-pathway selection, and aggregation permutation performance.

#### Fixed

- Improved handling of failed tight optimizations, SCF-convergence retries,
  single-geometry clustering, and user-name lookup on shared systems.

### 2017-2019 - 2017-12-01 to 2019-10-03

#### Added

- The initial PyAR aggregation and reaction codebase, XYZ I/O, molecule model,
  and package installation support.
- Interfaces for ORCA, Turbomole, xTB, Psi4, and Gaussian.
- AFIR restraints with covalent-radius data, site-specific orientation support,
  bond restraints, and optimization convergence controls.
- K-means/MeanShift clustering, binary and ternary aggregation, configurable
  maximum seed counts, and hydrogen-bond analysis.

#### Changed

- Improved proximity checking, Turbomole and xTB optimization behavior,
  reaction/aggregation command-line handling, and backend restart controls.
