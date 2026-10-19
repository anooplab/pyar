Scientific validation
=====================

This directory defines the evidence PyAR requires before making scientific
claims. The validation framework separates deterministic tests, interface
checks, reference comparisons, and full workflows. A passing mocked boundary
test establishes software behavior only; it is not a real-backend chemistry
result.

Architecture and tiers
----------------------

Tier 1 runs on every pull request. It contains deterministic analytical,
finite-difference, and structural invariants with no external quantum chemistry
program. CI selects it with ``python -m pytest -m "scientific and tier1"``.
The result is retained as ``scientific-tier1.xml`` and is required through the
existing ``R1 required checks`` dependency on the quality job.

Tier 2 is an inexpensive real-backend smoke run on a scheduled or release
runner. Record actual executable and version, input and reference hashes,
charge, multiplicity, method, units, raw logs, energy/gradient/optimization
evidence, and the candidate SHA. It must report ``NOT_TESTED`` when the backend
or reference is unavailable. The local GFN2-xTB scan records an actual backend
smoke; it is not a frozen release corpus.

Tier 3 is extended, explicitly budgeted validation for major releases or
scientific algorithm changes. It may require licensed datasets, proprietary
programs, model weights, dedicated runners, repeated trials, or frequency and
IRC analyses. The test plan must identify the responsible runner and store raw
and structured evidence. Missing prerequisites are ``BLOCKED`` or
``NOT_TESTED`` and block a required release gate.

Tests carry pytest markers ``unit``, ``integration``, ``regression``,
``scientific``, ``external_backend``, ``slow``, and ``tier1``/``tier2``/``tier3``.
Markers describe evidence type; the full suite still runs all configured tests.
``tests/conftest.py`` conservatively labels unknown legacy tests as integration.
Mocks and recorded fixtures must not be counted as external backend executions.
The current pytest suite has no test that performs a real quantum-chemistry
calculation; the xTB Tier 2 runner is a separate workflow. Long benchmark
campaigns are standalone runners, so the ``slow`` marker is available but no
current pytest case uses it.

Reference data and governance
-----------------------------

The current inventory is
``benchmarks/scientific_validation/catalogue.json``. It points to existing
reference material and states what each source can establish. None of its
entries is currently approved as the complete release corpus. In particular,
the structure-comparison dataset validates construction-based structural
labels, not energy surfaces; the RGD1 dataset is external and its record does
not grant redistribution permission; and historical scan results are smoke
evidence, not a candidate-bound release result.

Keep small permitted inputs and immutable reference metadata under
``benchmarks/scientific_validation/reference_data/``. Every case must specify
system, ordered input geometry and SHA-256, charge, multiplicity, method,
backend, relevant settings, reference provenance/version/hash, metrics,
units, acceptance tolerance and validation procedure. Do not overwrite a
reference. A change requires an explicit scientific rationale, source and
method/backend version, old/new hashes, and scientific-owner review.

``tools/scientific_validation.py`` validates an evidence file against a frozen
case manifest. It checks reference hashes, electronic state, method/backend,
input hash, backend version, numerical metrics, tolerances and raw log links.
It verifies that each log exists inside the evidence directory and records its
SHA-256. The expected candidate SHA is required and must match the results
record. It emits JSON and Markdown with ``PASS``, ``FAIL``, ``NOT_TESTED``,
``NOT_APPLICABLE`` or ``BLOCKED`` case states. The tool does not execute
chemistry, convert a missing result into a pass, or approve exceptions. Example
invocation:

.. code-block:: console

   python tools/scientific_validation.py cases.json results.json \
     --candidate-sha "$(git rev-parse HEAD)" \
     --json evidence/scientific.json --markdown evidence/scientific.md

Tolerance policy
-----------------

There is no universal cross-method chemistry tolerance. Each frozen case must
provide a tolerance in the metric's stated unit. Use these classes:

* **Exact invariants:** atom symbols/order, charge, multiplicity, number of
  frames, connectivity labels, selected backend/method, and required workflow
  stages. These are exact comparisons.
* **Analytical/numerical derivatives:** compare analytic gradients or forces
  with independent finite differences or autodiff. Record the displacement,
  norm/component metric, absolute and relative tolerance, and scale. Compare
  values close to zero with an absolute floor; do not use relative error alone.
* **Same-method numerical replay:** energy, gradient, and optimized-geometry
  tolerances are set per backend/version, precision, thread policy, and method.
  Compare energies in Hartree, gradients in Hartree/Bohr (or explicitly
  converted units), and coordinates in Angstrom. Save both absolute and
  relative energy where meaningful.
* **Structural equivalence:** declare atom mapping, alignment policy, fragment
  policy, RMSD in Angstrom, and identity/connectivity checks. A distance cutoff
  by itself does not establish chemical equivalence.
* **Reaction-path evidence:** declare barrier metric and unit, scan constraint
  error, path continuity, TS gradient norm, frequency method and cm\ :sup:`-1`
  threshold, imaginary-mode displacement, and endpoint identity. An imaginary
  frequency alone does not qualify a TS; IRC branches must connect to the
  intended minima.
* **Conformer/aggregate/microsolvation search:** declare seed, search budget,
  duplicate rule and energy-ranking comparison. Assess acceptable alternative
  basins separately from hard regressions; do not require identical stochastic
  discovery order across platforms.

The tolerances in Tier 1 tests are local mathematical test tolerances, not
release-wide backend thresholds. Release thresholds remain unset until the
scientific owner approves a candidate corpus and numerical conditions.

Backend validation matrix
-------------------------

``pyar/backend_capabilities.py`` and ``pyar/dependencies.py`` define current
routing and prerequisites. This matrix distinguishes implemented code paths
from demonstrated numerical chemistry. “Mock/fixture” never counts as a real
backend calculation.

.. list-table:: Current integration evidence
   :header-rows: 1

   * - Interface
     - Current evidence
     - Release validation gap
   * - xTB GFN1/GFN2
     - Adapter, parsing, energy/gradient and failure tests; recorded GFN2 scan smoke
     - Repeat against frozen candidate/input and record convergence and hashes
   * - xTB g-xTB
     - Capability gate rejects executables without ``--gxtb``
     - Real g-xTB runner and model version are unavailable in this environment
   * - ORCA
     - Input/provider/parser fixtures; recorded ORCA GFN2-xTB scan smoke
     - Candidate-bound licensed executable tests and broader method coverage
   * - Gaussian
     - Adapter/provider and subprocess behavior tested at mocked boundaries
     - Proprietary calculation not run
   * - Turbomole and xTB-Turbomole
     - Adapter and failure handling in unit/integration tests
     - Real installation matrix and versioned numerical runs not recorded
   * - Psi4, MOPAC, OpenBabel
     - Available adapter paths and dependency/error behavior tested
     - Release-specific real calculation case coverage remains incomplete
   * - AIMNet2
     - Python/model discovery, adapter contracts, and missing-asset diagnostics tested
     - Model asset absent from this checkout; numerical model calculation not run
   * - MLatom AIQM1, ANI, composite models
     - Optional dependency and routing tests where present
     - Model/licensing/backend-specific calculation validation is not complete
   * - ASE, geomeTRIC, Sella
     - Optimizer integration tests; scan validation includes analytic potential
     - Required real-backend combinations need a release support decision

Workflow coverage and known gaps
--------------------------------

.. list-table:: Implemented workflow validation
   :header-rows: 1

   * - Workflow
     - Existing coverage
     - R3 status/remaining gap
   * - Aggregate and molecular growth
     - Geometry, selection, state/restart, and mocked optimizer workflow tests
     - Interface covered; no frozen independent energy-ranking release case
   * - Conformer generation/ranking
     - Reproducibility, comparison harness, and external ensemble pilots
     - Pilot only; no approved release corpus or production reference-energy set
   * - Reaction exploration and adaptive AFIR
     - Bias finite differences, topology/state tests, failure/recovery tests
     - Tier 1 math covered; actual pathway corpus not frozen
   * - Bond scanning
     - Analytic constraint tests and recorded real xTB/ORCA smoke evidence
     - No release-candidate real scan campaign
   * - NEB, TS, IRC
     - Analytic/mocked stages, endpoint and path validation; RGD1 adapter
     - External dataset execution and full real TS/frequency/IRC evidence not run
   * - Microsolvation
     - Candidate policy tests and MS0 comparative runner/protocol
     - External comparators not guaranteed; no approved release baseline
   * - Optimization failure/recovery
     - Extensive negative/replay tests; 17 strict expected failures recorded
     - Those expected failures are open defects and release blockers for affected claims

Failure reporting and release gate
----------------------------------

Classify each case as scientific regression, numerical tolerance failure,
backend incompatibility, infrastructure failure, missing backend, unsupported
workflow, nondeterministic result, or reference inconsistency. Retain case ID,
system, method, expected/observed values, tolerance/units, backend/version,
candidate SHA, input/reference hashes, raw logs, and supported cause. Do not
label unavailable calculations PASS. Tier 1 is required on every PR. For each
release, explicitly select required Tier 2/Tier 3 cases and authorized
exclusions; any required unexecuted case blocks approval. A scientific
explanation cannot override a failed comparator.

To add a benchmark, first identify an independent permitted source; record its
license and provenance; freeze exact inputs and their hashes; define method,
settings, electronic state, expected metrics and tolerances; write the case
manifest without overwriting prior references; implement a deterministic
runner/comparator with positive and negative controls; test missing backend and
stale-hash behavior; then request scientific-owner review before assigning
release-required status.

.. toctree::
   :maxdepth: 1

   benchmark_catalogue
   numerical_tolerances
   reference_provenance
   backend_validation
   troubleshooting
