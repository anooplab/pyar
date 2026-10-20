PyAR scientific validation plan
===============================

Purpose and boundary
--------------------

This document is the standalone plan for validating PyAR's scientific
behavior. It defines what evidence to collect, how to judge that evidence, and
what is required before scientific claims are treated as release-qualified.
Software test success, backend interface compatibility, and scientific
validation are separate outcomes. A mocked backend test verifies software
behavior only; it never counts as a real calculation.

The plan covers implemented functionality only. A missing or incomplete
workflow is reported as ``NOT_APPLICABLE`` or ``NOT_TESTED`` with its reason;
it is never silently counted as a pass. No reference result or tolerance in
this plan is approved merely because it appears in an existing test, pilot, or
historical report. Scientific-owner approval is required for the frozen
release corpus, thresholds, backend exclusions, and exceptions.

Validation principles
---------------------

* Preserve the approved method, units, charge, multiplicity, settings, and
  workflow objective. Do not modify scientific algorithms to make a test
  pass.
* Bind every result to the exact source commit, input bytes, reference bytes,
  backend executable/model version, and relevant numerical settings.
* Retain raw outputs. A successful process exit does not prove convergence,
  valid geometry, or a correct chemical result.
* Distinguish analytical and structural invariants, mocked/interface tests,
  real backend calculations, and end-to-end workflow evidence.
* Use explicit budgets and bounded time for expensive jobs. Missing tools,
  licenses, model assets, or infrastructure are reported as unavailable, not
  as successful validation.
* Never overwrite reference data in place. A reference change needs source,
  rationale, old/new hashes, method and version details, and scientific-owner
  review.

Evidence tiers
--------------

Tier 1: deterministic scientific tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run on every pull request using ``python -m pytest -m "scientific and tier1"``.
These tests must be fast, deterministic, and independent of proprietary
software. Use analytical results, finite differences, conservation or
symmetry invariants, and independent implementations where suitable. The
selection must include, where implemented:

* coordinate shape, units, atom identity/order, mappings, and rigid transforms;
* molecular connectivity, fragments, bond changes, and structural comparison;
* energy/gradient and force-sign consistency, including bias components;
* finite-difference checks for potentials, collective variables, constraints,
  and restraints;
* optimizer and cache state invariants, including explicit convergence and
  rejection of malformed or stale results;
* reaction topology and endpoint consistency.

Tier 1 establishes deterministic software/science invariants. It does not
establish backend accuracy or validate a full chemical pathway.

Tier 2: bounded real-backend smoke tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run on a scheduled or manually dispatched runner with a freely available,
supported backend where possible. The current intended starting point is xTB.
Each run records the candidate SHA, PyAR and backend versions, exact inputs and
SHA-256 hashes, method, charge, multiplicity, numerical settings, raw output,
structured result, and validation verdict. A result without the raw artifacts
or exact candidate binding is incomplete.

The initial bounded set is:

1. Energy and analytical-gradient comparison against central finite
   differences, with step size, units, component/norm metric, and scale-aware
   tolerances recorded.
2. Geometry optimization from a perturbed structure, requiring explicit
   backend convergence evidence, finite results, preserved atom order, valid
   geometry, and a final gradient check.
3. One small aggregate/conformer or reaction-exploration workflow that exists
   in the repository, with fixed seed/budget where applicable and candidate
   structures independently checked after unbiased relaxation.
4. A small microsolvation proposal and relaxation, when the implementation and
   backend permit it, checking solvent integrity, placement, clashes, and
   finite energy/gradient.
5. Controlled failure/recovery tests that verify the input is preserved,
   failures are recorded, and failed candidates cannot be reported as
   successful. These are software reliability tests, not chemistry evidence.

Every unavailable case is reported ``NOT_TESTED`` with the exact missing
prerequisite. A Tier 2 workflow's green status alone does not establish
scientific correctness.

Tier 3: extended pathway and backend validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run before major scientific releases or changes to scientific algorithms.
Declare runner, backend, data/license access, compute and wall-time budgets,
repeats, and artifact retention before dispatch. Include representative
reaction pathways and selected backend comparisons only when the workflow
exists and the case has an independent scientific basis.

For NEB, transition-state, frequency, and IRC work, require as applicable:

* independently specified reactant and product minima and a stable atom map;
* endpoint identity and connectivity checks, continuous path geometry, and
  documented path settings;
* converged TS optimization with finite energy, forces, and coordinates;
* Hessian evidence for one significant imaginary frequency for an ordinary
  first-order saddle, plus a mode displacement consistent with the proposed
  reaction coordinate;
* both IRC branches followed by optimization, with each branch reaching the
  intended endpoint basin.

One imaginary frequency is insufficient to prove a reaction connection.
Report alternate minima, failed trials, and uncertain assignments rather than
selecting only favorable trajectories.

Scientific coverage matrix
--------------------------

The following matrix guides case selection. It describes coverage goals, not a
claim that all cases have been executed or approved.

.. list-table:: Planned coverage by implemented workflow
   :header-rows: 1

   * - Area
     - Minimum validation evidence
     - Evidence class
   * - Geometry and molecular topology
     - Atom identity/order, shape, mapping, connectivity, fragment count,
       bond changes, overlap checks, and structural equivalence with declared
       alignment policy
     - Tier 1 plus selected real examples
   * - Energies, gradients, and forces
     - Unit/sign checks, total/component consistency, finite differences or
       independent derivative reference, scale-aware tolerance
     - Tier 1 and Tier 2
   * - Optimization and recovery
     - Explicit convergence, finite final state, final gradient, request/cache
       identity, stale-output rejection, restart integrity, geometry restore
       on failure
     - Tier 1, backend interface tests, and Tier 2
   * - Aggregate growth and conformers
     - Valid structures, deterministic fixed-seed behavior where meaningful,
       duplicate handling, energy ordering and distinctness with independent
       reference evidence
     - Tier 1 and bounded Tier 2; extended search at Tier 3
   * - Reaction exploration and adaptive bias
     - Bias energy/gradient, physical and biased energy separation, topology
       changes, candidate validity, unbiased relaxation and failure records
     - Tier 1; representative real workflow at Tier 2 or Tier 3
   * - Bond scanning
     - Correct distance schedule and constraints, continuation behavior,
       finite energies, atom identity, and independently checked endpoints
     - Tier 1 plus real-backend smoke; pathway claims require Tier 3
   * - NEB, TS, frequency, and IRC
     - Endpoint/path continuity, converged saddle, Hessian/mode, IRC branches,
       and endpoint-basin connection
     - Tier 3 real backend
   * - Microsolvation
     - Site-directed placement, intact solvent geometry, clash handling,
       finite interaction/total energies, relaxation, and candidate provenance
     - Tier 1 and bounded Tier 2; bulk thermodynamics excluded unless separately
       validated
   * - Backend failures
     - Error classification, isolation, preservation of input state, no false
       success, restart/cache safety, and durable diagnostic artifacts
     - Tier 1/interface tests; real executable failure probes where practical

Optimization and cache safety case set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimization-recovery audit identified concrete negative-path assertions.
Keep each as an ordinary regression test once its root cause is repaired; do
not convert a reproduced defect to a skip or weaken its assertion. Before
using this inventory for a candidate, inspect current tests and results because
it records a historical audit baseline.

* ORCA must reject normal process termination without optimizer convergence
  evidence; missing or nonfinite final energy; stale output with mismatched
  symbols; and malformed singleton coordinates. Missing energy must never be
  replaced with zero.
* xTB must not let an old success marker override current-run SCC or execution
  failure, regardless of marker-file arrangement.
* Generic cache reuse must require request identity, method/backend and
  settings compatibility, atom count/symbol/order, valid geometry, finite
  results, and explicit completion/convergence evidence.
* geomeTRIC must not promote cycle-exceeded/provisional state to success. It
  must reject nonfinite energy or forces, frame/atom mismatch, missing explicit
  convergence evidence, and unrelated stale fallback XYZ output.
* A backend failure must preserve the original in-memory molecular geometry.
  Truncated or malformed cached output must become a typed backend error, not
  terminate the Python process.

The historical audit enumerated 17 parameterized cases across these behaviors.
For candidate qualification, map each current test identifier to its fix or
verified disposition, run it as a normal test, and verify both failure and
valid-success paths. Cache provenance must include backend executable or
model version where available; unknown versions must not silently authorize
reuse across incompatible requests.

Initial evidence inventory
--------------------------

Existing catalogue entries and pilot records are starting points for review;
they are not approved release references. The following historical inventory
must be rechecked against immutable source files before cases are frozen.

.. list-table:: Existing candidate material and its current evidentiary limit
   :header-rows: 1

   * - Candidate material
     - Source and expected observable
     - Current limit and disposition
   * - Structural comparison corpus
     - Constructed geometries and pair labels for structural equivalence and
       non-equivalence
     - Supports narrow structural invariants only; construction labels and
       exploratory RMSD sweeps are not calibrated energy-surface references.
       Not approved.
   * - MPCONF196GEN clustering pilot
     - Source conformer ensembles for selected macrocycles; operational basin
       and ranking comparisons under GFN2-xTB
     - Source license is recorded, but upstream energy descriptions and
       numeric comments require reconciliation. Small system count and
       operational basin labels limit generalization. Pilot only.
   * - RGD1 transition-state material
     - External dataset (Zenodo DOI ``10.5281/zenodo.20489312``, version 1.0)
       for TS structure/performance studies
     - Redistribution rights are not established; method/backend mapping and
       candidate-bound energy, gradient, frequency, geometry, and IRC criteria
       are unresolved. Do not redistribute or mark release eligible without
       rights review and scientific approval.
   * - Microsolvation MS0 comparison
     - Finite candidate pools for methane, ammonia, formamide, and acetate with
       water; common-rescore comparisons are described in the protocol
     - Pilot protocol, not a numeric tolerance or bulk-solvation reference.
       Repeats and external comparator availability are limited. Not approved.
   * - Backend scan continuation records
     - Analytic potential tests and historical xTB/ORCA scan smoke outputs
     - Evidence of interface/continuation execution only; no approved scan
       energy/path tolerance or reaction mechanism reference.
   * - H2 GFN2-xTB derivative check
     - Historical xTB 6.7.1 result for neutral singlet H2 at 0.74 Angstrom:
       energy ``-0.981983694723 Eh``; analytic component
       ``-0.021144501361635 Eh/Bohr``; finite difference
       ``-0.021144504999682 Eh/Bohr`` at ``1e-4 Bohr`` displacement
     - A narrow real calculation with observed difference about
       ``3.64e-9 Eh/Bohr``. The record is not independently frozen in the
       current manifest and does not define an approved tolerance.
   * - H2 constrained scan smoke
     - Historical three-point xTB scan with finite energies
       ``-0.971086678722``, ``-0.980467429012``, and ``-0.981983694723 Eh``
     - Interface smoke only; not a reaction pathway, transition state,
       frequency, IRC, or benchmark approval.

Before adopting any row, verify the input/reference bytes and rights, reproduce
or independently establish the expected observable, define the exact method
and settings, and obtain scientific-owner review. The current machine-readable
catalogue at ``benchmarks/scientific_validation/catalogue.json`` remains the
inventory authority until a reviewed case manifest explicitly supersedes it.

Tolerance and acceptance policy
-------------------------------

There is no universal cross-method chemical tolerance. Each frozen case must
declare metric, unit, expected value, absolute and (where meaningful) relative
tolerance, backend/method/version, precision, settings, and justification.
Tolerance must reflect the intended scientific decision and numerical
reproducibility, not be widened to accommodate a failure.

* **Exact invariants:** atom symbols/order, charge, multiplicity, method,
  backend selection, frame count, and required workflow stages are exact.
* **Derivatives:** specify displacement, finite-difference formula, units,
  component or norm comparison, absolute floor, and scale. Relative error
  alone is unsuitable near zero.
* **Same-method numerical replay:** set energy, gradient, geometry, and
  convergence tolerances per method, backend/version, precision, thread count,
  and scientific purpose.
* **Structure:** declare mapping, alignment, fragment policy, RMSD unit, and
  identity/connectivity checks. Distance cutoffs alone cannot establish
  chemical equivalence.
* **Reaction path:** define barrier and unit, constraint error, path continuity,
  TS gradient, frequency significance, mode displacement, and endpoint
  identity. Require IRC connection when the claim concerns a reaction path.
* **Stochastic search:** declare seed, budget, duplicate rule, and ranking
  comparison. Acceptable alternative basins must be distinguished from hard
  regressions; discovery order need not be identical across platforms.
* **Optimization:** record optimizer termination and explicit convergence
  evidence separately from process exit. A terminated calculation is not
  necessarily converged.

Numerical values in unit tests are local test thresholds only. They become
release thresholds only through a reviewed, candidate-bound case manifest.

Backend validation matrix
-------------------------

Build this matrix from ``pyar/backend_capabilities.py`` and
``pyar/dependencies.py`` at each candidate. For every backend and supported
composition, report separately: installed, interface tested with fixtures or
mocks, actual numerical calculation completed, and unavailable/unsupported.
Never infer a real calculation from import discovery or mocked tests.

.. list-table:: Current evidence categories to verify per candidate
   :header-rows: 1

   * - Backend or interface
     - Software/interface evidence
     - Real-calculation requirement
   * - xTB GFN1/GFN2
     - Adapter, parser, unit, and failure tests
     - Pinned executable/version energy-gradient and optimization cases; record
       model and raw outputs
   * - xTB g-xTB
     - Capability detection and actionable unsupported-executable behavior
     - Model/executable identity and representative real calculation
   * - ORCA
     - Input generation, parser, provider and subprocess fixtures
     - Licensed runner with energy/gradient/optimization plus selected pathway
       checks
   * - Gaussian
     - Adapter and process-boundary tests
     - Authorized licensed calculation environment
   * - Turbomole / xTB-Turbomole
     - Adapter and failure-handling tests
     - Versioned real runs for declared supported methods
   * - Psi4, MOPAC, OpenBabel
     - Dependency, routing, parse and error behavior as implemented
     - Real calculations only for declared numerical support claims
   * - AIMNet2 / MLatom models
     - Dependency/model discovery and missing-asset diagnostics
     - Exact model asset/version/hash and numerical calculation
   * - ASE, geomeTRIC, Sella
     - Optimizer interface, state, cache, and termination tests
     - Real combinations with each declared supported calculator; verify
       convergence and geometry independently

Update actual availability and versions at execution time. Do not install
proprietary software without authorization. Optional backend absence is
``NOT_TESTED`` for numerical validation, not a pass.

Failure classification and report schema
-----------------------------------------

Classify each failed or unavailable case as one of: scientific regression,
numerical tolerance failure, backend incompatibility, infrastructure failure,
missing backend, unsupported workflow, nondeterministic result, or reference
data inconsistency. Every case record should contain:

* stable test/case identifier and workflow stage;
* candidate commit, input/reference hashes, charge, multiplicity, method,
  settings, and units;
* backend executable/model and version;
* expected value, observed value, tolerance, and comparator result;
* raw log/artifact paths and hashes;
* status ``PASS``, ``FAIL``, ``NOT_TESTED``, ``NOT_APPLICABLE``, or ``BLOCKED``;
* likely cause only when supported by recorded evidence.

Machine-readable JSON is authoritative; a human-readable report is generated
from the same data. The comparator must fail on missing mandatory evidence,
stale hashes, mismatched candidate SHA, malformed values, failed metrics, and
unauthorized exceptions. An explanation cannot override a failed result.

Execution and approval sequence
-------------------------------

#. Inventory implemented workflows, backend routes, fixtures, datasets,
   licenses, and existing real-run artifacts. Verify current source state
   rather than relying on historical status reports.
#. Select a small initial corpus using independent analytical, published, or
   independently verified reference evidence. Record source rights and exact
   input/reference hashes.
#. Write immutable manifests with method, electronic state, settings,
   observables, units, tolerance rationale, validation procedure, and budgets.
#. Add or reuse deterministic runners and comparators. Test positive results,
   wrong hashes, missing backend, malformed result, wrong state/method, and
   out-of-tolerance values.
#. Run Tier 1 on each change. Run Tier 2 against the exact candidate and retain
   workflow URL, raw artifacts, versions, and hashes. Schedule Tier 3 only
   after the scope and budget are approved.
#. Have a scientific reviewer inspect references, tolerances, raw results,
   pathway assignment, and exclusions. The scientific owner records approval
   or rejection; implementers do not approve on their behalf.
#. Freeze the approved corpus and candidate-bound evidence. Any changed bytes,
   method, backend version, candidate commit, tolerance, or exception requires
   a new review and applicable reruns.

Exit criteria
-------------

Scientific qualification for a stated release scope requires all of the
following:

* mandatory deterministic scientific tests pass;
* each required real-backend case has exact candidate and raw-artifact
  provenance and passes its reviewed criteria;
* required workflow coverage is executed or has an explicit, authorized
  exclusion;
* optimization/cache/convergence safety checks prevent failed or stale
  results from entering successful downstream workflows;
* references and tolerances are scientifically reviewed and immutable;
* reports distinguish pass, failure, not tested, not applicable, and blocked;
* remaining limitations and excluded claims are explicit.

Unexecuted mandatory evidence, unapproved references/tolerances, missing raw
artifacts, or unresolved safety defects prevent scientific qualification for
the affected claims. This plan itself is not evidence that any gate has
passed.

Related documentation
---------------------

The framework overview, benchmark catalogue, tolerance notes, provenance
rules, backend matrix, and troubleshooting guidance are in this directory's
``index``, ``benchmark_catalogue``, ``numerical_tolerances``,
``reference_provenance``, ``backend_validation``, and ``troubleshooting``
pages. The executable evidence validator is ``tools/scientific_validation.py``.
