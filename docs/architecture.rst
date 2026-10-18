Architecture Roadmap
====================

This document defines the intended long-term design of PyAR as a
professional scientific Python package. The redesign should preserve verified
scientific behaviour while progressively clarifying APIs, workflow state,
backend boundaries, and reproducibility.

Design Goals
------------

* Provide a small, stable Python API, with the CLI as a thin frontend.
* Make geometry generation, optimization, selection, and restart behaviour
  reproducible and inspectable.
* Isolate external calculation engines behind explicit backend adapters.
* Make biased internal-coordinate optimization a core service for reaction
  workflows, without requiring Turbomole to optimize xTB or ML potentials.
* Persist structured run state so interrupted calculations can resume safely.
* Keep heavy optional dependencies tied to the functionality that requires
  them.
* Support tested documentation, versioned interfaces, and a clear
  deprecation policy.

Target Package Layout
---------------------

The package structure target separates domain data, scientific services,
workflows, persistence, backend adapters, and user interfaces. This section is
a design map, not a promise that every listed target module exists:

.. code-block:: text

   pyar/
     core/
       molecule.py
       geometry.py
       elements.py
       exceptions.py
       types.py
     io/
       xyz.py
       results.py
       run_directory.py
     sampling/
       sphere.py
       rotation.py
       trial_generator.py
       metrics.py
     selection/
       deduplication.py
       clustering.py
       diversity.py
       basin_memory.py
       reports.py
     structure_comparison/
       models.py
       identity.py
       equivalence.py
       legacy.py
     backends/
       base.py
       registry.py
       xtb.py
       orca.py
       gaussian.py
       psi4.py
       turbomole.py
       mopac.py
       aimnet2.py
       mlatom.py
     biases/
       base.py
       afir.py
     aggregation/
       request.py
       growth.py
       pathways.py
       reporting.py
     solvation/
       request.py
       growth.py
       reporting.py
     conformer/
       request.py
       generation.py
       selection.py
       reporting.py
     optimizers/
       base.py
       objective.py
       geometric.py
       berny.py
       ase_cartesian.py
       dlfind.py
     workflows/
       aggregate.py
       reaction.py
       scan.py
       solvation.py
       explore.py
     state/
       models.py
       store.py
       restart.py
     cli/
       main.py
       options.py
       reporting.py
   diagnostics/
     logging.py
     benchmark_sampling.py
     benchmark_selection.py

The current package has canonical modules under ``pyar.core``, ``pyar.io``,
``pyar.sampling``, ``pyar.state``, ``pyar.backends``, ``pyar.aggregation``,
``pyar.solvation``, and ``pyar.conformer``. Use the generated API reference and
module documentation for the current implementation; the target tree above is
illustrative.

Feature Package Boundaries
--------------------------

Feature-specific code should live in feature packages rather than becoming the
shape of the whole project. ``pyar.conformer`` is a conformer-search package,
``pyar.aggregation`` is an aggregation package, and ``pyar.solvation`` is a
solvation and ligand-growth package; none should become the main architecture
for every workflow. Reaction exploration and mechanism-finding, clustering,
similarity, bonding analysis, and backend interfaces should each keep their own
request models, services, and tests where they have domain-specific behaviour.

Structure comparison separates chemical identity, geometrical equivalence,
and diversity selection. Identity providers and geometry comparators live in
``pyar.structure_comparison``; selection continues to own candidate reduction
and clustering. See :doc:`structure_identity_and_selection` for the
behavior-preserving contract.

Workflow modules such as ``pyar.workflows.conformer`` should remain thin public
orchestration layers. They may translate CLI/API arguments into typed feature
requests, call feature services, persist workflow state, and return public
result objects, but feature-specific validation and policy should be owned by
the relevant feature package.

Public API
----------

The supported public API should be intentionally narrow:

.. code-block:: python

   from pyar import Molecule
   from pyar.sampling import generate_trial_vectors
   from pyar.selection import choose_geometries
   from pyar.workflows import aggregate, react
   from pyar.backends import get_backend_capabilities

Workflows and backend calls should exchange typed request and result objects:

* ``Molecule``
* ``OptimizationRequest`` and ``OptimizationResult``
* ``AggregateRequest`` and ``AggregateResult``
* ``ReactionRequest`` and ``ReactionResult``
* ``RunState``

Internal helper modules should remain internal rather than becoming accidental
user-facing interfaces.

Module Migration Map
----------------------

This migration map describes the intended paths; use the imports above and
generated API reference for the current implementation:

* ``pyar/Molecule.py`` becomes ``core/molecule.py``.
* ``pyar/molecule_geometry.py`` becomes ``core/geometry.py``.
* ``pyar/molecule_io.py`` becomes ``io/xyz.py``.
* ``pyar/orientation_sampling.py`` splits into ``sampling/sphere.py``,
  ``sampling/rotation.py``, and ``sampling/metrics.py``.
* ``pyar/trial_generation.py`` becomes ``sampling/trial_generator.py``.
* ``pyar/data_analysis/clustering.py`` splits into selection services.
* ``pyar/aggregator.py`` becomes ``workflows/aggregate.py``.
* ``pyar/reactor.py`` has moved to ``workflows/reaction.py``.
* Legacy ``pyar/checkpt.py`` has been replaced for reaction workflows by
  ``pyar/reaction_state.py``; aggregation now uses ``pyar/aggregate_state.py``
  and solvation uses ``pyar/solvation_state.py``. Future workflow migration
  should converge on the structured ``state`` package.
* ``pyar/interface/`` has moved to ``backends/``.
* ``pyar/afir/restraints.py`` has moved to ``biases/afir.py``.
* ``pyar/backends/xtb_turbo.py`` is replaced by a backend-neutral
  reaction-optimization service, not moved into the new package as a
  permanent adapter.

Core Domain Rules
-----------------

``Molecule`` should represent molecular geometry and minimal workflow
metadata:

* atom symbols and Cartesian coordinates
* charge, multiplicity, and optional energy
* fragment definitions needed by aggregate and reaction workflows
* lightweight provenance metadata

It should not own:

* restart or output-directory policy
* backend command construction
* optimization strategy
* clustering or selection policy

Geometry transforms should prefer explicit non-mutating APIs such as
``translated()``, ``rotated()``, and ``merged_with()``. Existing in-place
operations may remain during migration where workflow compatibility requires
them.

Workflow Model
--------------

Every workflow should follow the same visible lifecycle:

.. code-block:: text

   validate request
   create or resume run directory
   generate candidates
   evaluate candidates through a backend
   select survivors
   persist state and reports
   return a result object

Aggregate, reaction, scan, solvation, and exploration workflows should reuse
common services for optimization, logging, selection, and persistence rather
than reimplementing these behaviours independently.

Backend Contract
----------------

External calculation methods should provide energies and Cartesian gradients
independently of the algorithm that updates coordinates. A backend may retain
a direct optimization shortcut for unbiased calculations, but reaction
workflows must use the composable objective:

.. code-block:: python

   @dataclass(frozen=True)
   class EnergyGradientResult:
       energy_hartree: float
       gradient_hartree_per_bohr: np.ndarray

   class EnergyGradientProvider(Protocol):
       def evaluate(self, molecule, coordinates_bohr) -> EnergyGradientResult: ...

   class BiasPotential:
       def evaluate(self, molecule, coordinates_bohr, settings) -> EnergyGradientResult: ...

   class GeometryOptimizer:
       def minimize(self, objective, molecule, settings) -> OptimizationResult: ...

Capabilities should report whether a backend supports:

* loose, normal, and tight optimization stages
* method and basis configuration
* Cartesian energy and gradient evaluation
* charge and multiplicity
* parallel execution
* an external executable or an optional Python dependency

The CLI should validate requested options against capabilities. Workflow
modules should operate on the generic adapter contract rather than branch on
software names.

Biased Optimization Service
---------------------------

Reaction-path optimization is a first-class core service, not a backend
special case. The current xTB/AFIR route calls xTB for energies and gradients
but uses Turbomole to update coordinates and determine convergence. The new
architecture removes this dependency by composing an electronic-structure
gradient with an AFIR bias before passing it to an internal-coordinate
optimizer.

The detailed optimizer selection, objective contract, migration, and benchmark
requirements are defined in :doc:`reaction_optimization`.

Restart And Output State
------------------------

Restart reliability is the highest-priority architectural improvement.
Reaction workflows use ``reaction/state.json`` together with XYZ snapshots
instead of mutable pickle checkpoints. Aggregation workflows now use
``aggregates/state.json`` to validate the request, persist pathway order,
record completed pathways and selected outputs, and resume interrupted
pathways without relying on log parsing. Other workflows should converge on
the same principle: JSON metadata together with XYZ geometries and
human-readable summaries:

.. code-block:: text

   run_name/
     run.json
     pyar.log
     pathways/
       path_000/
         state.json
         C2H6/
           candidates/
           selected/
     selected/
       C2H6/
         structures/
         energies.csv
         selection.json

``run.json`` should store:

* package version and invocation
* inputs and backend parameters
* sampling configuration and deterministic sequence identifiers
* workflow stages and completed jobs
* paths to selected results

Restart rules:

* Resume only when the saved configuration is compatible with the request.
* A current pathway selection replaces the prior current selection atomically.
* Historical selections are stored separately from current selection inputs.
* Final cross-path selection reads only completed, current pathway outputs.

Trial Generation
----------------

Trial geometry generation is a primary scientific component. The intended
default is:

* spherical Fibonacci approach directions
* quaternion-based rotations for multi-atom monomers
* deterministic indexed variants when multiple populations are requested
* explicit metrics for sphere and rotational coverage

Candidate generation should remain separate from contact placement and
optimization. Future work may add symmetry-aware orientation reduction,
reactive-site constraints, and adaptive sampling based on discovered basins.

Selection Pipeline
------------------

Selection should be expressed as a composable, reportable pipeline:

.. code-block:: text

   validate structures
   filter disconnected candidates
   deduplicate geometry
   extract cluster minima
   trim excess cluster minima by max-min diversity
   update basin memory
   write report

The cluster minima are the candidate seed set. If clustering produces fewer
minima than the requested seed limit, selection keeps those minima without
filling the remaining slots from other geometries. Max-min diversity is used
only to trim the cluster minima when they exceed the limit.

Each selection report should record the input count, rejected disconnected
geometries, duplicates removed, clusters found, selected count, relative
energies, algorithm, and thresholds.

Command Line Interface
----------------------

The current command-line interfaces are ``pyar`` for task-oriented commands and
``pyar-cli`` for legacy syntax. See :doc:`cli_reference` for the registered
console scripts. Command-specific ``--help`` output is authoritative. The modern
interface provides read-only ``--check`` and ``--dry-run`` preflight controls;
not every legacy command supports every modern control.

Packaging And Quality
---------------------

Current packaging metadata is maintained in ``pyproject.toml``. Runtime
capabilities are declared as optional extras; checkout tools are declared as
dependency groups. CI runs pytest, scoped Ruff checks, strict Sphinx, archive
validation, and clean artifact installation. See :doc:`installation` for the
commands and :doc:`cli_reference` for installed entry points.

PyAR uses the optional external ``mlatom`` package for MLatom calculations.
The ``pyar/mlatom/__init__.py`` adapter preserves older package-level imports
by forwarding them to the external package. The vendored source, binaries,
and model assets have been removed; the adapter requires MLatom to be installed
with ``python -m pip install "pyar-chem[ml]"``.

Migration and implementation status
------------------------------------

This page began as an architecture roadmap. Its target designs are not a live
implementation or release-status ledger. The current public API and compatibility
paths are described in :doc:`migration`, :doc:`migration_2_0`, and the generated
:doc:`api` reference. Release stages and remaining engineering gates are tracked
in the repository's release-engineering plan.

Shared single-addition service
------------------------------

``pyar.growth.service.add_one`` is the shared orientation/optimization/selection
engine. ``pyar.workflows._growth`` remains a compatibility alias. Aggregate
retains composition and pathway orchestration and its existing restart schema;
grow owns fixed-seed sequential orchestration, a typed ``GrowRequest``, atomic
``GrowRunState``, intermediate pools, and ``GrowResult``. Solvation continues to
use the same service with unchanged frontend semantics, including its legacy
geometry-only propagation. Aggregate and grow geometry-only pools
are bounded using existing graph-first pruning and max-min diversity, without
energy assignment.

Run-layout persistence
----------------------

The :doc:`persistence` module provides version detection, artifact references,
JSON publication, and local writer locking. Workflow owners adopt these helpers
according to their documented continuation contracts. See :doc:`run_layout` for
current behavior; the architecture target does not imply that every workflow
uses the same state format.
