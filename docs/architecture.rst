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

The eventual package structure should separate domain data, scientific
services, workflows, persistence, backend adapters, and user interfaces:

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

Public API
----------

The supported public API should be intentionally narrow:

.. code-block:: python

   from pyar import Molecule
   from pyar.sampling import generate_trials
   from pyar.selection import select_geometries
   from pyar.workflows import aggregate, react
   from pyar.backends import get_backend

Workflows and backend calls should exchange typed request and result objects:

* ``Molecule``
* ``OptimizationRequest`` and ``OptimizationResult``
* ``AggregateRequest`` and ``AggregateResult``
* ``ReactionRequest`` and ``ReactionResult``
* ``RunState``

Internal helper modules should remain internal rather than becoming accidental
user-facing interfaces.

Current Module Mapping
----------------------

The current codebase already contains useful boundaries. These should be
stabilized before modules are physically moved:

* ``pyar/Molecule.py`` becomes ``core/molecule.py``.
* ``pyar/molecule_geometry.py`` becomes ``core/geometry.py``.
* ``pyar/molecule_io.py`` becomes ``io/xyz.py``.
* ``pyar/orientation_sampling.py`` splits into ``sampling/sphere.py``,
  ``sampling/rotation.py``, and ``sampling/metrics.py``.
* ``pyar/trial_generation.py`` becomes ``sampling/trial_generator.py``.
* ``pyar/data_analysis/clustering.py`` splits into selection services.
* ``pyar/aggregator.py`` becomes ``workflows/aggregate.py``.
* ``pyar/reactor.py`` becomes ``workflows/reaction.py``.
* Legacy ``pyar/checkpt.py`` has been replaced for reaction workflows by
  ``pyar/reaction_state.py``; aggregation now uses ``pyar/aggregate_state.py``.
  Future workflow migration should converge on the structured ``state``
  package.
* ``pyar/interface/`` becomes ``backends/``.
* ``pyar/afir/restraints.py`` becomes ``biases/afir.py``.
* ``pyar/interface/xtb_turbo.py`` is replaced by a backend-neutral
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
   complete or trim by max-min diversity
   update basin memory
   write report

Each selection report should record the input count, rejected disconnected
geometries, duplicates removed, clusters found, selected count, relative
energies, algorithm, and thresholds.

Command Line Interface
----------------------

The long-term CLI should use subcommands:

.. code-block:: bash

   pyar aggregate ...
   pyar react ...
   pyar optimize ...
   pyar cluster ...
   pyar benchmark sampling ...

The existing command entry points may remain as temporary wrappers during
migration. The CLI should print a run plan before execution, validate backend
compatibility early, report missing external programs clearly, and support
``--dry-run``, ``--resume``, ``--output``, and ``--log-level``.

Packaging And Quality
---------------------

The professional package baseline should include:

* ``pyproject.toml`` as the authoritative metadata source
* a later move to a ``src/pyar/`` layout in a major-version branch
* optional dependency groups such as ``selection``, ``xtb``, ``aimnet2``,
  ``docs``, and ``test``
* typed public APIs
* ``pytest`` for testing and ``ruff`` for formatting and linting
* optional static type checking of public modules
* Sphinx user guides and API references
* executed examples in continuous integration
* semantic versioning, a changelog, and a deprecation policy

Bundled code such as ``pyar/mlatom/`` should either become an explicitly
maintained vendored dependency or be replaced by an optional external adapter.

Migration Plan
--------------

The migration should be incremental and behaviour-preserving:

1. Freeze and commit the currently validated baseline.
2. Define structured request, result, capability, and exception types without
   moving existing modules.
3. Introduce energy/gradient provider, bias-potential, and optimizer
   interfaces; verify AFIR energy and gradients with unit and finite-difference
   tests.
4. Extract an xTB energy/gradient provider independent of Turbomole and
   integrate the TRIC optimizer for biased reaction jobs.
5. Extend the structured run-state approach now implemented for reaction and
   aggregation workflows to other long workflows before making new engines
   the default.
6. Benchmark the new optimizer against current reaction behaviour and
   Cartesian baselines, then remove the Turbomole requirement from the xTB
   reaction path.
7. Split the selection implementation into focused services while preserving
   selected outputs.
8. Convert workflows to return structured result objects.
9. Introduce a subcommand CLI while retaining current wrappers for one
   transition release.
10. Move modules to the target layout only in a major-version development
    branch.
11. Decide whether to vendor or externalize MLatom.
12. Publish migration documentation and reproducibility examples.

Immediate Priority
------------------

The next architectural work consists of two coupled priorities:

* continue expanding the versioned structured restart state implemented for
  reaction and aggregation runs to solvation and other long workflows
* Turbomole-free internal-coordinate reaction optimization, beginning with a
  tested xTB energy/gradient provider and AFIR objective

Both must be in place before expanding reaction functionality: the optimizer
must be scientifically appropriate, and a long calculation must be safely
restartable.
