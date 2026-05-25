Biased Reaction Optimization
============================

This design defines the replacement for the current Turbomole-dependent AFIR
optimization path in ``react``. Reaction optimization must be an independent
PyAR service that can use xTB, machine-learning potentials, or quantum
chemistry backends without assigning coordinate-update policy to a backend.

Current Problem
---------------

The current route through ``pyar.interface.xtb_turbo`` evaluates an xTB
gradient, adds the AFIR contribution, writes the modified energy and gradient
into Turbomole files, and asks Turbomole to update the coordinates and check
convergence. This has three design consequences:

* an xTB/AFIR reaction calculation requires a Turbomole installation
* optimization strategy is hidden inside a backend wrapper rather than owned
  by the reaction workflow
* the current ``xtb_turbo`` and Turbomole AFIR implementations set
  ``gamma = 100.0`` internally, instead of using the gamma schedule supplied
  by ``react``

``pyar.afir.restraints.isotropic`` already provides the important reusable
piece: additive AFIR energy and Cartesian force evaluation. The geomeTRIC
adapter uses that force directly instead of coupling it to Turbomole file
rewriting.

Target Objective
----------------

The reaction workflow constructs one objective from a physical energy/force
provider and zero or more additive bias potentials:

.. math::

   E_\mathrm{total}(x) =
   E_\mathrm{method}(x) + E_\mathrm{AFIR}(x)

.. math::

   F_\mathrm{total}(x) =
   F_\mathrm{method}(x) + F_\mathrm{AFIR}(x)

The optimizer receives only this objective. The implemented providers are xTB,
AIMNet2, ORCA, and Gaussian; other engines require energy/force
adapters before they can use this channel.
The AFIR gamma value belongs to the bias settings stored in the request and
run state; it must never be replaced by a wrapper-level constant.

Core Interfaces
---------------

The initial contracts should be small and testable:

.. code-block:: python

   @dataclass(frozen=True)
   class EnergyGradientResult:
       energy_hartree: float
       gradient_hartree_per_bohr: np.ndarray

   class EnergyGradientProvider(Protocol):
       def evaluate(
           self, molecule: Molecule, coordinates_bohr: np.ndarray
       ) -> EnergyGradientResult: ...

   class BiasPotential(Protocol):
       def evaluate(
           self, molecule: Molecule, coordinates_bohr: np.ndarray
       ) -> EnergyGradientResult: ...

   class GeometryOptimizer(Protocol):
       def minimize(
           self, objective: Objective, molecule: Molecule, settings: OptimizerSettings
       ) -> OptimizationResult: ...

``Objective`` sums the provider and bias results, counts evaluations, logs
physical and bias energies separately, and records all optimizer steps for
restart and diagnosis.

The current registry wires this contract through backend-specific provider
factories. The provider owns backend settings such as charge, multiplicity,
and executable details; the objective only sees the result object above.

Units are part of this contract:

* coordinates at the objective boundary: Bohr
* energies: Hartree
* Cartesian forces or gradients: Hartree/Bohr, with sign translated
  explicitly at the optimizer boundary
* stored XYZ and user-facing geometry coordinates: Angstrom

Conversions should occur once at the objective boundary. This aligns the new
objective with the units already expected by the current AFIR calculation.

Internal-Coordinate Optimizer Choice
------------------------------------

Cartesian ``LBFGS`` and ``FIRE`` are useful baselines and fallbacks, but they
do not satisfy the primary design requirement: efficient molecular
optimization in internal coordinates. The implemented xTB, AIMNet2, ORCA,
and Gaussian AFIR reaction channels use **geomeTRIC with TRIC
coordinates**.

TRIC adds explicit translation and rotation coordinates for each molecular
fragment alongside intrafragment internal coordinates. This matches an AFIR
reaction workflow in which independently oriented fragments approach and
rearrange while covalent bonding changes. GeomeTRIC also accepts a custom
energy/Cartesian-gradient engine, using Hartree and Hartree/Bohr units, which
matches the objective contract above.

Candidate implementations should be treated as follows:

.. list-table::
   :header-rows: 1
   :widths: 22 32 46

   * - Optimizer
     - Coordinate model
     - PyAR role
   * - geomeTRIC/TRIC
     - Delocalized intrafragment internals plus explicit fragment translation
       and rotation coordinates
     - Proposed default for biased reaction minimization and the first
       integration target
   * - PyBerny
     - Redundant internal coordinates with RFO, trust radius, and Hessian
       updates
     - Secondary Python benchmark for bonded molecular cases; not the default
       for fragment-association paths
   * - DL-FIND
     - DLC and HDLC, plus Cartesian modes
     - Strong reference implementation and optional later integration,
       particularly for transition-state extensions
   * - xTB ANCopt/L-ANCopt
     - Approximate normal coordinates internal to xTB
     - Native unbiased xTB benchmark; usable for AFIR only if an additive
       external-gradient integration is demonstrated
   * - ASE LBFGS/FIRE
     - Cartesian steps, optionally preconditioned
     - Simple fallback and performance baseline; not the AFIR production
       default

ASE remains valuable as a calculator interoperability layer. For example, an
AIMNet2 or other ASE-compatible potential can be adapted to the common
energy/gradient provider. It should not force the reaction workflow to use a
Cartesian ASE optimizer.

Proposed Package Boundary
-------------------------

The new core should be built under these logical boundaries before modules are
physically moved:

.. code-block:: text

   biases/
     base.py             # BiasPotential protocol
     afir.py             # AFIR energy and gradient implementation
   optimizers/
     base.py             # OptimizerSettings and GeometryOptimizer protocol
     objective.py        # Provider plus bias composition and logging
     geometric.py        # geomeTRIC/TRIC adapter
     berny.py            # optional benchmark adapter
     ase_cartesian.py    # LBFGS/FIRE fallback and benchmark adapter
     dlfind.py           # optional future integration
   backends/
     xtb.py              # xTB energy/gradient provider and unbiased shortcut
     aimnet2.py          # ML energy/gradient provider
     ...
   workflows/
     reaction.py         # gamma schedule, product tests, state, reporting

The intended execution path is:

.. code-block:: text

   reaction request
     -> trial geometry and gamma schedule
     -> provider = xTB/AIMNet2/QC energy and gradient
     -> objective = provider + AFIR(gamma)
     -> geomeTRIC/TRIC optimizer
     -> product identity analysis
     -> structured state and reports

Current Restart Foundation
--------------------------

Reaction calculations now persist versioned state in ``reaction/state.json``
and geometry snapshots under ``reaction/state/geometries/``. The state records
the numeric gamma schedule, pending and retained current-cycle geometries,
completed jobs, products, and invocation settings. Resume is rejected if the
new invocation changes scientific inputs or backend settings.

Legacy ``jobs.pkl`` checkpoints are migrated only when their old gamma labels
map unambiguously to the requested numeric schedule. This gives the later
internal-coordinate optimizer an inspectable and atomic restart foundation.

Implemented Channel And Next Steps
----------------------------------

For ``react`` calculations using ``--software xtb``, ``--software aimnet_2``,
``--software orca``, or ``--software gaussian``, PyAR now selects the
geomeTRIC/TRIC optimizer and
evaluates the physical energy/forces together with the AFIR force. When a
bonded candidate is found, product relaxation is repeated without AFIR bias
(``gamma=0.0``). Native optimizers remain the default for aggregate and
ordinary backend optimization workflows.

The remaining development work is to store optimizer step and objective
component details in restart state, benchmark this channel on controlled
reaction cases, add quantum-chemistry force-provider adapters, and implement
transition-state optimization as a separate validated workflow.

Validation And Benchmarks
-------------------------

The default cannot be chosen solely from algorithm descriptions. The
integration must report:

* AFIR analytic gradient agreement with finite differences
* final energy and product identity for each successful optimization
* energy/gradient evaluation count, optimizer step count, wall time, and
  convergence reason
* results under multiple gamma values and multiple initial orientations
* restart equivalence after an interrupted optimization
* failure behavior for disconnected fragments, failed backend evaluations,
  and non-convergent trials

The initial benchmark matrix should compare:

* geomeTRIC/TRIC with xTB plus AFIR
* PyBerny with the same objective where fragment handling permits
* ASE LBFGS and FIRE with the same objective as Cartesian baselines
* native xTB ANCopt for unbiased reference optimizations
* DL-FIND/DLC or HDLC as a later reference if installation and integration are
  practical

Algorithm References
--------------------

The implementation should cite the algorithm that is actually used:

* L.-P. Wang and C. Song, *Geometry optimization made simple with translation
  and rotation coordinates*, Journal of Chemical Physics **144**, 214108
  (2016), DOI ``10.1063/1.4952956``. This defines TRIC and geomeTRIC.
* J. Kastner, J. M. Carr, T. W. Keal, W. Thiel, A. Wander, and P. Sherwood,
  *DL-FIND: An Open-Source Geometry Optimizer for Atomistic Simulations*,
  Journal of Physical Chemistry A **113**, 11856-11865 (2009), DOI
  ``10.1021/jp9028968``. This documents DLC/HDLC-capable optimizer
  infrastructure.
* C. Peng, P. Y. Ayala, H. B. Schlegel, and M. J. Frisch, *Using redundant
  internal coordinates to optimize equilibrium geometries and transition
  states*, Journal of Computational Chemistry **17**, 49-56 (1996), DOI
  ``10.1002/(SICI)1096-987X(19960115)17:1<49::AID-JCC5>3.0.CO;2-0``.
  This is a direct basis for redundant-internal-coordinate Berny-style
  approaches.

The new implementation should integrate and test an established
internal-coordinate optimizer rather than reimplement the Pulay/Schlegel
coordinate transformation machinery from scratch.
