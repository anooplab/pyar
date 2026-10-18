ORE implementation status after the R1–R4 review
================================================

Scientific scope: **frozen, version 1.1, 2026-10-18**, by project direction.
This freeze covers the scientific contract below, including the current
review fixes and explicit workflow opt-in. It does not assert that uncommitted
code has been merged or that every feature in the original broader plan exists.

Version 1.1 amendment: the project owner explicitly requested single-molecule
bond scans after the initial freeze. One input and two indices now select an
intramolecular formation/breaking scan, deterministic distance/grid planning,
unconstrained final relaxation and canonical product comparison. Opt-in
scan-seeded path following is restricted to identified changed products and
retains every existing NEB, TS, frequency, IRC and endpoint gate. The scientific
contract is specified in :doc:`bond_scan`; this amendment adds no automatic
recovery strategy or benchmark requirement.

The R1–R4 review found two implementation defects and a mismatch between the
original full plan and the delivered pilot. The reporting defects are repaired:
ordinary diagnostic I/O failures preserve the owner's result and original
exception, a missing report is represented by null, and cancellation propagates
through execution, inspection and report writing.

Scope
-----

This document tracks implementation and software checks. Benchmark campaigns,
paired cost/outcome measurements and chemistry-performance assessment are a
separate follow-up project, excluded from implementation completion criteria.
Their absence does not leave this implementation scope open.

Frozen scientific contract
--------------------------

* Automatic recovery is explicitly opted in and restricted to unbiased,
  closed-shell GFN2-xTB minima through the verified xTB 6.7.1 / geomeTRIC 1.1.1
  route. Defaults remain off. The eligible workflow roles listed below remain
  fixed; other routes retain diagnosis and their existing owner gates.
* SAFE permits bounded cycle continuation from verified accepted frames when
  energy and gradient evidence show progress. GUIDED additionally permits the
  existing plateau-triggered rollback, fresh Hessian and reduced trust radius.
  No additional recovery strategy is part of this scope.
* Hamiltonian, charge, spin/reference state, atom ordering, physical objective,
  constraints and requested convergence thresholds remain invariant. A retry
  cannot weaken an acceptance criterion or substitute a minimum for a TS.
* Attempt, cycle, evaluation and wall-time budgets remain finite. Immutable
  accepted evidence, request fingerprints and explicit compatible resume own
  restart decisions. Diagnostic files and failed observed paths do not.
* Workflow owners retain composition, molecular/fragment identity, ranking,
  shell/site acceptance, unbiased-energy, scan-distance, saddle/Hessian,
  frequency, IRC and endpoint-connection gates. Provisional acceptance remains
  distinct from numerical convergence and final chemistry acceptance.
* Unsupported controls, unverified electronic-root continuity and missing
  accepted-step contracts stop recovery. Diagnose-only observations cannot
  certify convergence.

Bug fixes, implementation cleanup, documentation and software regression checks
may continue within this contract. Any change to methods, recovery actions,
eligibility, scientific thresholds, objectives or chemistry gates requires a
separate proposal and explicit project-owner approval of a new scope version.
The freeze is a scope decision, not a performance or chemistry-validation claim.

Delivered code
--------------

* R1 retains the verified ORCA/xTB/geomeTRIC parser matrix. Other advertised
  backends now have explicit diagnose-only observations: artifacts, original
  exception chains, process failures and backend-specific negative markers.
  Neither a zero exit nor a composite method's first-stage success certifies
  convergence. Unknown formats remain unvalidated.
* R2 retains bounded accepted-frame cycle continuation, immutable attempt
  records, objective fingerprints, evaluation/cycle/wall budgets and explicit
  resume. Only independently converged outcomes publish final XYZ files.
* R3 delivers the guided rollback/trust/Hessian restart and the bounded
  finite-difference diagnostic helper. Broader actions are future extensions
  outside the frozen scope.
* R4 records each owner's explicit permission to use incomplete results,
  separately from numerical convergence and final workflow chemistry gates.
  The modern optimize report distinguishes validated and provisional counts;
  its historical optimized count retains its legacy meaning.
* Constrained native/generic scans now have typed execution/exception evidence
  and independent distance, atom mapping, finite-energy and point-count checks,
  including cached scan adoption. Invalid claimed successes are rejected using
  the existing 0.001 angstrom constraint tolerance.
* NEB, TS, frequency, IRC and endpoint summaries contain separate owner gates
  and original exception chains on failure or cancellation. Available failed
  band/geomeTRIC/Sella observations are retained as diagnostic XYZ files.
  These can contain unaccepted trials and are never restart authority. A TS
  optimization boolean does not pass the saddle-Hessian or IRC connection gate.

Explicit workflow opt-in
------------------------

The Python workflow adapter accepts ``workflow_recovery`` in its QC settings.
For example, an owner of an unbiased xTB/geomeTRIC GFN2 minimum may supply:

.. code-block:: python

   qc_params["workflow_recovery"] = {
       "mode": "safe", "max_attempts": 2, "max_cycles": 200,
       "max_evaluations": 400, "wall_seconds": 600,
   }

The adapter removes these execution controls before invoking a backend and
uses the existing recovery controller. Only a converged accepted frame updates
the owner's molecule. This is a Python API opt-in; no workflow CLI default or
selection ranking is changed. The original input and explicit ``resume: True``
are required for resume. Checkpoints preserve the workflow role and reject
cross-owner reuse. Read-only preflight uses ``preflight_recovery`` with that
same workflow role and policy before invoking the adapter.

Eligible roles are conformer refinement, growth preselection/refinement,
microsolvation unbiased screening/minima/postrelaxation, reaction unbiased
relaxation and bond-scan unbiased postrelaxation. They retain the pilot's
closed-shell, GFN2, optimizer-version and physical-control restrictions.
Walls, AFIR/release retries, constraints, TS/path objectives and shared
microsolvation evaluation ledgers have no recovery permission. An explicit
unsupported request fails before a backend starts; it never silently falls
back to an unbudgeted legacy execution. Omit the option to use existing behavior.

Future extensions outside the frozen scope
------------------------------------------

These original-plan features are excluded from version 1.1. They are future
proposals, not open requirements for completion of the frozen implementation:

* R3 SCF correction ladders: backend/version-specific controls and verifiable
  electronic-root continuity are required before enabling solver changes.
* R3 optimizer-family switching and internal-coordinate rebuild strategies:
  compatible accepted-frame/Hessian, constraints and objective contracts are
  required for each route. The current fresh geomeTRIC restart is narrower.
* R3 automatic gradient-consistency diagnosis: the existing 6N-evaluation
  helper still needs integration with each provider's immutable cost ledger.
* R3 guided collision repair/preconditioning: workflow authorization, identity
  checks and return to the original Hamiltonian must precede enablement.
* R4 constrained scan recovery and NEB/image/Sella/IRC recovery: current
  integration diagnoses failures and retains observations; it cannot restart
  from those observations. Certified accepted-step/band checkpoints and
  compatible owner-controlled replay are still required.
* R4 remaining backend full validation/recovery: diagnose-only evidence does
  not replace versioned native convergence, gradients, accepted frames or
  electronic checkpoint compatibility.

The original broad R3/R4 plan remains a historical roadmap. Completion of this
frozen scope is assessed against the version 1.1 contract, delivered code and
software checks above. Wider recovery permissions are not implemented and must
not be described as available.
