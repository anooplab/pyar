Optimization diagnosis (ORE R1)
===============================

The opt-in Python service diagnoses one minimum-optimization execution using
independent native output checks. Existing CLI commands and workflow APIs keep
their current behavior. This stage adds no retries, fallback methods, cache
adoption, resume API, or change to scientific defaults.

R0 review
---------

The R0 audit at commit ``ecaf7d5`` was reviewed against the private ORE plan,
its source inventory and adversarial fixtures. No blocking issue was found for
starting observational R1. All 52 reported R0 check runs passed at that head.
The audit identified separate native, geomeTRIC, Sella, NEB/IRC, scan, RDKit,
ASE/ML and workflow boundaries; a wrapper around ``optimise`` cannot cover all
of them. The 17 strict expected failures continue to characterize legacy risks.
They are not silently fixed or removed in R1.

The review accepts typed distinctions between numerical convergence, accepted
incomplete geometry, execution termination, cache provenance, and workflow
scientific validity. The new facade validates fresh executions separately;
legacy status helpers, continuation semantics and committed records stay intact.
This implements R1 only. R0b correctness changes and R2 recovery remain separate
review decisions. The earlier uncommitted retry prototype is excluded.

Single-execution Python API
---------------------------

Use resolved backend settings and perform the selected workflow's normal input,
spin and dependency preflight before execution. For a validated neutral closed-
shell XYZ input, for example:

.. code-block:: python

   from pyar.core.molecule import Molecule
   from pyar.optimization_request import preflight
   from pyar.optimization_service import optimize_once

   molecule = Molecule.from_xyz("input.xyz")
   molecule.name = "molecule"
   settings = {
       "software": "xtb", "xtb_model": "gfn2",
       "geometry_optimizer": "native", "opt_target": "minimum",
       "opt_threshold": "normal", "nprocs": 1,
   }
   preflight(settings, [molecule])
   result = optimize_once(molecule, settings)
   print(result.to_json())
   if result.converged:
       frame = result.accepted_frame
       print(frame.physical_energy_hartree)

``optimize_once`` calls the existing optimizer exactly once on a copy of the
molecule and settings, inside its existing ``job_<name>`` convention. The caller's
molecule is unchanged. Noninteger or conflicting molecule/QC charge, spin or
xTB unpaired-electron requests are rejected before execution. It refuses any existing job directory, obtains the
existing single-writer ``RunLock`` and never adopts a cached XYZ as newly
converged. No additional process-wide working-directory layer is introduced;
the legacy optimizer still changes cwd and execution is not thread-safe.

The execution receipt contains observed subprocess return codes, exact request
fingerprint, artifact hashes, wall time, legacy status and any caught native
process exception. Original wrapper exceptions and their chains are retained
on ``result.exception``; native caught exceptions are also accessible through
``attempt.execution.processes``. No diagnostic hook changes command arguments,
exit status, file naming, or the legacy return value. When no observer is active,
the process hooks perform no recording or I/O.

Readonly inspection
-------------------

.. code-block:: python

   from pyar.optimization_service import inspect_optimization

   attempt = result.attempts[0]
   checked = inspect_optimization(
       attempt.request, "job_molecule", execution=attempt.execution,
   )
   assert checked.converged == result.converged

Inspection reads files and hashes them without creating directories, locks,
reports or markers. Artifacts are checked before and after parsing; malformed,
truncated, modified, missing or symlinked evidence fails closed. Repeat inspection
of unchanged evidence serializes deterministically.

Without a matching fresh-execution receipt, inspection can still report observed
failures from historical logs, but marks the context unverified and cannot certify
a current converged result or publish an accepted frame. A receipt supplied by a
Python caller is an explicit attestation by that execution owner, not a signature
or permission to trust arbitrary user-supplied reports. R1 does not deserialize
reports into restart authority or adopt old workflow files.

Typed results and evidence
--------------------------

``OptimizationRequest`` freezes ordered atoms, input coordinates, charge, spin,
calculator/optimizer/target, workflow role, requested objective/constraints and
settings. Its scientific fingerprint excludes only parallelism and iteration
budgets. The exact request fingerprint includes those numerical controls;
requested convergence thresholds remain in both identities. Known model/wrapper
files use hashes captured at request creation rather than machine paths in these
fingerprints. The raw settings retain their execution provenance.

``OptimizationResult`` and ``OptimizationAttempt`` retain diagnostics, original
request, versions, execution receipt, source paths and SHA-256 artifact refs.
``to_json()`` uses canonical finite JSON with ``schema_version: 1`` and a machine-
readable reason. Units and meanings are explicit: physical and effective energies
are Hartree, coordinates are Angstrom, effective forces are eV/Angstrom, and
optimizer gradient observations are Hartree/Bohr. xTB's total gradient norm is
recorded as a norm, not mislabeled as an RMS gradient.

The result records backend termination, SCF observations, optimizer convergence,
stationarity, physical/effective objective validity, constraints and workflow
scientific gates separately. An unknown gate stays unknown. SCF root identity,
frequency/Hessian index, conformer identity, shell/site acceptance and reaction
connection remain the owning workflow's responsibility. ``converged`` means
validated numerical convergence for the requested objective;
``workflow_accepted`` remains false in R1. It is not proof of a physical minimum,
a first-order saddle or a connected reaction pathway.

* ``CONVERGED`` requires genuine optimizer evidence, finite matching frame and
  energy, original gradient thresholds, and observed normal process completion.
* ``INCOMPLETE_USABLE`` requires independently verified accepted geometry. It
  never counts as numerical convergence. A matching adaptive accepted-step
  checkpoint must also pass its existing owner's configuration/state checks.
* A cycle-limited ordinary geomeTRIC trajectory can end in an unaccepted trial.
  With no reliable accepted checkpoint, R1 retains it only as ``observed_frame``
  and emits an incomplete-artifact diagnostic. No automatic restart is permitted.
* ``FAILED_RECOVERABLE`` identifies a failure class that might be recoverable;
  it grants no action permission. ``FAILED_FINAL`` and ``INVALID_RESULT`` retain
  evidence and reasons. ``RECOVERED_CONVERGED`` is reserved and never emitted by
  R1.

The legacy adapter reports bool/string status evidence without treating a legacy
``True`` or ``CycleExceeded`` result as independent proof. Explicit legacy
rejection or cycle-limit evidence cannot be promoted by contradictory success
text in native output. ``analyze_accepted_steps``
in ``pyar.recovery.trajectory`` excludes rejected and unknown steps, and reports
only measured accepted energy/gradient trends. It does not choose remedies.

Pilot support and limits
------------------------

The output protocols tested in this stage are native ORCA **6.1.1**, native xTB
**6.7.1**, and geomeTRIC **1.1.1 with xTB**. Other output versions and diagnostic
routes fail closed. Existing backend execution capabilities are unchanged;
``optimization_diagnostic_capabilities`` exposes only parsed-convergence and
frame-validation capabilities for these pilot routes. It does not advertise
compatible orbital/Hessian restart or automatic recovery.

ORCA normal termination alone is insufficient. The parser requires actual
optimization convergence, finite final energy, ordered final geometry matching
the log and original gradient thresholds. Missing energy is invalid; an explicitly
reported finite zero is valid. Failure/limit evidence takes precedence over
success text.

xTB inspection requires genuine convergence and normal-termination text, actual
Hamiltonian agreement, finite frame-associated energy, and the requested gradient
norm threshold. SCC failure markers take precedence over a success marker. The
existing native adapter consumes ``.xtboptok`` on success, so an absent marker
alone cannot invalidate an otherwise verified fresh successful execution.

geomeTRIC inspection requires matching state/XYZ/optimizer energies and positions,
finite force arrays, consistency between stored forces and logged gradients,
objective identity, genuine convergence and original gradient thresholds. Physical
energy is never replaced by biased/wall energy. Arbitrary fallback XYZ files are
not accepted. Unknown acceptance of intermediate steps is preserved explicitly.

Constrained objectives and non-minimum targets are diagnose-only/unsupported by
this pilot and are rejected before execution. Native adapters cannot silently
ignore requested AFIR/wall controls through this API. Other workflows, Sella,
NEB, scan and IRC retain their current independent paths until their own R4
integration. No new external program, heavy optional dependency, convergence
loosening or scientific fallback is introduced.

Artifact refs in diagnostic receipts are deterministic derived observations, not
new scientific commit IDs. Existing layout-v2 stage/selection commits and writer
ownership remain authoritative. This stage adds no persisted recovery state.

Parser references
-----------------

Marker and threshold checks were verified against the installed pilot outputs and
primary documentation:

* `ORCA geometry optimization manual <https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html>`_
* `xTB geometry optimization documentation <https://xtb-docs.readthedocs.io/en/latest/optimization.html>`_
* `geomeTRIC convergence documentation <https://geometric.readthedocs.io/en/latest/how-it-works.html>`_

Live pilot calculations validate these three execution/diagnosis paths only;
fixtures and analytic tests are not an external chemistry recovery benchmark.
Measured recovery rates and automatic corrective actions belong to later stages.

R1 validation record
--------------------

R1 is based on the reviewed R0 commit ``ecaf7d5`` and changes only the opt-in
service, typed/read-only diagnostic modules, minimal process observation hooks,
observational capability query, tests and this guide. No legacy R0 hazard is
silently declared fixed.

Validation used Python 3.14.7 with the provisioned test environment:

* New R1 contracts: **70 passed**.
* Full ``python -m pytest -q -rs``: **1543 passed, 1 skipped, 17 xfailed**,
  691 warnings and 194 subtests passed. The skip expects a built wheel in
  ``dist/``; package build outputs were kept in a temporary directory.
* Strict Sphinx and wheel/source builds passed. Archive inspection passed the
  existing distribution-content checker. The installed wheel was imported from
  a separate environment outside the source tree and passed typed, readonly
  historical-failure inspection. That environment reused provisioned dependencies;
  it is not evidence of a fresh dependency resolver installation.
* The installed delivery wheel executed each live pilot once, on neutral singlet
  H2 starting at a 1.0 Angstrom bond: xTB GFN2 native and geomeTRIC/xTB at normal
  thresholds; ORCA HF/STO-3G native at normal thresholds. All independently passed
  numerical validation and retained the caller's input geometry. The observed
  versions were xTB 6.7.1, geomeTRIC 1.1.1 and ORCA 6.1.1. Final physical energies
  were respectively -0.982686062615, -0.982686171118 and -1.118735189131 Hartree.
  The two Hamiltonians are different; their energies are not a performance or
  accuracy comparison.
* Fixtures cover contradictory success/failure evidence, current versus stale
  artifacts, accepted versus rejected trials, original versus rounded force
  thresholds, one atom, nonfinite/missing energies and vectors, request/asset
  identity, electronic-state conflicts, absent optional imports, original error
  chains, concurrent writer refusal and zero-write inspection/preflight paths.

These runs establish pilot execution/diagnosis and compatibility, not verified
recovery rates, scientific minima by Hessian, or all-workflow coverage. Other
external backends are not run in this stage. **R2 remains unstarted** and requires
a separate review of accepted-frame policy, action applicability and cost bounds.
