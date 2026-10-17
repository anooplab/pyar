Standalone optimization recovery pilot
========================================

R2 adds explicit, bounded cycle continuation to ``pyar optimize``. R3 adds an
opt-in guided optimizer restart for one evidence-backed geomeTRIC plateau. Recovery is
**off by default**: existing callers, workflows and legacy continuation keep
their existing behavior. This pilot neither enables recovery in other workflows
nor establishes Hessian, chemical identity, TS, IRC or pathway acceptance.

Supported execution
-------------------

Only unbiased, unconstrained, closed-shell singlet GFN2-xTB minima through
geomeTRIC are enabled. The execution requires geomeTRIC **1.1.1** and verifies
xTB **6.7.1**, its actual Hamiltonian and successful electronic evaluation on
every backend call. Unknown versions and unsupported routes stop safely.
Native ORCA/xTB continuation is unsupported because R1 has no verified accepted
incomplete-frame contract for those optimizers. Orbital/SCF restart is also
unsupported: electronic state compatibility has not been established. No
orbital file is copied or adopted. These exclusions are capability limits,
not recovered failures.

Example::

    pyar optimize input.xyz --backend xtb --geometry-optimizer geometric \
        --opt-cycles 20 --recovery safe --recovery-max-attempts 3 \
        --recovery-max-cycles 140 --recovery-max-evaluations 200 \
        --recovery-wall-seconds 600 --check

Remove ``--check`` to execute. Preflight validates the entire input batch before
starting any input. ``--check`` creates no directories, locks, reports or
checkpoints. It does not launch a backend or promise that its installed version
will pass the runtime checks. Each input has its own bounds, rather than sharing
a hidden batch budget. Unsupported recovery requests fail preflight with
``recovery_unsupported``.

SAFE and GUIDED policy
----------------------

The first attempt uses the requested cycle count. A cycle limit permits a new
attempt only when at least two explicitly accepted observations show both an
energy reduction above 1e-9 Hartree and an RMS gradient reduction above 1 percent.
These are conservative action predicates, not relaxed convergence thresholds.
The next attempt starts at the last accepted geometry, with up to twice the
previous cycle allocation and a freshly initialized ordinary optimizer. Requested
Hamiltonian, charge, spin, objective, atom ordering and thresholds remain fixed.
A limited final trial is never adopted: geomeTRIC checks its iteration limit
before deciding whether that trial was accepted.

``--recovery safe`` retains the R2 cycle-extension rule and stops on stagnation.
``--recovery guided`` adds one optimizer action after four accepted frames show
an energy span no greater than 1e-7 Hartree while the RMS gradient remains above
the requested threshold and varies by no more than five percent. It rolls back
one accepted step, starts a fresh geomeTRIC 1.1.1 process (therefore rebuilding
the approximate Hessian), and sets a smaller 0.05 Bohr initial trust radius with
a 0.15 Bohr maximum. The action is versioned in its immutable attempt intent.
The next attempt keeps the original GFN2 Hamiltonian, charge, spin, objective,
atom ordering and convergence thresholds. This is an optimizer restart, not a
claim that the plateau is a minimum.

Neither policy changes coordinates through collision repair, changes methods,
lowers thresholds, adds restraints, switches optimizer families or repeats an
interrupted process. SCF failure remains a stop condition: current evidence
does not prove electronic-root continuity for an automatic solver change.
Guided collision repair and preoptimization remain disabled because this pilot
has no workflow-owned chemical-identity authorization. An incomplete usable
result remains provisional and makes the command fail; only independent
numerical convergence publishes ``result_<name>.xyz`` and counts as optimized.
A successful retry is labeled ``RECOVERED_CONVERGED``.

The ``pyar.recovery.advanced.finite_difference_gradient`` helper provides a
bounded central-difference comparison for callers with an energy evaluator and
an analytic Cartesian gradient. It charges two energy evaluations per Cartesian
component, rejects nonfinite data and refuses work above its explicit cap. It
reports consistency evidence only; it never changes a threshold or authorizes
convergence. The helper is not invoked automatically by the standalone pilot.

Bounds are explicit and persisted:

* ``--recovery-max-attempts``: default 3, maximum 8, including the initial attempt.
* ``--recovery-max-cycles``: default 700, maximum 10000; allocations are reserved
  before launch and clipped to the remaining total.
* ``--recovery-max-evaluations``: default 1000, maximum 10000; the calculator
  durably counts physical calls before executing them and rejects calls beyond
  the remaining allowance, including rejected optimizer trials.
* ``--recovery-wall-seconds``: default 3600, maximum 86400; an absolute persisted
  deadline covers the run. The parent kills the complete process group on timeout.

``backend_evaluations`` reports charged evaluations. For sealed attempts it is
the counter of reserved physical calls; if an attempt loses its receipt, the entire remaining
allowance is conservatively charged. It is not then an observed call count.
Cycle allocations are upper bounds, not measured optimizer iterations or CPU
seconds. A backend evaluation can contain many electronic iterations. No recovery
rate or performance advantage follows from these accounting bounds.

Persistence and explicit resume
-------------------------------

New opt-in runs use the established ``job_<name>`` root and sibling ``RunLock``::

    job_<name>/
      recovery/
        run.json                 # immutable request identity, bounds, deadline
        request.json             # immutable original request
        attempts/attempt_001/
          input.json             # immutable intent, parent and action rationale
          job_<name>/            # existing adapter's native files, untouched
          worker.out
          result.json            # execution receipt and observations
          commit.json            # immutable validated attempt publication
        attempts/attempt_002/...
        outcome.json             # derived summary; never restart authority
      result_<name>.xyz           # derived, only after validated convergence

The worker process uses an explicit subprocess working directory. The controller
does not add process-wide ``chdir`` calls. Inside each native job, a version-pinned
standard geomeTRIC driver records accepted steps as immutable
``ore_accepted_NNNNNN.json`` files and an index. It observes the initial evaluated
frame and subsequent accepted steps; it does not change the optimizer's
objective, thresholds, coordinate setup or acceptance decision. Original
inputs, raw backend output, force/energy units, exact settings and SHA-256 artifact
references remain available. These are standalone attempt records, not layout-v2
workflow/stage commits, and do not replace workflow-owned persistence.

Resume requires the same original input, scientific settings and policy::

    pyar optimize input.xyz --backend xtb --geometry-optimizer geometric \
        --opt-cycles 20 --recovery safe --recovery-max-attempts 3 \
        --recovery-max-cycles 140 --recovery-max-evaluations 200 \
        --recovery-wall-seconds 600 --recovery-resume

The deadline is not reset. Every native artifact and accepted checkpoint is
revalidated against its recorded hash and exact intent. A request's original
scientific fingerprint retains its initial geometry; each attempt additionally
records its starting geometry and exact request fingerprint. The objective
fingerprint excludes starting geometry and numerical allocation, enabling an
explicit check that continuation changed only the permitted controls.

A crash after receipt publication can be recovered by validating that receipt
and sealing the attempt; a crash after commit can reconstruct the derived
outcome without another backend call. A started attempt without a valid receipt
stops as interrupted, rather than repeating uncertain work. Corrupt records,
symlink paths, legacy caches and mismatched requests are rejected. Stale lock
metadata requires the separate explicit ``--recovery-recover-lock`` option,
which still requires acquiring the operating-system lock. No active writer can
be displaced. Derived reports and final XYZ files never authorize a retry.

Review and validation
---------------------

R2 starts from merged R1, ``b77e59a``. The preserved uncommitted prototype in the
primary checkout is excluded. Focused fixtures cover recovery then convergence,
finite repeated failures, disabled recovery, progress rejection, corrupt and
incompatible checkpoints, receipt/commit crash boundaries, partial-worker replay,
cache rejection, objective identity, read-only preflight and whole-process-group
timeout. A synthetic executable exercises the real worker and geomeTRIC driver
without a chemistry executable; it is not a chemistry benchmark.

A real local neutral singlet H2 smoke used xTB 6.7.1 / geomeTRIC 1.1.1 with GFN2,
normal requested thresholds and an initial three-cycle allowance. The first
attempt reached its cycle limit. Continuation from the accepted checkpoint
converged independently. This is a small protocol check, not a measured recovery
rate. R3 adds deterministic synthetic evidence for a guided plateau restart,
unchanged objective identity, a new initial trust radius and noisy-gradient
rejection. No real chemistry run has exercised the guided action.

R2 validation used 48 focused SAFE cases and the full suite:
**1591 passed, 1 skipped, 17 expected legacy failures, 694 warnings and 194
subtests passed**. The skip requires a wheel in the checkout's ``dist/``;
built distributions were instead inspected separately in a temporary directory.
Strict Sphinx and source/wheel builds passed. The installed wheel was tested in
a separate environment sharing provisioned dependencies; this does not claim
fresh dependency resolution. Its real recovery used two attempts, nine reserved
cycles and nine physical calls, with energy -0.982686174565 Hartree. The ordinary
run at the same requested objective gave -0.982686171118 Hartree. The difference
is about 3.45e-9 Hartree. A two-call hard-cap run stopped without a completed XYZ.
These observations establish a small smoke protocol, not broader chemistry
coverage or performance statistics. R3 currently enables only the xTB 6.7.1 /
geomeTRIC 1.1.1 minimum route. It does not add an ORCA SCF ladder because the
available output contract cannot establish electronic-root continuity; it also
does not add optimizer-family switching or collision repair. Those actions
remain blocked on verified backend controls and workflow-owned identity gates.

R3 branch validation ran 57 focused recovery cases and the complete suite:
**1598 passed, 3 skipped, 17 expected xfails, 694 warnings and 194 subtests
passed**. The suite used a temporary editable install so the installed-entry-point
test read this checkout's metadata. Strict Sphinx and source/wheel builds passed.
These are software and synthetic-fixture checks; the new guided action has not
been exercised with a real chemistry run.
