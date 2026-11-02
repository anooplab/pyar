Architecture and ownership
==========================

The sole command entry is ``pyar.cli:main``. It loads task adapters lazily and
observes request preparation/execution/result through ``cli_runtime``. Parsers
translate user intent; request/service owners validate scientific settings;
workflows sequence sampling, execution, selection and persistence; backend
adapters translate physical values and interpret native evidence.

Shared electronic-state/preflight checks live in ``cli_validation`` independently
of the retired parser. ``cli_inputs`` builds task inputs. ``optimiser.optimise``
copies settings and injects the current molecule's charge/multiplicity/SCF type
and xTB unpaired-electron count before fingerprinting/backend construction.
Every backend attempt rewrites its trial XYZ from that request. Result-cache
reuse requires validated input/model/artifact provenance.

Scan orchestration owns route validation, chemistry preparation, immutable
request/recovery, orientation execution, result publication and finalization.
NEB separates ordered request validation, physical settings, cache adoption,
stage-specific execution, publication and failure receipts. Stage handlers
return result/input/artifact/dependency evidence to a single publication owner.
Recursive dependency/hash validation precedes acceptance of upstream evidence.

Current persistence owns supported series/layout identity, lock ownership,
immutable artifact references and workflow checkpoints. Derived reports are not
restart authority. Atomic mutable publication is distinct from an immutable
stage commit; neither creates a multi-file transaction by itself. Workflow and optimization entry points pass explicit run directories; production
code does not mutate the process working directory. Native command helpers use
``cwd=`` and resolve redirected streams within that directory. Relative-file
Python calculators and geomeTRIC execute in a child process whose directory is
set at creation. Recovery captures the actual native process observations, not
just the Python worker exit status.

Keep scientific kernels and owner acceptance rules explicit. A helper is useful
when it clarifies responsibility and dependencies. Moving branches between
functions alone is not a reduction in aggregate scientific decisions.

CLI invocation ownership
------------------------

``cli.invoke`` preserves the argument-based dispatch API and sequences console
controls, task configuration, observed execution and response publication.
Profiles and command arguments still use ``run_config`` resolution. Task
adapters remain the owners of scientific requests and preflight checks.

Each invocation owns its prepared request and a one-time configuration
snapshot. Repeated preparation callbacks can update the request recorded in run
provenance without overwriting that first saved configuration. Nested calls
have separate invocation state and context-local execution observations.

Execution captures workflow failures, including interrupts, so a started run
can attempt provenance publication before the failure is re-raised. A handled
provenance-write failure preserves an existing workflow failure; if the workflow
succeeded, it causes an invocation failure while retaining the scientific
result. Response rendering follows stream/logger cleanup. Check and dry-run
calls retain their read-only workflow behavior and existing output routing.

The scan-bond adapter resolves input-mode defaults and backend options before
continuation checks. It reuses ``scan_request`` input/capability validation and
the workflow's signed-request and restart checks. Those checks finish before
``prepared_run`` publishes the request or execution starts. A private prepared
value carries the same invocation mapping and preflight evidence into reporting
and execution; scientific settings remain owned by the existing services.
Scan results are published before a structured failure produces its existing
nonzero exit status. Lock errors retain their CLI diagnostic boundary.

Reaction request resolution separates calculation controls from loaded
reactants. Bias schedules use the existing controller and gamma-schedule rules;
fragment and combined states use the shared electronic-state validators. The
same backend/method restrictions run before the existing reaction workflow
builds its canonical restart request. Private stage values carry resolved data;
``ReactionRequest`` remains the supported mutable envelope, and preflight and
restart verification retain their existing owners. Resolution itself performs
no backend calculations or state publication.

Microsolvation request construction validates each molecule's geometry and
state before moving to the next molecule. Sampling/surface checks precede the
existing shallow backend-settings copy; objective/wall, placement, resampling
and site checks follow in their established order. Every cumulative cluster
state is checked with the shared multiplicity-combination rule and backend
capabilities. These validation helpers operate on the existing frozen request;
the public fields, defaults and serialized request remain unchanged.

Scientific constants and backend errors
---------------------------------------

``pyar.constants`` owns radii, connectivity factors and unit conversion values.
Different radius models and historical conversion precisions retain distinct
names and exact values. Existing data modules re-export their original tables;
centralization does not authorize replacing scientific profiles or thresholds.

Backend ``optimize`` and energy/gradient ``evaluate`` operations raise
``BackendExecutionError`` subclasses for failures: ``BackendSCFError``,
``BackendGradientError``, ``BackendGeometryError``, ``BackendOutputError`` and
``BackendConvergenceError``. Unexpected adapter exceptions retain their cause.
A convergence error's ``usable_geometry`` flag records an incomplete-frame
observation; it cannot certify convergence or chemistry. The optimization owner
validates that frame before returning ``OptimizationOutcome.CYCLE_EXCEEDED``.
The release-candidate enum likewise represents evidence for further validation.
Their serialized values and current workflow result schemas remain unchanged.
Internal predicates, optional metadata and low-level process return codes retain
their ordinary contracts; they are not scientific success flags.

Directory execution uses private one-shot IPC in an automatically removed,
mode-0700 temporary directory. It never reads public pickle inputs or stores IPC
as restart/cache authority. It returns result snapshots, typed exceptions,
backend logs, native-process observations and sequential Python/NumPy RNG state.
Model/calculator objects remain in the worker. Workers import the same source or
installed package as their caller. Existing relative artifact references,
scientific fingerprints and acceptance rules remain owned by their workflows.
This removes shared CWD mutation; it does not make mutable workflow objects,
registries or shared random streams safe for concurrent use. Process startup and
IPC add overhead per optimization, scan or NEB invocation.

Mocked in-process scientific tests explicitly use ``inline_directory_worker``
with test-only workspace adapters. Separate tests execute real child processes,
parallel isolated workspaces and native xTB; mock fixtures do not establish
external-backend scientific validity.

Growth settings and stage ownership
-----------------------------------

Aggregation and fixed-addend growth construct a ``GrowthPolicy`` once for an
addition. ``SelectionOptions`` owns its comparison/selection controls. The
policy copies and freezes the outer backend-settings mapping; attempts receive
ordinary parameter dictionaries with the existing loose/normal thresholds.
``GrowthRun`` and ``SeedRun`` identify filesystem owners independently of those
settings. Preparation, the unchanged ten-round retry budget, resumed selection,
refinement and final composition-group publication have explicit stage outputs.
``GrowthStopped`` belongs to the workflow owner; it never becomes a candidate
pool or authorizes publication of a partial step.
