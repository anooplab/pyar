# Optimization Recovery Engine: R0 audit

Date: 2026-10-17. Status: **R0 audit and characterization only; R1 unstarted.**
No production optimizer, scientific default, CLI option, cache rule, or persisted
workflow format changes in this audit. Expected failures below are reproducible
safety gaps, not recovery successes. This report proposes decisions for review.

## Baseline and scope

- Audited commit: `d2c38ec280ec01fec5302b8b5dacea038c054ea4` (master,
  PR #61 merge). GitHub `commits/master` returned this same SHA during R0.
  All returned check runs were completed/success, including Python 3.10–3.13,
  installation extras, Sella, documentation and distribution installation checks.
- Package metadata and installed distribution: `pyar-chem 1.5.1`. Provisioned
  interpreter: Python 3.14.7, NumPy 2.5.3, ASE 3.29.0, geomeTRIC 1.1.1,
  Sella 2.6.0, pytest 9.1.1. Python 3.14 local results do not replace supported
  Python CI coverage.
- A separate `codex/ore-r0-audit-20261017` worktree starts at that commit.
  The original checkout's uncommitted retry prototype remains untouched. Its
  modified files are `docs/{architecture,modern_cli}.rst`,
  `pyar/scripts/optimize.py`, `tests/test_modern_optimize.py`, and the untracked
  `pyar/optimization_recovery.py` / `tests/test_optimization_recovery.py`.
  Those changes are **not included** in R0 or its test baseline.
- Confirmed the imported package comes from the R0 worktree, rather than the
  editable installation's original checkout. Baseline full suite before adding
  fixtures: **1461 passed, 1 skipped, 691 warnings, 194 subtests passed**.
- xTB 6.7.1, ORCA and OpenBabel executables are discoverable. Gaussian, Psi4,
  MOPAC and Turbomole `define` are not discoverable. No new real external
  chemistry run was performed for R0. Fake-executable tests are not chemistry
  benchmarks. ORCA/OpenBabel runtime versions were not probed by this audit.

Read sources: repository `AGENTS.md`; the private ORE plan dated 2026-10-17;
master remaining-work plan dated 2026-10-15; run-layout-v2 plan; microsolvation
upgrade plan dated 2026-10-17; and the repository's
[run-layout contract](run_layout_v2/contract.md). Private plans were read from
Dropbox `/PyAR-2.0/doc/`; they are not copied into this public audit.

Compatibility reconciliation:

- Master M11/M12 requires measured optimizer evidence, preservation of failures
  and finite cost; it does not authorize changing optimizer defaults. Existing
  M9/M10 association references, stationarity/Hessian, IRC and energetic gates
  must remain intact (`pyar/workflows/{association_reference,reaction_report}.py`).
- Microsolvation MS0–MS4 is marked complete in its closeout and is present in
  source: contacts, optional walls, automatic shell policy, and budgeted whole-
  cluster resampling already exist. ORE must preserve these owners and must not
  reimplement them or claim improved basin recovery. MS5/MS6 remain separate.
- Layout v2 is implemented. Its immutable commits, artifact hashes, native
  executor subtrees and sibling writer locks already own restart authority.
  A recovery report cannot replace any of these. Legacy v1 remains readable
  under its original rules, without retroactive guarantees.

## Entry-point and call graph coverage

[Static inventory](optimization_recovery_inventory.tsv) lists source, line,
owner, and call/definition/status kind for optimizer boundaries and relevant
status helpers. It was obtained by AST traversal of all 228 package Python
files and the repository benchmark Python scripts, then checked against text
search for optimize/optimise, optimizer run calls, native subprocess execution,
RDKit minimizers, IDPP, geodesic and geometry preparation. It records 166 rows.
It includes definitions not routed by current public dispatch; those are
identified below. It is a snapshot, not a runtime tracer or proof of paths in
external packages. Source references below refer to the audited commit.

| User/API entry | Actual execution chain and source anchors | Ownership / separate bypass |
|---|---|---|
| `pyar optimize` | `modern_cli.py:49` → `scripts/optimize.py:35,72` → `optimiser.py:129,195` → native `optimize()` or geomeTRIC subprocess | Standalone per-input native `job_<name>`; modern loop catches backend/file-not-found per input |
| `pyar-optimiser` | `scripts/optimiser.py:24,88` → `optimiser.bulk_optimize:250` → `optimise` | Legacy bulk list comprehension can abort on an exception; no independent retry owner |
| `pyar-cli` aggregate/react/solvate and legacy wrappers | `cli.py:942,1025`, `scripts/react.py:202`; workflow APIs below | Compatibility aliases in `aggregator`, `reactor`, `solvator`, `_growth` forward to canonical workflows; not additional engines |
| Aggregate | `workflows/aggregate.py:311,496` → `growth/service.add_one` → `optimise:582,690` | v1 pathway directories; v2 `aggregation/layout_v2.py` commits selections, not solver evidence |
| Grow | `workflows/grow.py:39,149` → `add_one` | v2 `growth/layout_v2.py`; ten-round incomplete pool and loose→normal refinement belong to growth |
| Legacy solvation | `workflows/solvation.py:155` → `add_one` | Existing `solvation/state.json` and aggregate directories; no layout-v2 migration here |
| Conformer generation | `workflows/conformer.py:270,276,360,372` → RDKit MMFF/UFF multi/single conformer minimizers | **Bypasses generic wrapper**. Force-field status and kcal/mol energy are separate from later QC results |
| Conformer backend refinement | modern/legacy conformer scripts → `conformer_search` → `_refine_with_backend:887,896` → `optimise` | v2 `conformer/layout_v2.py`, `refinement/execution`; a failed call can abort refinement |
| Microsolvation | `scripts/modern_microsolvate.py:152` → `workflows/microsolvation.py` → `_optimize_candidates:149` → `optimise:168,200,218` | Unbiased screen, optional wall search, unbiased re-relaxation are distinct physical stages; per-candidate exceptions collected |
| Microsolvation resampling | `_resample_cluster_pool:571` → `_optimize_candidates` under persistent physical-evaluation budget | Resampling proposal/RNG/attempt ownership already exists; interrupted proposals are not replayed blindly |
| React/AFIR | `workflows/reaction.optimize_all:780,840` → `optimise` → geomeTRIC `PyarGeometricCalculator` + AFIR/softmin | Reaction cycles, gamma schedule, trace and release policy stay workflow-owned |
| AFIR free release / adaptive re-probe | `relax_without_afir_bias:271,292`, `_adaptive_release_probe:460,511` → `optimise` | Native free objective differs from biased search; bounded margin retries already exist and cannot be replaced by ORE |
| Scan-bond native ORCA | `workflows/scan_bond.py:490` → `backends/orca_scan.run_orca_bond_scan:176` | **Bypasses generic wrapper**; scan-native input/output/profile/trajectory; scan is not TS proof |
| Generic constrained scan | same dispatch → `backends/bond_scan._run_bond_scan:83,145` → geomeTRIC `Optimizer.optimizeGeometry` | **Bypass** with per-point distance constraint, finite result and optimizer-state gates |
| Free relaxation after scan | `workflows/scan_bond.py:551` → `optimise`, native ORCA or geomeTRIC | Separate unconstrained relaxation and product identity; does not change committed scan |
| Scan/react characterization | `workflows/scan_path.py:120`, `reaction_characterization.py:396` → `scan_path.run_path_stages` → `neb.run_neb` | Separate cumulative path executor; `PathLayout` wraps immutable stage commits |
| Standalone NEB/TS/IRC | `scripts/modern_path_stage.py:226,230,238,242,246`, `pyar-neb` → `neb.run_neb:1554` → `_execute_stage:1010` | **Bypass**; standalone roots retain their current format |
| Endpoint preparation, geomeTRIC TS, IRC forward/backward, endpoint relax | `neb.py:801,904,1081,1198,1243,1286` → in-process geomeTRIC `Optimizer` | **Bypass**; fresh scratch, actual optimizer-state bool, retained path; numerical convergence is not TS/connection proof |
| Sella TS | `neb.py:860,875,877` → `Sella(order=1, internal=...)` | **Bypass**; preserves exception chain, validates finite trajectory, returns convergence bool and distribution version; no minimum-optimizer substitution |
| NEB band | `neb.py:1094,1131` → geomeTRIC `ElasticBand` / `OptimizeChain` | **Bypass**; fixed endpoints, band gradients, climbing-image/interior maximum gates |
| Separated association reactants | `workflows/association_reference.py:27,91` → `_optimize_geometry`, then `_frequency` | **Bypass**, isolated fragment electronic states, identity and minimum checks; energy reference is sum of fragments |
| Path initialization | `neb.py:206,217` ASE `idpp_interpolate`, `273,295` geodesic `smooth` | Auxiliary pair-distance path objectives, not molecular minimum convergence; preserve endpoints and waypoint |
| Native ORCA OptTS benchmark | `benchmarks/orca_optts.py:151,223` → direct subprocess, then `run_neb` validation | **Bypass**; own fresh-output rejection, parser, atom-order, TS/frequency/IRC checks |
| TS/conformer benchmarks | `benchmarks/ts_optimizer.py:585`, `benchmarks/conformer.py:626` | Delegate to existing NEB/conformer owners; no new production optimizer |
| AIMNet2 optimizer executable | `AIMNet2/calculators/aimnet2_ase_opt.py:10,19,99,137` → ASE LBFGS | **Bypass**; writes results after `run` without checking returned convergence bool |
| ASE/MLatom and dimer helper | `backends/ASE.py:14,34,40,60` → selected ASE optimizer / `MinModeTranslate` | **Bypass**, trajectory returned, convergence return discarded; helper not selected by `build_geometry` |
| ANI direct API/module | `backends/ani.py:107,112,146` → ASE BFGS | **Bypass**, returns final energy, no typed convergence; deliberately excluded from modern bulk route |
| MLopt direct executable | `backends/mlopt.py:8,29`, `mlatom/__init__.py:23` → external `mlatom.optimize_geometry` | **Bypass** / external delegated optimizer; no in-package accepted checkpoint contract |
| Legacy MLatom direct class | `backends/mlatom_aiqm1.py:53,62` → external MLatom + Gaussian log | Alias normalizes to `aiqm1_mlatom` before dispatch, so the old `mlatom_aiqm1` branch is unreachable through normal `build_geometry`; direct import remains possible |
| Legacy ORCA-AFIR file | `backends/orca-afir.py:54` | Script with hyphenated filename; not dispatched by registry; runs `engrad` yet uses termination as success. Record as legacy direct surface, not supported optimization |
| Placement orientation minimization | `scan.py:9,18` SciPy differential evolution (`global_opt`) | Placement/angle search, not electronic relaxation or an ORE converged result |
| Scientific clustering reference optimizations | `scripts/scientific_clustering_benchmark.py:97,140` → direct xTB subprocess | **Bypass**; fresh attempt directories, executable/input/settings/output/log hashes, anchored normal termination, finite energy, final gradient norm ≤1e-3 and ordered symbols. Distinct benchmark policy, not generic cache |
| Dataset / pilot scripts | `benchmarks/{clustering_distances,structure_comparison}/generate_dataset.py` RDKit; `benchmarks/microsolvation/ms4_perturbation_pilot.py:220` → `optimise` | Research-only dataset preparation/pilot, status use included in inventory; no authority to change runtime defaults |

### Native backend routing and evidence

`optimiser.build_geometry:62–126` dispatches geomeTRIC first when selected, else
native adapters. `backend_capabilities.py:11–196` describes calculator/native/
charge/spin/options support but has **no parsed-convergence, accepted-checkpoint,
or recovery capabilities**. Do not infer safe recovery from `native_optimization`.

| Adapter / source | Present success evidence | Native artifacts / deficiency |
|---|---|---|
| ORCA `backends/orca.py:78–128` | exit 0 + `ORCA TERMINATED NORMALLY` → True | `trial_N.{inp,out,xyz,gbw}`, `result_N.xyz`; no geometry-convergence check; absent energy becomes 0; raw loadtxt lacks symbol/count/finite/freshness checks |
| xTB `backends/xtb.py:57–142` | exit 0; unknown-option text rejected; `.xtboptok` checked **before** SCC failure markers | `trial_N.xyz`, `xtb.out`, `energy`, `xtbopt.xyz`, `xtbopt.log`→`traj_N.xyz`; stale success/native files not invalidated before invocation; no finite/shape validation |
| geomeTRIC `backends/geometric.py:634–715` | nonzero exit except exact iteration-limit text is error; reads XYZ and physical energy; limit returns `CycleExceeded`, otherwise True | `geometric.out`, `trial_N_optim.xyz`/`opt.xyz`/newest other XYZ, `pyar_geometric_state.json`; validates energy presence, not finite value/frame correspondence/convergence evidence. Generic state is latest calculator evaluation, not automatically an accepted checkpoint |
| Adaptive geomeTRIC `backends/adaptive_geometric.py:122–179` | release-candidate / cycle-exceeded / optimizer state handled explicitly | `pyar_bias_controller_state.json` saved on accepted geometry; on limit restore accepted frame, clear calculator cache and re-evaluate; release candidate is workflow event, not necessarily stationary minimum |
| Gaussian `backends/gaussian.py:35–140` | exit 0 + optimization completed + SCF Done + final normal termination | `.com`, `.log`, `.chk`; native input omits Opt keyword, hence modern resolver requires geomeTRIC. Stronger markers than ORCA but no independent finite/shape check |
| Psi4 `backends/psi4.py:35–137` | exit 0 + `Optimization is complete!` | `.in/.out`; native method hardcoded B97-D/def2-SVP rather than arbitrary QC controls; missing energy can leave local unbound, coordinate parser does not require exact count |
| Turbomole `backends/turbomole.py:92–290,665` | native jobex marker / internal statpt convergence; differentiated SCF/gradient/update/limit strings | `coord`, `control`, `energy`, `gradient`, `converged`/failure files, logs; accepted/effective versus physical AFIR energy needs explicit contract |
| xTB/Turbomole `backends/xtb_turbo.py:61–113` | statpt convergence / failures / cycle limit; hardcoded 250 loops | xTB gradients and Turbomole updates; rewrites biased energy/gradient; gamma==0 dispatches directly to standalone xTB, otherwise this hybrid route |
| MOPAC `backends/mopac.py:61–127` | exit 0 + `.arc` exists | `.mop/.arc/.log`; no independent geometry stationarity marker, stale arc risk; heat of formation converted to atomic units is not automatically the same energy convention as an ab initio potential |
| OpenBabel `backends/babel.py:53–116` | minimizer exit status, result file and parsed obenergy | `result_N.xyz`, logs/energy output; presence/finite parsing is not force-residual convergence proof; force-field convention retained |
| AIMNet2 `backends/aimnet_2.py:85–137` | helper process success + output reads | `aimnet2.out`, optimized XYZ; helper ignores LBFGS convergence return, converts output eV→Hartree; require model/version/force/accepted-frame proof |
| AIQM1 MLatom `backends/aiqm1_mlatom.py:47–94` | subprocess success + optimized XYZ/energy reads | `mlatom.out`, minimum XYZ; external optimizer diagnostics not normalized; no SCF remedies inferred for other ML models |
| xTB→AIMNet2 / AIQM1 `backends/xtb_aimnet2.py:63–127`, `xtb_aiqm1.py:43–102` | sequential process/output checks → True | stage-specific optimized files; distinct preparatory/final models; final model completion needs independent validation |
| ORCA-AIQM1 `backends/orca_aiqm1.py:48–84` | normal termination position check → True | `.inp/.out/.xyz`; missing energy becomes zero; has same optimizer-evidence gap as ORCA |

No new capability flags are added in R0. R1 should mark a route diagnose-only
unless its convergence and accepted-frame evidence is independently tested.

## Status taxonomy and consumer compatibility

The current public status vocabulary is bool/strings, not an auditable typed
result. `optimiser.py:28–51` recognizes True/`converged`, and three spellings of
cycle limit; `is_usable` includes both. Python set equality also accepts integer
1 as True; avoid truthiness/equality coercion in a new typed API.

| Observation | Meaning required in ORE | Present mismatch / owner |
|---|---|---|
| Process exited normally | Backend execution ended | Not a convergence or finite-result certificate |
| Optimizer satisfied requested criteria | Numerical convergence at the requested objective | Requires solver evidence plus validated frame/energy/gradient |
| Accepted finite frame, limit reached | `INCOMPLETE_USABLE` with reason `OPT_MAX_STEPS` | May seed authorized exploration, cannot be labeled converged |
| Workflow gate passed | Owner accepts identity/constraints/stationarity/path connection | Independent from solver convergence; cannot be overridden by generic controller |
| Cached/reused file | Provenance source, not a convergence status | Must verify request + accepted frame + convergence; current generic cache returns True |
| Backend/error/invalid artifact | Diagnosed failure or `INVALID_RESULT` | No fabricated energy, silent spin/method changes or guessed remedy |

**Every `is_usable` consumer:** wrapper application at `optimiser.py:197`, bulk
count/filter `267,273`, modern optimize `scripts/optimize.py:76`, conformer
`workflows/conformer.py:901`, microsolvation `workflows/microsolvation.py:177,
201,219`, AFIR retry `workflows/reaction.py:512` and discovery `843`, and research
pilot `benchmarks/microsolvation/ms4_perturbation_pilot.py:229`. Definitions and
`is_success`/`is_cycle_exceeded` call sites are in the static inventory.

Proposed migration decisions (review required; no changes applied):

- Standalone optimize currently counts `CycleExceeded` as optimized/completed.
  The migrated facade must report incomplete separately with non-success final
  status; retain frame as evidence. Preserve the legacy callable until approved.
- Conformer refinement and final microsolvation ranking currently accept usable
  incomplete results. Preserve old records/readers; migrated stages must expose
  provisional counts and require genuine convergence for final minimum claims.
  Wall-search may pass an incomplete accepted frame to its explicitly owned
  unbiased stage, but final shell/site/composition/energy gates still apply.
- Growth `service.py:578–611` intentionally retains cycle-exceeded structures
  for up to ten rounds, then selects successful ones; refinement `690–717`
  accepts `is_success`, can fall back to the prior loose set. Preserve thresholds
  and selection behavior during diagnosis; a cache that turns provisional into
  True can defeat this distinction and requires separately reviewed correction.
- AFIR discovery and adaptive biased retries explicitly need provisional frames.
  Free release `reaction.py:292,295,321,481` requires `is_success`, plus product/
  survival evidence. A released product is not a TS; biased objectives never
  become physical barriers. Do not replace existing margin/release policy.
- Scan free relaxation requires `is_success` (`scan_bond.py:563`) and the result
  artifact; constrained scan validity is separate from free product identity.
- Direct bool checks in xTB hybrids, Turbomole, helper executables and NEB/Sella
  belong to their specific adapter contracts (table above). Sella's `bool(run)`
  is only numerical saddle-optimizer status. NEB `_STAGE_GATES:54–64` and
  `_stage_gate_passed:126`, TS frequency index/stationarity, both IRC branches,
  endpoint minima and connection remain required. Association reference must
  still validate each isolated fragment before summing physical energies.

Stable proposed failure classes: optimizer limits/stagnation/divergence/internal
coordinates/Hessian; SCF limits/oscillation/electronic-state mismatch; invalid
energy/gradient/frame; constraints/identity/path scientific rejection; unsupported
capability/input; resource/time/interruption; incompatible checkpoint/incomplete
artifact; UNKNOWN. Parsing must retain evidence and uncertainty. A generic
BackendExecutionError or log substring alone cannot safely choose a remedy.

## Artifacts, restart authority and path ownership

| Surface | Producer / consumers | Restart/cache authority and risk |
|---|---|---|
| `job_N/trial_N.xyz`, `result_N.xyz` | `backends.SF:60–78`, adapters, wrapper, growth/conformer/reaction | Input only written if absent; generic wrapper reads existing result before backend construction. No request, ordered-symbol, convergence, or method signature is required |
| geomeTRIC state / output trajectory | calculator writes state on evaluations (`geometric.py:312–343,431`); wrapper reads physical energy and XYZ | Direct wrapper lacks frame comparison. Generic cached route checks finite physical energy, shape and coordinate tolerance but not input identity, method, bias/constraint settings or convergence status |
| Adaptive controller checkpoint | `accept_geometry:275`, `_save_controller_checkpoint:261`, `_checkpoint_configuration:227`; adaptive runner restores | Configuration compatibility already checked; accepted-step state is distinct from latest trial calculator state. Generic ORE must not assume all state.json files mean accepted geometry |
| Native ORCA / xTB checkpoints | `.gbw`, `xtbrestart`, native solver outputs | Existing code does not provide a validated electronic identity/checkpoint compatibility certificate. Observe files, do not automatically resume orbitals |
| Generic scan | `bond_scan.py:24–78,83–180` request + summary + artifact hashes | Checks physics/request, trajectory completeness, finite energies, target distance and optimizer state. Do not replace with wrapper cache policy |
| ORCA scan | `orca_scan.py:114–198`, workflow scan request `scan_bond.py:97–153` | Native recovery checks termination, all symbols/frames, full profile/grid/distances; does not independently prove all unconstrained DOFs reached tolerance at each point. Keep this uncertainty explicit |
| Path stages | `neb.py:660–779`, `PathLayout` | Native summaries hash inputs/artifacts/upstream summaries and physical settings; stage parameters and scientific dependencies checked. Failed/running summaries are not success. v2 immutable commits are owner authority |
| Reaction discovery products | `state/reaction_v2.py:148–172,200,238–270` | Cycle/discovery commits and bounded product source snapshots; characterization hashes immutable discovery commit rather than mutable reports |
| Selected pools / final copies | growth/aggregation/microsolvation/conformer layout owners | Selection commits/refs own identity. Manifests and presentation copies are derived, not standalone native restart certificates |

Concrete fresh-v2 native execution locations already in code:

- Grow: `steps/step_NNN/attempts/attempt_000/execution/` (`workflows/grow.py:140`).
- Aggregate: `pathways/pathway_NNN/steps/step_NNN/attempts/attempt_000/execution/`
  plus legacy native subtree (`workflows/aggregate.py:493`).
- Microsolvation: `steps/step_NNN/attempts/attempt_000/execution/candidates/`,
  optional `wall_search/` and `unbiased_relaxation/`; resampling separately owns
  `resampling/attempts/attempt_NNNNN/` and its persistent evaluation budget.
- Reaction discovery: `discovery/cycles/cycle_NNN/attempts/attempt_000/execution/`
  containing native `gamma_*`, orientation, `job_*`, trace and release evidence
  (`state/reaction_v2.py:148–159`).
- Scan: `orientation_NNN/attempts/attempt_000/execution/{scan,relax}/`
  (`workflows/scan_bond.py:422–427`).
- Scan/react characterization: pathway `attempts/attempt_NNN/execution/`;
  commits in `stages/STAGE/attempts/attempt_NNN/commit.json`
  (`workflows/path_layout.py:33–80`). Changed path request selects another attempt.
- Conformer: `refinement/execution/job_N/` (`workflows/conformer.py:889`).
- Standalone optimize/stage APIs and legacy v1 use existing native roots.

### ADR-ORE-001: proposed ownership boundary

**Decision proposed for review:** one optimization execution is the unit of ORE.
The workflow supplies immutable physical request identity, task permissions,
existing writer ownership and a fresh native execution directory. ORE may own
only its uncommitted child evidence namespace, proposed as
`<native-execution>/recovery/<stable-job-id>/attempt_NNN/` where the parent does
not already identify a job. Within an owned job directory use
`recovery/attempt_NNN/` directly. The enclosing workflow attempt number and ORE
attempt number are separate; no second enclosing workflow attempt hierarchy.
R1 only observes; R2 must approve exact per-route paths before execution changes.

Keep native scratch and accepted snapshots within that child attempt. Use
`pyar.persistence` canonical JSON, ArtifactRef/SHA-256 and immutable publication.
Workflow request/stage/selection records remain authoritative; diagnostics and
outcome reports are derived. Do not overwrite an artifact referenced by an
existing commit. Reuse the workflow RunLock; standalone pilot must acquire an
appropriate single writer before creating output. No new process-wide cwd
patterns, implicit v1 adoption, symlink escapes, unknown schema acceptance or
report-driven resume. Cross-host locking and thread safety are not implied.

Alternatives rejected: global retry loop around `optimise` (misses path engines
and workflow semantics); recovery reports as new workflow state (competing
restart authority); moving existing committed native files; accepting all
energy-bearing XYZ files as checkpoints.

## Prioritized findings and regression evidence

All fixtures below are in `tests/test_optimization_recovery_audit.py`; fake
executables are generated under pytest temporary directories. Strict xfails
assert the proposed safe boundary. They characterize legacy defects; they do
not authorize silently changing a v1 API. R1 may apply the safe assertion to
an opt-in facade while retaining legacy behavior until its migration is approved.

| ID / priority | Reproducer and observed behavior | Required direction |
|---|---|---|
| R0-01 / P0 | ORCA exit 0 + normal termination + finite XYZ/energy returns True **without** optimization convergence | Independent native convergence evidence |
| R0-02 / P0 | ORCA missing energy becomes 0.0; NaN energy accepted with convergence text | Missing is null/invalid, finite scalar required; genuine physical zero is not globally forbidden |
| R0-03 / P0 | ORCA reuses old wrong-symbol XYZ after current process omits geometry | Fresh output identity and ordered atom/mapping validation |
| R0-04 / P1 | Single atom np.loadtxt shape `(3,)` crashes writer with TypeError | Normalize and validate `(N,3)` without weakening atom-count/order checks |
| R0-05 / P0 | Stale xTB `.xtboptok` + old energy/XYZ returns True, even with current `.sccnotconverged` | Fresh-attempt provenance, failure precedence, current SCC/optimizer evidence |
| R0-06 / P0 | Generic result cache succeeds without convergence metadata, ignores conflicting method record and wrong atom symbols | Compatible request/frame/convergence proof in opt-in path; v1 policy review |
| R0-07 / P0 | First geomeTRIC call returns CycleExceeded; next wrapper call returns True from same provisional result and matching physical-energy state | Cached provenance must retain incomplete status |
| R0-08 / P0 | geomeTRIC accepts NaN physical energy, mismatched state/XYZ, zero exit with no convergence evidence, or unrelated fallback XYZ | Validate objective/accepted frame/finite values and actual solver convergence |
| R0-09 / P1 | BackendExecutionError clears original in-memory geometry | Preserve immutable input and independently validated accepted checkpoint; raw failure artifacts remain evidence |
| R0-10 / P1 | Truncated cached XYZ raises SystemExit through legacy reader | New facade must return/raise typed invalid-artifact failure, not terminate a batch process |

Additional inspected risks, not newly proven external-backend defects in R0:
AIMNet2/ASE helpers ignore optimizer run result; MOPAC arc presence; composite
model provenance; ORCA scan unconstrained stationarity; mutable generic state
written on trial evaluations; resource failures collapsed to False/None; normal
success not guaranteeing fixed constraints/identity. These need route-specific
fixtures before R1 claims support. A valid converged geometry can still be the
wrong scientific structure; generic bond-connectivity rejection would also be
wrong for some reactions/dissociating aggregates.

Positive controls: ORCA nonzero exit does not promote frame; xTB SCC-only marker
and zero-exit unknown option reject; geomeTRIC limit retains finite provisional
frame without `is_success`; absent final state raises BackendExecutionError;
Sella exception preserves original cause; EnergyGradientResult rejects nonfinite
energy/gradient; incompatible adaptive checkpoint rejects; path-stage reuse
rejects changed GFN model, tampered geometry hash and running summary. Existing
NEB, scan and persistence suites cover scientific gates, stale dependencies,
constraints, artifact hashes and active/stale lock behavior.

## Proposed R1 contract and review decisions

Keep the public `optimise(molecule, qc_params)` status API intact initially.
Add an opt-in facade around one execution and compact frozen types, using
existing persistence types rather than a second storage layer:

- `OptimizationRequest`: ordered atoms/mapping, starting geometry ref/hash,
  charge/multiplicity/reference/root where known, calculator/model/basis/assets,
  optimizer/target, constraints, solvation/wall/AFIR objective and workflow role;
  distinct physical/request and numeric-solver fingerprints.
- `OptimizationAttempt`: ID/parent, input and accepted frame refs, process exit,
  optimizer termination, versions, settings, diagnostic evidence, evaluations/
  steps/time and explicit unknowns. No trial geometry masquerading as accepted.
- `FailureDiagnostic`: stable code, evidence source/path/offset, severity,
  confidence and root-cause precedence; no guessed action from UNKNOWN.
- `OptimizationResult`: `CONVERGED`, `RECOVERED_CONVERGED` (reserved for R2),
  `INCOMPLETE_USABLE`, `FAILED_RECOVERABLE`, `FAILED_FINAL`, `INVALID_RESULT`;
  converged boolean, failure reason, accepted geometry and physical/effective
  energy with units, force residual/threshold evidence, attempt list and
  workflow gate evidence separately. A numerical success is not an index/IRC gate.

R1 parsers should observe only ORCA/xTB/geomeTRIC pilot outputs first. Validate
N×3 finite coordinates, ordered symbols, scalar finite energy, N×3 finite
gradients where available, matching frame provenance, requested tolerance and
physical/effective objective separation. A zero energy is not rejected solely
for being zero; a missing value filled with zero is invalid. Unsupported output
versions/evidence remain unknown/diagnose-only. Extend capabilities only with
proven parsed-convergence, accepted-frame, restart and allowed numeric-control
support. No retry controller or automatic corrective actions in R1.

Review decisions before R1:

1. Approve typed distinctions and the migration treatment of incomplete results
   in standalone optimize, final conformer/microsolvation selection versus AFIR/
   growth intermediate pools. Do not globally redefine is_usable.
2. Approve independent ORCA convergence/energy validity and cache-safe facade
   behavior; decide whether any legacy correctness change belongs in a separate
   R0b, rather than mixing it into observational R1.
3. Approve physical/numerical fingerprint fields and accepted-frame provenance,
   storage/lock ownership in ADR-ORE-001, and v1 uncertainty handling.
4. Defer recovery policy, safe action selection, interruption replay and budgets
   to R2 after R1 acceptance. The prior prototype's original-geometry reruns,
   finite-only checks and retry CLI are not an accepted implementation of R2.

## Validation ledger

Run from the isolated baseline checkout using the provisioned interpreter
(the `python` below denotes that exact interpreter; package import location was
checked before tests). Complete process exits and log tails were inspected.

- Before new fixtures: `python -m pytest -q` — 1461 passed, 1 skipped,
  691 warnings, 194 subtests passed, exit 0.
- New fixtures: `python -m pytest -q tests/test_optimization_recovery_audit.py`
  — 12 passed, 17 strict xfailed, exit 0. Expected exceptions are restricted
  so unrelated import/runtime errors cannot be silently hidden by xfail.
- Unmasked diagnostic run: `python -m pytest -q --runxfail
  tests/test_optimization_recovery_audit.py` — 17 failed, 12 passed, exit 1
  intentionally. Inspected each failure: unsafe success assertions, unraised
  BackendExecutionError, one-atom TypeError, lost input assertion, and truncated
  cache SystemExit. This is reproduction evidence, not a passing gate.
- Final full suite: `python -m pytest -q -rs` — **1473 passed, 1 skipped,
  17 xfailed, 691 warnings, 194 subtests passed**, exit 0 (16.58 seconds).
  Skip: packaging metadata test expects a wheel in `dist/`; R0 builds go to a
  temporary output directory. The built wheel was inspected separately below.
- `python -m sphinx -E -b html -W --keep-going docs /tmp/ore-r0-docs` —
  exit 0, build succeeded. The new audit is Markdown, which this RST Sphinx
  configuration does not render; its tables, source anchors and inventory were
  reviewed directly. No new Sphinx extension/dependency was added for the audit.
- `python -m build --outdir /tmp/ore-r0-dist` — exit 0,
  source archive and wheel both built for 1.5.1. Build output is outside Git.
- Wheel/source contents inspected for expected package modules, absence of
  vendored model binaries, and R0 test inclusion in the sdist. The new Markdown
  audit and TSV are repository review artifacts; the existing sdist manifest
  does not include these docs. No manifest change was made. `git diff --check`
  passed. Sphinx-generated `generated_api.rst` restored after inspection; the
  test-generated trial XYZ was moved outside the checkout.

External negative-path tests run generated Python executables, not real ORCA,
xTB or Sella chemistry. Existing analytic Sella tests, where run, are analytic
calculator checks. The full-suite skip reason and final command counts are
recorded with final validation. No recovery-rate or cost improvement is claimed.

**Next unstarted stage: R1. Stop for R0 review; no merge or release is authorized
by this audit.**
