ORE R5 benchmark and validation plan
====================================

Status and boundary
-------------------

R5 benchmark execution and chemistry-level validation are deferred by project
direction (2026-10-18). This document records the work required to resume that
stage. No recovery-rate, false-success, diagnostic-accuracy, cost-saving, or
chemistry-correctness result is claimed here.

R4 added observational records to conformer/growth, microsolvation,
reaction/AFIR, bond-scan, NEB/TS/frequency/IRC/endpoint, and standalone modern
optimization routes. Workflow owners still decide whether results are usable.
The only automatic recovery route currently available remains the explicitly
opted-in SAFE/GUIDED xTB + geomeTRIC minimum route, with its existing policy
gates. R4 diagnostics do not widen backend, optimizer, target, or workflow
recovery permissions.

R4 software pull requests ran the focused tests and repository CI appropriate
to each change, including the full test matrix, installation checks and
documentation build. Those checks verify software behavior; they do not replace
the R5 benchmark or validate chemistry outcomes.

Benchmark corpus to freeze before execution
-------------------------------------------

Create a versioned, immutable manifest. Each case must identify its input
structure and SHA-256, requested objective and controls, backend and version,
software environment, expected baseline outcome, whether failure is
recoverable, expected diagnostic class, and workflow-owned chemistry gates.
Keep provenance, licenses, and source references with the manifest. Do not
classify a case as recoverable after seeing the ORE result.

Include representative cases from each enabled route and hard negative
controls:

* Minimum optimizations: ordinary convergence controls; iteration-limit and
  trust/line-search failures; flat or noisy surfaces; SCF nonconvergence and
  wrong-spin/root warnings; malformed or non-finite energies, gradients and
  geometries; truncated output; termination signals; stale or incomplete
  artifacts.
* Conformer and aggregate/growth workflows: controls for ordinary candidates,
  weakly bound or dissociating structures, duplicate basins, and structures
  whose selection/ranking must remain unchanged. Include the established
  CH4(H2O)n microsolvation sizes and representative noncovalent clusters.
* Microsolvation: solvent-count and fragment identity controls, preserved
  original solute indices, shell membership and site/cluster objective checks,
  biased search followed by unbiased relaxation, and negative controls that
  fail shell or physical-energy acceptance.
* Reaction and scan: AFIR orientation and gamma cases, adaptive release
  evidence, unbiased relaxation that retains or loses forming bonds, product
  identity controls, scan-distance adherence, and a scan whose relaxed result
  must be rejected by the workflow gate.
* Pathway and TS work: failed and converged NEB paths, endpoint connectivity,
  interior maximum/climbing-image selection, TS optimizer cycle limits,
  wrong-index or non-first-order saddles, frequency failures, IRC failures,
  and endpoint identity/connectivity failures. A minimum optimizer must never
  count as TS recovery.
* Infrastructure: interrupted workers, crash/resume boundaries, wall-time and
  evaluation limits, missing executables/dependencies, malformed checkpoints,
  stale locks, disk exhaustion, output tampering, and unsupported backend or
  output versions.

For every recoverable case include a matched no-recovery run. Include
non-recoverable negative controls specifically designed to catch false
success, changed physical objectives, wrong molecular identity, cache reuse,
and invalid TS/pathway outcomes. Unsupported or unrun routes must remain
``NOT_RUN``/unknown; they must never be counted as ``PASS``.

Paired execution and evidence
-----------------------------

Freeze the protocol before calculating: PyAR commit, dependency lock and
backend executable versions; operating system and hardware; environment
variables and thread limits; seeds; requested thresholds; per-case attempt,
cycle, evaluation, and wall-time limits; and rules for reruns and exclusions.
For each case, run the legacy/no-recovery control and ORE under matched
conditions. Preserve raw stdout/stderr, native output, state/checkpoint files,
diagnostic reports, run manifests, and exact commands. Hash inputs and result
artifacts. Record every attempt, including failed, interrupted and excluded
runs, with the reason for exclusion. Do not discard costs from unsuccessful
attempts.

Repeat stochastic cases enough times to estimate variability and report paired
differences with uncertainty intervals. Use the same seeds in each paired arm
where the workflow supports them. Have chemistry reviewers assess a blinded
subset and all apparent recovered successes against the case's predefined
identity, connectivity, electronic-state, objective, and pathway criteria.

Metrics and definitions
-----------------------

Publish counts, denominators, missingness, and uncertainty for every metric.
Use the frozen corpus labels as the reference truth and retain ``unknown`` as a
separate outcome.

* **Verified recovery rate:** recoverable reference cases that end in an
  independently verified, workflow-accepted result after ORE, divided by all
  reference recoverable cases attempted. Also report the conditional rate
  among cases where the baseline actually failed.
* **False-success rate:** ORE runs reported as successful but rejected by the
  independent objective, electronic-state, geometry, identity, constraint, or
  workflow gate, divided by all ORE runs reported successful. Report the raw
  false-success count and each gate-specific cause.
* **Diagnostic accuracy:** confusion matrix and per-class precision/recall for
  the frozen failure labels, with ``unknown`` and unsupported routes shown
  explicitly rather than forced into a known class.
* **Added cost:** paired CPU/GPU time, backend energy/gradient calls,
  optimizer steps/cycles, process attempts, wall time, peak memory where
  available, and monetary cost when a documented price source exists. Report
  totals and distributions for both successful and unsuccessful recoveries.
* **Chemistry correctness:** objective and settings identity; charge, spin and
  electronic-root checks; finite, ordered and chemically identified geometry;
  preserved workflow constraints; product/solute/fragment identity;
  noncovalent shell and unbiased-energy gates; TS index, frequency, IRC and
  endpoint gates. Report any gate changed or bypassed as a failure.
* **Reliability and compatibility:** repeat-run determinism; crash/resume
  equivalence to uninterrupted execution; checkpoint rejection when provenance
  differs; legacy status/API compatibility; and behavior for unsupported
  backends and output versions.

Report results by workflow, backend, optimizer, diagnostic class and action;
do not present an aggregate average that hides a false success or a failing
chemistry gate.

Acceptance and release evidence
-------------------------------

Before a release decision, provide the frozen corpus and protocol, raw
machine-readable results, hashes and software versions, analysis code, a
human-readable report, and review notes for false positives and excluded
cases. A proposed recovery action passes only if it demonstrates verified
recovery on its declared recoverable subset, no false success in the negative
controls, preserved workflow chemistry gates, bounded added cost, deterministic
restart behavior, and backward-compatible outcomes. Set numeric minimums for
recovery and maximums for cost before unblinding results; do not choose them
after seeing the measurements.

Run the final release suite after benchmark-driven changes: complete tests,
strict Sphinx build, source and wheel builds, clean installation of the built
artifacts, supported Python-version checks, and relevant CLI ``--help`` and
workflow smoke checks. Record exact commands, exit statuses, versions, and
artifact locations. Chemistry smoke runs supplement these checks; they do not
replace the corpus or justify claims beyond the tested cases.

Resume sequence
---------------

#. Review and freeze the corpus manifest, licenses, labels and protocol.
#. Set predeclared numeric acceptance thresholds and cost ceilings.
#. Run matched baseline and ORE arms, retaining every attempt and raw artifact.
#. Compute the metrics above and independently review all apparent recoveries.
#. Resolve every false success, unknown, unsupported route and reproducibility
   discrepancy before any recovery permission is expanded.
#. Run full release validation and make a documented go/no-go decision.

Until those steps are completed, R5 remains pending and R4 diagnostics must be
read as observations rather than evidence that a general recovery action is
chemically safe or cost-effective.
