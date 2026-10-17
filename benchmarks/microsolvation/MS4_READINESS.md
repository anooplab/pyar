# MS4 readiness review

Date: 2026-10-17  
Repository baseline: `anooplab/pyar` `master` at `9cfb8a2` (MS3 PR #58).  
Decision: **REVISE** — do not start production MS4 yet.

This is a readiness and feasibility review, not an implementation of the MS4
search algorithm. The only production-code edits in this review restore the
historical default shell tie-breaker and separate automatic confinement from
candidate selection. The new benchmark drivers are opt-in external experiments.

## Stage acceptance matrix

| Stage | Acceptance question | Evidence reviewed | Status and remaining gap |
|---|---|---|---|
| MS0 baseline | Can fixed-composition cluster-growth results be recorded and compared reproducibly? | PR #54 and its merged benchmark harness/report; this review's two independent growth trials per solute; common GFN2-xTB final energies and basin analysis | **Partial pass.** Harness/provenance works, but two trials per case are a small feasibility sample. Historical ORCA/CREST/QCG comparisons are singleton or failed runs and are not an uncertainty estimate or global-minimum reference. |
| MS1 objective-aware search | Are shell and whole-cluster objectives explicit, and are cluster survivors ranked by fixed-composition physical energy? | PR #55; `microsolvation` request/selection code and tests; v1/v2 workflow paths | **Pass with correction in this review.** Cluster candidates retain composition and are energy ranked in Hartree. Default shell ranking had lost its historical energy tie-breaker; restored below. |
| MS2 site contacts | Are contact/saturation diagnostics explicit and kept separate from physical energy? | PR #56; contact-policy code, tests, documentation, and persisted diagnostics | **Pass for diagnostics.** Contacts are heuristic placement/surface diagnostics, not a thermodynamic objective. MS4 still needs solvent-fragment identity and atom-order-safe deduplication tests. |
| MS3 optional wall | Is confinement opt-in/case-aware, with correct gradient, unbiased relaxation, and physical-energy reporting? | PRs #57–58; wall and automatic-policy docs; tests; this review's all-coordinate finite-difference regression; prior live MS3 pilots | **Pass with selection-policy correction.** Explicit wall behavior remains on its existing fingerprinted policy. Automatic wall activation was improperly coupled to energy-first survivor selection even when no wall was applied; automatic policy is now independently versioned. |
| MS4 readiness | Does whole-cluster perturbation recover additional low-energy basins at matched backend cost, with exact restart/provenance and failure handling? | This review's two-trial cluster growth and external 40-gradient-call-per-arm perturbation pilot; existing restart and artifact tests | **Not demonstrated.** No novel basin appeared in the small perturbation pilot, and each perturb arm exhausted its budget inside an optimization. No global-minimum recovery claim is made. Production MS4 acceptance tests below remain open. |

PRs [#54](https://github.com/anooplab/pyar/pull/54),
[#55](https://github.com/anooplab/pyar/pull/55),
[#56](https://github.com/anooplab/pyar/pull/56),
[#57](https://github.com/anooplab/pyar/pull/57), and
[#58](https://github.com/anooplab/pyar/pull/58) are merged to current master.
Their merge commits are respectively `134ede4`, `9b38d8f`, `5f6dba2`,
`6e6a3a9`, and `9cfb8a2`; each had 52 successful checks at merge.

## Correctness findings and scoped fixes

The historical `_quality_key()` (before MS1) sorted shell candidates by
descending surface coverage, then ascending finite Hartree energy, then name.
MS1 removed the energy key from default geometrical shell ranking. This review
restores that tie-breaker and adds a direct regression test. Energy remains a
tie-breaker only after coverage, so a lower-coverage structure cannot displace
a higher-coverage shell candidate. Cluster ranking remains fixed-composition
energy-first.

Before this review, `wall="auto"` was treated as an energy-ranked wall request
by `_quality_key()` and the selection frontier, even when no candidate wall
was activated. Confinement is now separate from selection: auto shell ranking
is surface coverage, then physical energy, then diversity, independent of
whether the per-candidate wall activates; auto cluster ranking is fixed-
composition energy. Auto requests use `objective-aware-v2` selection
fingerprints. Existing automatic-wall requests created under the previous
policy will fail request-fingerprint matching and must be restarted as new
runs; this avoids silently changing an in-progress selection trajectory.
Default requests retain their historical serialization, and explicit
`outer-sphere` selection serialization is regression-tested as unchanged. No
active auto-wall state/request files were found in the searched
`Dropbox/PyAR-2.0` tree.

The wall gradient test now compares every Cartesian component for solute atoms
and atoms in multiple solvent fragments against central finite differences.
The workflow removes the wall contribution before final unbiased relaxation
and stores/ranks the final physical energy. The wall is applied as an optimizer
potential; it is not added to the final reported physical energy.

Other reviewed invariants: backend energy records and prescreen thresholds are
Hartree; every compared cluster candidate has the same solute plus four waters;
shell coverage filtering is distinct from cluster objective selection; the
solute's original atom sequence remains the prefix of the combined structure;
solvent fragments map to contiguous appended atom indices; fragment charges
sum to the requested total and electron-parity/multiplicity checks happen
before backend work; and optimization failures are retained in attempt/state
records rather than ranked as successful structures. MS4 should add explicit
permutation/reordering and solute-identity preservation checks: the existing
prefix/provenance invariants do not by themselves prove that optimization
preserves the solute's chemical connectivity.

Restart review found v1 continuation remains owned by its legacy path; v2 uses
request fingerprints, hashed relative artifact references, immutable stage
records, locking/recovery, and completed-run idempotence. Existing tests cover
v1 fixtures, input mismatch, artifact tampering, completed runs, and read-only
preflight. Result copies are not restart authority. A live read-only preflight
for cluster methane and acetate succeeded without creating an output run.
The auto selection fingerprint change above is the intentional compatibility
exception for in-progress auto runs.

## GFN2-xTB cluster-growth trials

Commands and raw run directories are under `/home/anoop/tmp/ms4-growth-*`;
reproduction scripts and fixtures are `run_ms4_trials.py` and
`benchmarks/microsolvation/fixtures/`. The runner uses GFN2-xTB 6.7.1-6.fc44,
Python 3.14.7, geomeTRIC, one thread, fixed charge/multiplicity, four waters,
four placement orientations, two deterministic rotation seeds (1103, 2207),
two independently initialized trials per solute, and a 300-cycle optimization
limit. Each invocation of `xtb --grad` is counted. All six trials completed
four growth steps and their optimizer logs contained no iteration-limit
warning. Composition, solute element prefix, Hartree energy records, and
artifact hashes were checked.

Basins are clustered among the four final selected structures per solute by
heavy-atom graph-constrained Kabsch RMSD at 0.5 Å. “Near minimum” means within
1.0 kcal/mol of the lowest energy observed in that small pool.

| Solute + 4 waters | Charge / multiplicity | xTB gradient calls (trial 1, 2) | Wall time (s) | Best observed energy (Eh) | Distinct observed basins | Basins within 1 kcal/mol |
|---|---:|---:|---:|---:|---:|---:|
| Methane | 0 / 1 | 2474, 2340 | 87.298, 85.695 | -24.504075049700 | 2 | 1 |
| Ammonia | 0 / 1 | 2247, 2063 | 71.267, 65.971 | -24.760206404524 | 2 | 2 |
| Acetate | -1 / 1 | 2354, 2353 | 87.742, 85.917 | -34.543113497987 | 4 | 1 |

These are observed minima, not validated global minima. Two deterministic
rotations are not an independent statistical ensemble; variation, confidence
intervals, and between-method convergence cannot be inferred. Energies are
compared only within each fixed-composition solute case, never across rows.
The first acetate launch failed before any xTB call because the harness
incorrectly gave acetate's -1 charge to water as well. It was preserved at
`/home/anoop/tmp/ms4-growth-acetate-300`, diagnosed from its preflight error,
and rerun with fragment-specific charges at
`/home/anoop/tmp/ms4-growth-acetate-300-corrected`. The corrected trials
completed. A separate 100-cycle methane smoke run emitted geomeTRIC iteration
limit warnings and is excluded from the numerical table; 300 cycles were used
for all reported growth results.

## External perturbation feasibility pilot

`ms4_perturbation_pilot.py` is a scratch driver outside production search code.
For each case it used the same best starting structure(s) for two arms:
unperturbed reoptimization and one-water rigid rotation plus 1–2 Å translation,
followed by collision screening (1.3 Å cutoff) and unbiased GFN2-xTB
reoptimization. Each arm had a hard budget of exactly 40 actual `xtb --grad`
calls, one thread, and 250 optimizer cycles. Selection is simply the lowest
energy seen, with no Monte Carlo acceptance. Basin/RMSD thresholds match the
growth comparison.

| Case | Arm | Calls | Attempts | Collisions rejected | Failed/incomplete optimizations | New basins beyond starts |
|---|---|---:|---:|---:|---:|---:|
| Methane | Baseline | 40 | 20 | 0 | 0 | 0 |
| Methane | Perturb | 40 | 4 | 2 | 1 | 0 |
| Ammonia | Baseline | 40 | 20 | 0 | 0 | 0 |
| Ammonia | Perturb | 40 | 2 | 0 | 1 | 0 |
| Acetate | Baseline | 40 | 16 | 0 | 0 | 0 |
| Acetate | Perturb | 40 | 2 | 1 | 1 | 0 |

In each perturb arm, the 40-call cap stopped a geomeTRIC optimization before
completion (one attempted backend call beyond the budget was blocked and
logged). Consequently, perturbation arms yielded fewer completed relaxations;
the comparison is matched in actual backend calls but not in completed
optimization count. It shows the need for resumable/budget-aware optimizer
steps and a principled proposal/acceptance design, not that rigid moves are
ineffective or effective. This one-water pilot does not test whole-cluster,
multi-solvent, exchange, exact restart, or global basin recovery.

## Remaining MS4 blockers and acceptance gates

Before production implementation, specify and test:

1. A budget-aware proposal/optimizer contract that can stop cleanly without
   losing partial work, distinguish rejected collisions from backend failures,
   and charge every backend evaluation to the declared budget.
2. Solvent identity and fragment mapping through translation, rotation,
   exchange, atom permutation, reordering, optimization, serialization, and
   restart. Verify that the original solute's identity/connectivity is
   preserved or explicitly detect and classify changes.
3. A seed-independent proposal stream and exact RNG/proposal-state recovery
   after interruption, with immutable attempt provenance and versioned request
   fingerprints.
4. Fair baseline comparisons from the same starting structures and matched
   actual GFN2-xTB evaluation budgets, across more independent trials and
   chemically distinct basins. Report cost, failures, diversity, and uncertainty;
   do not describe the lowest observed energy as a global minimum.
5. Tests for symmetric solvent exchange, multi-solvent moves, fragment
   permutation, collisions, charged systems, failed optimizations, budget
   exhaustion, duplicate detection after atom reordering, and v1/v2 restart
   compatibility.

**Recommendation: REVISE before MS4.** The foundation and two narrowly scoped
selection fixes are in good shape, and the external perturbation experiment is
feasible. The current pilot does not support a scientific benefit claim, and
its hard-budget interruption plus missing identity/restart acceptance evidence
should be resolved in the MS4 design and benchmark protocol before production
search code is authorized.
