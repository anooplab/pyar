# MS4 implementation smoke check

Date: 2026-10-17
Purpose: verify the opt-in implementation, durable attempt records, and strict
backend-evaluation cap. This is a software smoke test, not a basin-recovery
benchmark.

## Successful perturbation

System: neutral methane plus one neutral water, cluster objective, GFN2-xTB
6.7.1-6.fc44, geomeTRIC, one thread, one resampling move, 200-call cap, 600 s
limit, seed 23. The worktree was explicitly placed on `PYTHONPATH` so the
geomeTRIC subprocess imported the code being tested.

The growth survivor energy was -9.247282609368 Hartree. One rigid-fragment
perturbation and reoptimization completed using 26 backend energy/gradient
calls; the final survivor energy was -9.247286435935 Hartree. The energy change
was -0.000003826567 Hartree (about -0.0024 kcal/mol). This is one trajectory
and provides no evidence of a newly recovered basin or a search advantage.

## Hard-cap and clash checks

With the same starting structure and seed but a one-call resampling budget, the
optimizer attempted a second evaluation after the first reservation. The
calculator stopped it before the backend call; the persistent counter remained
exactly 1, the attempt was recorded as `budget_exhausted`, and the completed
growth survivor remained available as the final result.

A second seed generated a proposal rejected by the interfragment clash filter.
It used zero backend calls and was recorded as `rejected_clash`.

Raw run records and optimizer logs are retained outside the repository under
`/home/anoop/tmp/ms4-live-cluster-run2`,
`/home/anoop/tmp/ms4-live-budget-run2`, and
`/home/anoop/tmp/ms4-live-budget-run`.

## Scope and limitation

These checks exercise one-water rigid-body perturbation, physical-energy
reporting, the evaluation cap, clash rejection, and result fallback. They do
not test multi-water basin recovery, independent trial statistics, or global
minimum recovery. Those require a separate matched-budget benchmark before
making a scientific benefit claim.
