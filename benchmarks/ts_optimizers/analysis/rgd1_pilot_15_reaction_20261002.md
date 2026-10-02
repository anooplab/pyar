# RGD1-TSopt-GFN2 15-reaction pilot analysis

## Run integrity

- Run output: `pilot_committed_2591ed0` (raw job artifacts are stored locally, not in the repository)
- Benchmark manifest SHA256: `4bd95c23742b60194f82809a6d408b8ea38d8d8f6869c0e189427095113339f1`
- PyAR commit: `2591ed0569ff2d20b350d73b0c0844a73a15a810`
- Runtime: Python 3.14.7 (main, Aug 10 2026, 00:00:00) [GCC 16.1.1 20260515 (Red Hat 16.1.1-2)]; geomeTRIC 1.1.1; Sella 2.6.0; xtb version 6.7.1-6.fc44.fc44 compiled by Fedora project on Sun Jan 18 23:32:34 UTC 2026
- Shared TS settings: `ts_fmax=0.02 eV/angstrom`, `ts_max_cycles=200`
- Completed jobs: 90 / 90; incomplete: 0
- Paired inputs with matching hashes: 45 / 45
- Experimental units: 15 reactions × 3 correlated difficulty tiers (45 starts).

## Reaction-connected success

Success requires PyAR endpoint-connection validation, not only optimizer convergence or a first-order saddle.

| Optimizer | Successes | Rate | Reaction-cluster bootstrap 95% interval |
|---|---:|---:|---:|
| geometric | 10/45 | 22.2% | 4.4%–42.2% |
| sella | 10/45 | 22.2% | 4.4%–42.2% |

| Paired outcome | Cases |
|---|---:|
| Both succeed | 9 |
| geomeTRIC only succeeds | 1 |
| Sella only succeeds | 1 |
| Neither succeeds | 34 |

Reaction-cluster paired success-rate difference (geomeTRIC − Sella): +0.0%; 95% cluster-bootstrap interval -6.7% to +6.7%; exact reaction-level sign-flip p = 1.000 (n = 15 reaction clusters).

| Tier | Optimizer | Reaction-connected successes | Other outcomes |
|---|---|---:|---|
| easy | geometric | 4/15 | endpoint_relaxation_exception: 1, first_order_saddle_wrong_connection: 8, optimizer_not_converged: 2 |
| easy | sella | 4/15 | endpoint_relaxation_exception: 1, first_order_saddle_wrong_connection: 10 |
| hard | geometric | 2/15 | endpoint_relaxation_exception: 1, first_order_saddle_wrong_connection: 7, optimizer_not_converged: 2, stationary_not_first_order_saddle: 3 |
| hard | sella | 3/15 | first_order_saddle_wrong_connection: 8, optimizer_not_converged: 1, stationary_not_first_order_saddle: 3 |
| med | geometric | 4/15 | endpoint_relaxation_exception: 1, first_order_saddle_wrong_connection: 5, optimizer_not_converged: 2, stationary_not_first_order_saddle: 3 |
| med | sella | 3/15 | endpoint_relaxation_exception: 1, first_order_saddle_wrong_connection: 10, stationary_not_first_order_saddle: 1 |

## Cost and reference agreement

Cost values are descriptive across all 45 jobs (mean / median). Validation cost is shown separately because it can dominate TS optimizer work and depends on whether earlier scientific gates passed.

| Measure | geomeTRIC | Sella |
|---|---:|---:|
| TS provider evaluations (mean / median) | 154.533 / 126.000 | 106.711 / 64.000 |
| TS optimizer steps (mean / median) | 81.933 / 44.000 | 44.711 / 26.000 |
| TS wall time (s) (mean / median) | 6.345 / 4.034 | 2.941 / 2.137 |
| Total provider evaluations (mean / median) | 502.844 / 357.000 | 561.422 / 590.000 |
| Validation provider evaluations (mean / median) | 401.897 / 285.000 | 465.045 / 387.000 |
| Validation wall time (s) (mean / median) | 23.401 / 9.299 | 25.045 / 15.702 |

Reference geometry and energy comparison is available for successful runs only. The dataset reference uses tblite 0.6.0; this pilot used the xTB executable.

| Optimizer | Successful reference RMSD, median / max (Å) | Absolute reference energy difference, median / max (Eh) |
|---|---:|---:|
| geometric | 2.45697e-05 / 0.00068444 | 9.6725e-09 / 6.1248e-08 |
| sella | 0.00314857 / 0.00727965 | 2.92967e-06 / 6.80786e-06 |

## Interpretation and limits

- The pilot shows equal validated success (10/45 each). Only two pairs were discordant, one in each direction.
- Sella used fewer TS evaluations and less measured TS wall time on average, while mean end-to-end provider evaluations were higher. These are estimates from one host and this selected subset.
- The three starts per reaction are correlated. Chemical diversity is 15 reactions, not 45 independent reactions; per-tier rates are descriptive.
- Success intervals resample whole reactions. The paired test flips optimizer labels at the reaction level.
- These results do not justify selecting a default. Continue to the broader paired benchmark and native ORCA reference first.
- GFN2-xTB executable results are not claimed numerically identical to the tblite reference implementation.
