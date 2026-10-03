# ORCA OptTS comparison on the RGD1-TSopt-GFN2 pilot

## Integrity and protocol

- Paired starts: 45; reaction clusters: 15.
- All three methods used the exact same starting geometry for every case (SHA-256 checked).
- ORCA OptTS used Program Version 6.1.1  -  RELEASE   - and the ORCA GFN2-xTB external-method interface with the same xTB executable used by PyAR: `6.7.1-6.fc44.fc44`.
- PyAR revision recorded by the optimizer runs: `57fb99a6d086eba3a091c5be1a368c4964ed7f15`.
- A same-geometry single-point control for `MR_293227_1_easy` differed by approximately `9e-12 Eh` between ORCA's GFN2-xTB interface and PyAR's xTB provider.
- ORCA used `InHess XTB2`, Bofill updates, and ORCA's default OptTS convergence thresholds. PyAR used its shared `ts_fmax=0.02 eV/angstrom` protocol. Therefore optimizer convergence thresholds are not identical.
- Every converged optimizer result was passed through the same PyAR frequency, IRC, and endpoint validation stages. A converged geometry alone does not count as reaction success.
- ORCA does not expose provider calls in the same way as the PyAR energy-gradient adapter. ORCA optimizer timing-record counts are descriptive only and are not equated to PyAR provider evaluations.

## Reaction-connected outcomes

| Method | Validated starts | Rate | Reaction groups with at least one validated start |
|---|---:|---:|---:|
| ORCA OptTS | 9/45 | 20.0% | 4/15 |
| geomeTRIC | 10/45 | 22.2% | 4/15 |
| Sella | 10/45 | 22.2% | 4/15 |

| Paired validated-start result | geomeTRIC | Sella |
|---|---:|---:|
| Both / ORCA only / PyAR only / neither (geomeTRIC) | 9 / 0 / 1 / 35 | — |
| Both / ORCA only / PyAR only / neither (Sella) | 9 / 0 / 1 / 35 | — |

The paired reaction-cluster bootstrap compares each method's success fraction within each reaction (three difficulty starts per reaction), then resamples the 15 reactions. The cluster-level exact sign-flip test is also reported.

| ORCA minus PyAR | Mean paired difference | Reaction-cluster bootstrap 95% interval | Exact sign-flip p |
|---|---:|---:|---:|
| geomeTRIC | -2.2% | -6.7% to +0.0% | 1.000 |
| Sella | -2.2% | -6.7% to +0.0% | 1.000 |

## Outcomes and cost

| Method | Outcome counts | Median optimizer cycles | Median optimizer wall time (s) | Median ORCA gradient-log events / PyAR provider evaluations |
|---|---|---:|---:|---:|
| ORCA OptTS | `{'converged_not_stationary': 4, 'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 25, 'max_iterations': 3, 'reaction_connected_success': 9, 'stationary_not_first_order_saddle': 1}` | 46.000 | 1.326 | 45.000 |
| geomeTRIC | `{'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 20, 'optimizer_not_converged': 6, 'reaction_connected_success': 10, 'stationary_not_first_order_saddle': 6}` | 44.000 | 4.034 | 126.000 |
| Sella | `{'endpoint_relaxation_exception': 2, 'first_order_saddle_wrong_connection': 28, 'optimizer_not_converged': 1, 'reaction_connected_success': 10, 'stationary_not_first_order_saddle': 4}` | 26.000 | 2.137 | 64.000 |

## Conclusion

ORCA OptTS did not improve validated success on this small pilot: it validated 9 starts, compared with 10 for either PyAR optimizer. All methods found at least one validated start for the same 4 of 15 reaction groups. All three methods agreed on success for 9 starts and failure for 34. There were two discordant starts: geomeTRIC alone succeeded on one and Sella alone on another; ORCA had no ORCA-only success. The reaction-cluster comparison has very wide uncertainty and does not establish that the methods are equivalent or that one is better.

ORCA's optimizer converged for 42/45 starts, but many resulting saddles did not connect the target endpoints. This reinforces that an ORCA OptTS convergence message is not sufficient: frequency, IRC, and endpoint checks remain necessary. ORCA's median optimizer wall time was lower on this host, but the convergence criteria differ and the ORCA energy/gradient counters are not directly comparable to PyAR provider-call counts.

This is a 15-reaction pilot, not enough evidence to choose a default TS optimizer. Broader and more chemically diverse validation is still needed.
