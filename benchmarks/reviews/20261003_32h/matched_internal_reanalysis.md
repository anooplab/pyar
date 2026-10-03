# ORCA OptTS comparison on the RGD1-TSopt-GFN2 pilot

## Integrity and protocol

- Paired starts: 45; reaction clusters: 15.
- All three methods used the exact same starting geometry for every case (SHA-256 checked).
- The first ORCA record reports ORCA Program Version 6.1.1  -  RELEASE   - and xTB `6.7.1-6.fc44.fc44`.
- PyAR revision recorded by the optimizer runs: `57fb99a6d086eba3a091c5be1a368c4964ed7f15`.
- Provenance incomplete; unverified fields: irc_endpoint_rmsd_tolerance, product hash, product_relaxation_fmax, product_relaxation_max_steps, reactant hash, reference_ts hash.
- All three optimizers used internal coordinates: ORCA redundant internal coordinates, geomeTRIC delocalized internals, and Sella internal coordinates.
- All used a 200-cycle limit and a common maximum-gradient target of `0.02 eV/angstrom`. ORCA's maximum and RMS gradient cutoffs were converted from that target; its native energy and displacement criteria remained enabled, so the full convergence definitions are not identical.
- ORCA used `InHess XTB2` and Bofill updates. PyAR used its existing Hessian setup for geomeTRIC and Sella.
- Wall-time comparisons describe the supplied runs; contemporaneous execution and equal machine load have not been verified.
- Reaction success requires independent frequency, IRC, and endpoint validation. Missing provenance prevents verification that these checks used an identical protocol.
- ORCA does not expose provider calls in the same way as the PyAR energy-gradient adapter. ORCA optimizer timing-record counts are descriptive only and are not equated to PyAR provider evaluations.

## Reaction-connected outcomes

| Method | Validated starts | Rate | Reaction groups with at least one validated start |
|---|---:|---:|---:|
| ORCA OptTS | 9/45 | 20.0% | 4/15 |
| geomeTRIC | 10/45 | 22.2% | 4/15 |
| Sella | 11/45 | 24.4% | 4/15 |

| Paired validated-start result | geomeTRIC | Sella |
|---|---:|---:|
| Both / ORCA only / PyAR only / neither (geomeTRIC) | 9 / 0 / 1 / 35 | — |
| Both / ORCA only / PyAR only / neither (Sella) | 9 / 0 / 2 / 34 | — |

The paired reaction-cluster bootstrap compares each method's success fraction within each reaction, then resamples the reaction groups present in these results. A cluster-level sign-flip test is also reported.

| ORCA minus PyAR | Mean paired difference | Reaction-cluster bootstrap 95% interval | Sign-flip p (method) |
|---|---:|---:|---:|
| geomeTRIC | -2.2% | -6.7% to +0.0% | 1.000 (exact) |
| Sella | -4.4% | -11.1% to +0.0% | 0.500 (exact) |

Direct paired optimizer comparison (geomeTRIC minus Sella):
mean reaction-cluster difference -2.2%; bootstrap 95% interval -6.7% to +0.0%; exact sign-flip p=1.000.

## Outcomes and cost

| Method | Outcome counts | Median optimizer cycles | Median optimizer wall time (s) | Median ORCA gradient-log events / PyAR provider evaluations |
|---|---|---:|---:|---:|
| ORCA OptTS | `{'converged_not_stationary': 4, 'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 25, 'max_iterations': 3, 'reaction_connected_success': 9, 'stationary_not_first_order_saddle': 1}` | 46.000 | 1.980 | 45.000 |
| geomeTRIC | `{'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 20, 'optimizer_not_converged': 6, 'reaction_connected_success': 10, 'stationary_not_first_order_saddle': 6}` | 44.000 | 4.034 | 126.000 |
| Sella | `{'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 28, 'irc_exception': 1, 'reaction_connected_success': 11, 'stationary_not_first_order_saddle': 2}` | 14.000 | 3.177 | 31.000 |

## Conclusion

Validated reaction-connected outcomes were: ORCA OptTS 9/45, geomeTRIC 10/45, Sella 11/45. The corresponding counts of reaction groups with at least one validated start were: ORCA OptTS 4/15, geomeTRIC 4/15, Sella 4/15. All three methods agreed on success or failure for 43/45 paired starts. The paired reaction-cluster estimates and their intervals above describe this pilot; they do not establish that the methods are equivalent or that one is generally better.

ORCA's optimizer converged for 42/45 starts, while full validation determined whether each saddle connected the intended endpoints. An optimizer convergence message alone is not sufficient; frequency, IRC, and endpoint checks remain necessary.

Optimizer wall times are descriptive for these recorded runs. ORCA's gradient-log event count is not equivalent to PyAR provider evaluations, so the count columns must not be read as directly comparable electronic-structure call counts.

This pilot alone does not establish a general default optimizer; broader validation is still needed.

## Provisional default and recovery recommendation

Retain the current default pending broader validation. Internal-coordinate Sella validated 11/45 starts versus 10/45 for geomeTRIC. Reaction-group coverage is reported above.
Keep internal-coordinate Sella available as an explicit retry after geomeTRIC fails end-to-end validation. It added one validated start among the 35 geomeTRIC failures here. Retry from the original TS guess, then run frequency, IRC, and endpoint checks again. The gain is small, so automatic retries should wait for broader validation.
The comparison does not isolate the effect of P-RFO from Hessian initialization, updates, or energy-gradient call overhead. ORCA timing counters are not directly comparable to PyAR provider calls.
