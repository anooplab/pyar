# ORCA OptTS comparison on the RGD1-TSopt-GFN2 pilot

## Integrity and protocol

- Paired starts: 45; reaction clusters: 15.
- All three methods used the exact same starting geometry for every case (SHA-256 checked).
- ORCA OptTS used Program Version 6.1.1  -  RELEASE   - and the ORCA GFN2-xTB external-method interface with the same xTB executable used by PyAR: `6.7.1-6.fc44.fc44`.
- PyAR revision recorded by the optimizer runs: `57fb99a6d086eba3a091c5be1a368c4964ed7f15`.
- A same-geometry single-point control for `MR_293227_1_easy` differed by approximately `9e-12 Eh` between ORCA's GFN2-xTB interface and PyAR's xTB provider.
- All three optimizers used internal coordinates: ORCA redundant internal coordinates, geomeTRIC delocalized internals, and Sella internal coordinates.
- All used a 200-cycle limit and a common maximum-gradient target of `0.02 eV/angstrom`. ORCA's maximum and RMS gradient cutoffs were converted from that target; its native energy and displacement criteria remained enabled, so the full convergence definitions are not identical.
- ORCA used `InHess XTB2` and Bofill updates. PyAR used its existing Hessian setup for geomeTRIC and Sella.
- The geomeTRIC result batch was reused from an earlier run with the same starts, force target, cycle limit, and delocalized internal coordinates. Sella and ORCA were run in a later batch. Wall-time comparisons are therefore descriptive across separate batches, not a same-session timing experiment.
- Every converged optimizer result was passed through the same PyAR frequency, IRC, and endpoint validation stages. A converged geometry alone does not count as reaction success.
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

The paired reaction-cluster bootstrap compares each method's success fraction within each reaction (three difficulty starts per reaction), then resamples the 15 reactions. The cluster-level exact sign-flip test is also reported.

| ORCA minus PyAR | Mean paired difference | Reaction-cluster bootstrap 95% interval | Exact sign-flip p |
|---|---:|---:|---:|
| geomeTRIC | -2.2% | -6.7% to +0.0% | 1.000 |
| Sella | -4.4% | -11.1% to +0.0% | 0.500 |

Direct paired optimizer comparison (geomeTRIC minus Sella):
mean reaction-cluster difference -2.2%; bootstrap 95% interval -6.7% to +0.0%; exact sign-flip p=1.000.

## Outcomes and cost

| Method | Outcome counts | Median optimizer cycles | Median optimizer wall time (s) | Median ORCA gradient-log events / PyAR provider evaluations |
|---|---|---:|---:|---:|
| ORCA OptTS | `{'converged_not_stationary': 4, 'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 25, 'max_iterations': 3, 'reaction_connected_success': 9, 'stationary_not_first_order_saddle': 1}` | 46.000 | 1.980 | 45.000 |
| geomeTRIC | `{'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 20, 'optimizer_not_converged': 6, 'reaction_connected_success': 10, 'stationary_not_first_order_saddle': 6}` | 44.000 | 4.034 | 126.000 |
| Sella | `{'endpoint_relaxation_exception': 3, 'first_order_saddle_wrong_connection': 28, 'irc_exception': 1, 'reaction_connected_success': 11, 'stationary_not_first_order_saddle': 2}` | 14.000 | 3.177 | 31.000 |

## Sella coordinate-system comparison

The same 45 starting geometries were also run with Sella in Cartesian coordinates. Validated successes were 10/45 for Cartesian and 11/45 for internal coordinates; both succeeded on 10, Cartesian alone on 0, internal alone on 1, and neither on 34. The full outcome label changed for 8/45 starts. This shows that coordinate choice can affect outcomes, but this pilot does not establish a universal setting.

## Conclusion

Validated reaction-connected outcomes were: ORCA OptTS 9/45, geomeTRIC 10/45, Sella 11/45. The corresponding counts of reaction groups with at least one validated start were: ORCA OptTS 4/15, geomeTRIC 4/15, Sella 4/15. All three methods agreed on success or failure for 43/45 paired starts. The paired reaction-cluster estimates and their intervals above describe this pilot; they do not establish that the methods are equivalent or that one is generally better.

ORCA's optimizer converged for 42/45 starts, while full validation determined whether each saddle connected the intended endpoints. An optimizer convergence message alone is not sufficient; frequency, IRC, and endpoint checks remain necessary.

Optimizer wall times are descriptive for these recorded runs. ORCA's gradient-log event count is not equivalent to PyAR provider evaluations, so the count columns must not be read as directly comparable electronic-structure call counts.

This is a 15-reaction pilot. It is not enough evidence by itself to establish a general default optimizer; broader validation is still needed.

## Provisional default and recovery recommendation

Keep geomeTRIC as the default for now. Internal-coordinate Sella validated 11/45 starts versus 10/45 for geomeTRIC, a one-start difference; both reached at least one validated result for the same number of reaction groups. This pilot does not justify changing the default.
Keep internal-coordinate Sella available as an explicit retry after geomeTRIC fails end-to-end validation. It added one validated start among the 35 geomeTRIC failures here. Retry from the original TS guess, then run frequency, IRC, and endpoint checks again. The gain is small, so automatic retries should wait for broader validation.
ORCA was faster by median optimizer wall time in these recorded batches, but it had no validated success unique to ORCA and its timing counters are not comparable to PyAR provider calls. Since the coordinates were matched, this timing result cannot be explained by coordinate choice alone; the comparison also does not isolate the effect of P-RFO from Hessian initialization/updates or energy-gradient call overhead. It does not support replacing the current default with ORCA.
