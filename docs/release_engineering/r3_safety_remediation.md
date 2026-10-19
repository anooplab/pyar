# R3 optimization safety remediation

This remediation addresses the 17 strict expected-failure instances recorded by
the R3 closure review. Each expectation was removed only after the corresponding
test passed as an ordinary test. The tests use synthetic backend processes and
fixtures; they validate failure handling and artifact contracts, not chemical
accuracy of ORCA, xTB, or geomeTRIC calculations.

## Defect-to-fix record

| Review ID | Defect and regression coverage | Root fix |
| --- | --- | --- |
| R0-01 | ORCA normal process exit without optimization convergence; `test_orca_normal_termination_without_geometry_convergence_is_rejected` | Require normal termination and explicit optimization convergence before accepting a result. |
| R0-02a | Missing ORCA energy; `test_orca_missing_or_nonfinite_energy_is_rejected[]` | Missing energy remains absent and fails result validation; no zero energy is fabricated. |
| R0-02b | Nonfinite ORCA energy; same parameterized test (`nan`) | Reject nonfinite final energy. |
| R0-03 | Stale/wrong-symbol ORCA geometry; `test_orca_does_not_reuse_old_wrong_symbol_xyz` | Remove the canonical result before execution, preserve the prior frame as an input diagnostic, then validate fresh atom symbols/order/count and coordinates. |
| R0-04 | Single-atom ORCA coordinates flattened; `test_orca_one_atom_has_matrix_coordinates` | Parse XYZ into an explicit `(N, 3)` coordinate array. |
| R0-05a | Stale xTB `.xtboptok` before an SCC failure; parameterized `test_xtb_stale_success_cannot_authorize_current_run[False]` | Archive fixed-name markers/results before a run and give current-run failure evidence precedence. |
| R0-05b | Stale success plus stale SCC marker; same test `[True]` | Same isolation; only newly produced, mutually consistent current-run output can succeed. |
| R0-06a | Cache without verification metadata; `test_generic_cache_needs_valid_convergence_and_request_identity[missing]` | Require an atomic cache record tied to the current request and parsed result. |
| R0-06b | Cache from a different method; same test `[wrong_method]` | Fingerprint backend, method, settings, electronic state, input symbols/order, and initial coordinates. |
| R0-06c | Cache atom identity mismatch; same test `[wrong_symbols]` | Validate exact atom identity/order and geometry before reuse. |
| R0-07 | Provisional geomeTRIC result promoted on reuse; `test_geometric_provisional_cache_is_not_promoted` | Preserve explicit cycle-exceeded status through cache metadata and return it as non-converged. |
| R0-08a | Nonfinite geomeTRIC energy; `test_geometric_invalid_result_is_not_success[nan_energy]` | Reject nonfinite physical energy. |
| R0-08b | State/result geometry mismatch; same test `[frame_mismatch]` | Require matching symbols, shape, finite coordinates, and state/result geometry agreement. |
| R0-08c | Missing explicit convergence evidence; same test `[no_convergence_evidence]` | Require recognized optimizer convergence evidence; successful process exit alone is insufficient. |
| R0-08d | Arbitrary stale XYZ fallback; same test `[fallback_xyz]` | Read only named geomeTRIC result artifacts; validate finite forces as well as energy and geometry. |
| R0-09 | Backend exception destroys original in-memory molecule; `test_backend_failure_does_not_destroy_original_input` | Snapshot and restore coordinates, energy, and optimized coordinates on typed backend failure. |
| R0-10 | Malformed cached XYZ exits Python with `SystemExit`; `test_truncated_cache_raises_backend_error` | Parse through the typed XYZ parser and raise `BackendExecutionError`, preserving the process and diagnostics. |

The suite also tests newly discovered edge cases: a fresh ORCA frame with wrong
symbols, cache reuse for an identical request, cache rejection after a method
change, nonfinite geomeTRIC forces, and a geomeTRIC release candidate that stays
usable without being classified as converged.

## Cache and status behavior

Optimization cache records use schema 1 and include a request fingerprint, a
digest of the parsed result, and an explicit status. Legacy result-only caches
are unverified and are archived under `stale_cache/` before recomputation.
Changed input geometry, atom order, charge/multiplicity, backend/method, or
settings invalidates reuse. For geomeTRIC, the state JSON is also checked for
matching geometry, finite energy and forces, and status. A cycle-exceeded or
release-candidate result remains distinguishable from convergence. Valid
backend results remain usable if cache bookkeeping itself cannot be written;
the result is then recomputed on a later call without a verified cache record.

Fixed-name xTB outputs and markers are retained under `stale_attempts/` before
each invocation. ORCA's old result frame is retained as `.input.xyz` while its
canonical result filename is cleared, so stale coordinates cannot authorize a
new calculation.

## Validation and limitations

The tests prove software behavior at mocked process and file boundaries. They
do not qualify scientific reference energies, convergence tolerances, reaction
pathways, or backend version compatibility. Real xTB smoke evidence, the
complete repository test suite, and CI results for the opened remediation PR
are recorded with the PR. Proprietary ORCA was not installed or executed as
part of public validation. Cache schema 1 intentionally rejects legacy
unproven results and recomputes them; the legacy files are retained for
inspection rather than deleted.
