Numerical tolerance policy
==========================

Every numerical benchmark case must state the metric, unit, expected value,
and nonnegative tolerance in its frozen manifest. There is no universal
cross-method tolerance. Keep absolute and relative tolerances separate when
near-zero values make relative error unstable. Record finite-difference step
size and norm/component comparisons for derivative checks.

The policy classes are exact invariants, analytical/numerical derivatives,
same-method numerical replay, structural equivalence, reaction-path evidence,
and stochastic search outcomes. Their required metadata and comparison
interpretation are specified in :doc:`index` under “Tolerance policy”.
Release thresholds remain unset until the scientific owner approves the
corpus, precision, thread policy, and backend versions.
