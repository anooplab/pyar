:orphan:

Adaptive bias and reaction-path review
======================================

Review date: 2026-09-25. Scope: adaptive AFIR/soft-min forces, controller updates,
release and product gates, trajectory analysis, NEB, TS optimization, frequencies,
IRC, and final endpoint optimization. Adaptive findings below are review findings;
the implementation changes in this review address the reaction-path driver and
the requested final endpoint verification.

Adaptive findings addressed
---------------------------

High priority, fixed: cycle-limit recovery now restores the accepted controller
checkpoint, reevaluates that geometry, and replaces the terminal progress frame
before writing the optimizer result. A rejected final trial is not exported.

High priority, fixed: the initial adaptive evaluation is recorded as the explicit
reactant topology baseline, without incrementing release persistence. Restarted
traces continue from their existing records rather than inserting a second input
baseline.

Medium priority, fixed: ``max_force`` remains the total-objective diagnostic;
candidate strain filtering now uses the separately recorded maximum per-atom
backend force. Existing traces without that field use backend force vectors;
frames lacking physical force data cannot pass an enabled force filter.
Outlier rejection uses unbiased ``backend_energy_hartree``;
zero-MAD traces now classify energies distinct from the median as outliers rather
than silently disabling the filter. A TS guess itself need not be stationary.

Medium priority, fixed: the adaptive controller now rejects non-default alpha
smoothing. The loading rule is the instantaneous force-cancellation estimate
plus a positive margin, with a monotonic increment floor; smoothing that estimate
would undermine the stated force-cancellation guarantee. The CLI retains the
option for compatibility but documents that only ``1`` is supported.

Additional termination detail, fixed: the driver checks the segment budget before
advancing the controller, so checkpoint strength always describes a segment that
was actually evaluated.

Follow-up corrections (2026-09-26): missing early bond-order values no longer
permanently prevent geometry-based release probes. Endpoint validation explicitly
requires distinct observed minima, and stage dependencies bind the upstream
validation summary as well as numerical artifacts. Trace CSVs distinguish current
resistance estimates from segment-start estimates. Filtered topology events are
reported unavailable when no eligible event geometry remains. The restart
regression test now uses genuinely different controller settings.

Scientific interpretation
-------------------------

The reviewed force signs and AFIR/soft-min analytical gradients pass the existing
finite-difference tests. Trial evaluations retain the segment's alpha; loading
is reestimated after optimization of a biased segment. The additive offset
preserves energy continuity at a strength change. Use ``backend_energy_hartree``
for physical energy comparisons: the offset-containing total objective is not
a barrier profile.

The projected cancellation strength is a local force coefficient in Hartree/Bohr,
not a barrier energy. It cancels resistance along the selected contact coordinate;
it does not ensure that this coordinate follows the desired chemical channel.
Soft-min can concentrate on an unproductive closest contact, and monotonic loading
can continue compressing an already unfavorable contact. The release gate and
unbiased product-identity check remain necessary. Topology-change frames are TS
guesses whose relevance must be established by stationary-point and IRC validation.

Reaction-path corrections implemented
-------------------------------------

* Removed undefined endpoint variables left by the previous NEB stage refactor.
* Replaced the coupled TS/IRC helper with independent ``ts``, ``frequency``,
  ``irc``, and ``endpoints`` stages. Explicit TS guesses need no endpoint files.
* Bound stage artifacts to geometry hashes, physical settings, and validated
  dependencies. Fresh scratch avoids reusing an unrelated Hessian. Failed
  reruns cannot reuse an old successful stage summary.
* Required a converged band and interior maximum before automatic TS optimization.
  Added finite XYZ and iteration/tolerance checks before calculation.
* Required stationarity as well as one significant imaginary frequency for IRC.
  Frequencies and Hessians are calculated on the unbiased physical surface.
* Traced IRC directions independently. In geomeTRIC 1.1.1, the combined-direction
  driver can continue after the forward iteration limit and finally report
  backward convergence; that final state does not establish both branches.
* Added final geomeTRIC optimization and frequency checks for both IRC endpoints.
  Topology is allowed to change during this optimization; classification uses the
  resulting minima and requires a one-to-one match to distinct reference endpoints.
* Resolved sub-1e-5 angstrom numerical departures from linearity before evaluating
  frequency gradients and Hessians, retaining both bending modes of linear HCN.
  The correction and exact evaluated geometry are persisted.

Validation
----------

The numerical integration example is HCN to HNC using the installed xTB provider
(``xtb --gxtb``), not HCN dimerization. NEB converged in 34 cycles; the TS has one
significant imaginary frequency, -1426.1094 cm⁻¹. Both IRC branches converged.
The reoptimized endpoints have zero significant imaginary frequencies and match
HCN/HNC with mapped RMSDs of 4.56e-7 and 8.63e-7 angstrom, respectively.

Final energies are -5.504066223643 and -5.472159887922 Hartree. The TS energy is
-5.387373532811 Hartree. The final HCN/HNC analyses each retain four vibrational
modes. This validates the small-molecule execution and acceptance sequence; it
does not establish the chemical adequacy of this model for HCN dimerization or
constitute a new adaptive reaction-search benchmark.

Regression coverage exercises the complete stage sequence, separate execution,
failed prerequisites, different-method and modified-geometry rejection, failure
of either IRC branch, two branches reaching the same minimum, nonstationary
frequency inputs, and final minimum verification. Existing adaptive tests cover
analytical forces, fixed trial segments, checkpoint restart, release persistence,
and rejection of failed/dissociating release probes; those tests do not resolve
the outstanding adaptive findings above.

See :doc:`neb` for commands and artifact conventions. Methodology references:
`geomeTRIC transition-state documentation <https://geometric.readthedocs.io/en/latest/transition.html>`_
and `geomeTRIC IRC documentation <https://geometric.readthedocs.io/en/latest/irc.html>`_.
