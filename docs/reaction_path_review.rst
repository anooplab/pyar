:orphan:

Adaptive bias and reaction-path review
======================================

Review date: 2026-09-25. Scope: adaptive AFIR/soft-min forces, controller updates,
release and product gates, trajectory analysis, NEB, TS optimization, frequencies,
IRC, and final endpoint optimization. Adaptive findings below are review findings;
the implementation changes in this review address the reaction-path driver and
the requested final endpoint verification.

Outstanding adaptive findings
-----------------------------

High priority: cycle-limit recovery is not necessarily the last accepted geometry.
``pyar/backends/adaptive_geometric.py:130`` writes ``optimizer.progress`` after
``GeomOptNotConvergedError``. geomeTRIC appends the evaluated trial before its
acceptance decision, and iteration-limit termination precedes rejection handling.
The adaptive driver returns immediately on ``FAILED``, without updating its
accepted checkpoint. A one-cycle harmonic-pair reproduction exported coordinates
that differed from that checkpoint by 0.002645886 angstrom. Recovery should use
the accepted checkpoint and reevaluate its physical energy before writing the
parent-process result. Do not call that trial an accepted geometry.

High priority: the initial reactant frame is absent from adaptive traces.
``PyarGeometricCalculator.calculate`` records adaptive frames only when
``_record_accepted_geometry`` is set after a step. The first evaluation initializes
the controller but writes no trace frame. ``_persistent_transition_index`` in
``pyar/reaction_analysis.py`` consequently uses the first accepted step as its
reference. A topology change on that first step can be missed, and there is no
pre-change geometry to export. Record an explicit initial reference without
advancing release persistence; preserve this convention during restart.

Medium priority: optional force filtering uses the biased total force.
``max_force`` in ``pyar/backends/geometric.py`` is derived from the physical plus
bias force, and the candidate filter in ``pyar/reaction_analysis.py`` uses that
field. Force cancellation can therefore let a physically strained geometry pass
a low-force filter. Preserve the total-force diagnostic but expose a separate
physical-force criterion for assessing candidate strain. A TS guess itself need
not be stationary; the distinction belongs in explicit selection criteria.

Low priority: adaptive smoothing is currently ineffective.
``BiasController.select`` makes ``target`` at least ``previous_alpha + margin``
and takes the maximum of that target and its convex interpolation with the
previous alpha. For every allowed smoothing value, the maximum is the target.
Reproductions with smoothing values 0.1, 0.5, and 1.0 give identical sequences.
The CLI should describe this option as inactive/deprecated for monotonic loading,
or a scientifically explicit smoothing rule should replace it.

Additional termination detail: after exhausting ``adaptive_max_segments``, the
driver has already selected the next segment's alpha and written its checkpoint,
although that segment is never optimized. The controller checkpoint and last
calculator-state metadata can consequently describe different strengths. Final
state serialization should correspond to the last actually evaluated segment.

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
