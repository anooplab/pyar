React
=====

Use reaction search when you want PyAR to explore possible bond formation,
bond rearrangement, or close-contact reaction candidates between two input
reactants.

Basic commands
--------------

.. code-block:: bash

   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

The reaction workflow uses AFIR-style biased optimisation and then checks
whether the relaxed structure is a new product.

Adaptive bias control
---------------------

For registered Cartesian energy-gradient backends, enable the accepted-step
adaptive controller with ``--bias-controller adaptive``. For example, test an
HCN + HCN orientation with xTB:

.. code-block:: bash

   pyar-cli react HCN_A.xyz HCN_B.xyz -N 1 --bias-max 100 \
     --software xtb --bias-controller adaptive

The maximum force scale is set by ``--bias-max``. The adaptive controller
holds alpha fixed within each optimization segment and updates it after a
segment converges. Optional
``--bias-alpha-min``, ``--bias-alpha-margin``, and ``--bias-alpha-epsilon``
flags control its lower bound, safety margin, and denominator regularization.
Alpha values and the safety margin are in Hartree/Bohr.
The margin defaults to 0.001 and must be positive: a zero
margin merely cancels physical resistance and can stall before contact.
Each segment applies the instantaneous force-cancellation estimate plus the
positive margin, with a monotonic increment floor. Alpha smoothing is
incompatible with that loading rule; the deprecated ``--bias-alpha-smoothing``
option accepts only ``1``.
Adaptive runs use a single ceiling from ``--bias-max``; ``--bias-min`` is ignored.
Convergence below the ceiling starts another segment with increased loading,
until release criteria or an iteration/segment limit stops the search.
For a constant scheduled scale, use ``--bias-controller scheduled
--bias-scheduled-alpha VALUE``.
Controller settings and per-step decisions are written to the reaction trace
and ``path_summary.csv``. ``pyar-react`` accepts the same flags.
CSV ``bias_alpha_critical`` and ``bias_alpha_target`` describe the current
accepted geometry; ``bias_segment_alpha_critical`` and
``bias_segment_alpha_target`` retain the estimates used to start its segment.
``bias_alpha`` is the applied strength.
For ``pyar-cli`` geomeTRIC runs, ``--opt-cycles`` and ``--opt-threshold``
control the external optimizer even if the energy-gradient backend does not
support these options for its native optimizer.

Products are accepted only when their InChI differs from the separated
reactants after unbiased relaxation. SMILES are canonicalized for reporting.
Bond orders omitted by backend print thresholds remain unknown. Incomplete
electronic history permits a geometry-based free-relaxation probe once contacts
are persistent and stable; complete electronic history must show growth and
stabilization. This fallback does not itself confirm a product.
Candidate selection and result XYZ energy labels use physical backend
energies. Biased energies and continuity offsets remain in calculator state
and traces. Older reaction states using biased energies cannot be resumed;
start a fresh calculation in a separate directory.

Supported AFIR energy-gradient providers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``xtb``
* ``aimnet_2``
* ``orca``
* ``gaussian``

Useful outputs
--------------

* ``reaction/state.json`` for restart and provenance
* ``reaction_trace/trace.jsonl`` and ``reaction_trace/steps/`` for trace data
* ``path_summary.csv`` for a compact path summary
* ``candidate_ts/`` for candidate geometries to inspect further

Trace analysis
~~~~~~~~~~~~~~

.. code-block:: bash

   pyar-reaction-trace reaction/gamma_0100/orientation_xxxxx --plot
   pyar-reaction-trace . --plot

See also
--------

* :doc:`quickstart`
* :doc:`reaction_optimization`
* :doc:`workflows`
