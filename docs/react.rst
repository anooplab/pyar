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

   pyar-cli react HCN_A.xyz HCN_B.xyz -N 1 --bias-min 100 --bias-max 100 \
     --software xtb --bias-controller adaptive --bias-alpha-margin 0.001 \
     --bias-alpha-smoothing 0.5

The maximum force scale is set by ``--bias-max``. The adaptive controller
chooses a scale up to this ceiling after accepted optimizer steps. Optional
``--bias-alpha-min``, ``--bias-alpha-margin``, ``--bias-alpha-smoothing``, and
``--bias-alpha-epsilon`` flags control its lower bound, safety margin, temporal
filtering, and denominator regularization. Alpha values and the safety margin
are in Hartree/Bohr; the smoothing fraction must be between zero and one.
For a constant scheduled scale, use ``--bias-controller scheduled
--bias-scheduled-alpha VALUE``.
Controller settings and per-step decisions are written to the reaction trace
and ``path_summary.csv``. ``pyar-react`` accepts the same flags.

Product validity is determined by molecular identity and bond-change logic,
not by a simple connected-versus-disconnected connectivity test.

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
