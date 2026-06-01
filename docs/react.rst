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
