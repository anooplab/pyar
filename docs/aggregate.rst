Aggregate
=========

Use aggregation when you want PyAR to build and screen low-energy structures
from fragments or a formula. This is the main workflow for clusters,
noncovalent complexes, and other build-up problems.

Disconnected covalent graphs are expected for molecular aggregates and other
noncovalent complexes. Use ``--connectivity-policy auto`` for the default
chemistry-aware choice, or override with ``off``, ``prefer``, or ``strict``
when you need a different selection rule.

Basic commands
--------------

.. code-block:: bash

   pyar-cli aggregate C H -as 1 4 -N 8
   pyar-cli -a C H -as 1 4 -N 8
   pyar-cli --aggregate --formula C5H4 -N 8

What it does
------------

* generates trial geometries
* evaluates and ranks candidates
* removes near-duplicates
* persists restart state in ``aggregates/state.json``

Useful outputs
--------------

* ``aggregates/state.json`` for restart and provenance
* ``selected/`` for the chosen structures
* the energy table for quick inspection of relative energies

See also
--------

* :doc:`quickstart`
* :doc:`installation`
* :doc:`workflows`
