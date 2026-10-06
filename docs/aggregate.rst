Aggregate
=========

Modern interface
----------------

.. code-block:: bash

   pyar aggregate A.xyz B.xyz
   pyar aggregate water.xyz --size 6 --backend xtb
   pyar aggregate water.xyz ammonia.xyz --size 4 2 --backend xtb
   pyar aggregate --formula CH4 --check

Multiple inputs default to one of each. A single input requires ``--size``;
explicit sizes contain one positive count per input. Inputs describe building
block types for a final composition, with existing alternative build pathways;
no component is a permanently privileged seed.

Backend refinement is optional. Geometry-only search uses the existing bounded
structural selection without artificial energies. Orientations default to eight,
and the survivor budget remains eight. Selection/connectivity/pathway overrides
retain their established policies. Formula input specifies composition, not bonds.

``--check`` validates inputs, electronic states, backend requirements when
selected, and the existing ``aggregates/state.json`` restart request without
creating output or modifying state. Legacy aggregate syntax remains supported.

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

Aggregate versus grow
---------------------

Aggregate searches the requested final composition without a permanently
privileged seed and may explore alternative build pathways. For a specified
seed with one repeatedly added species and retained intermediate stages, use
:doc:`grow`. Both use the same single-addition service.
