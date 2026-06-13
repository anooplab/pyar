Conformer Search
================

Use conformer search when you want PyAR to build and rank a set of plausible
starting geometries for a single molecule before optional backend refinement.
The workflow is RDKit-based and is aimed at flexible molecules where a single
local minimum is not enough to describe the search space.

Basic commands
--------------

.. code-block:: bash

   pyar-conformer CCO
   pyar-cli conformer CCO
   pyar-conformer input.sdf --software xtb

What it does
------------

* embeds multiple RDKit conformers with ETKDGv3
* optionally starts from multiple seeds for broader coverage
* can keep a larger backend-refinement pool than the final output set
* can use random coordinates or a tighter prune threshold for more breadth
* optionally refines the selected pool with an existing PyAR backend

Useful options
--------------

* ``--num-conformers`` controls how many conformers are embedded per seed
* ``--num-seeds`` repeats embedding with deterministic seed offsets
* ``--backend-top-n`` widens the set sent to backend refinement
* ``--diversity-fraction`` controls how much of that backend pool comes from
  geometric diversity instead of pure RDKit energy ranking
* ``--use-random-coords`` starts RDKit from random coordinates instead of the
  default distance-geometry eigenvector start
* ``--rms-threshold`` or ``--prune-rms-threshold`` sets RDKit's greedy
  pruning threshold during embedding

Useful outputs
--------------

* ``conformers/state.json`` for restart/provenance-style metadata
* ``conformers/summary.csv`` for the ranked conformer table
* ``conformers/rdkit/`` for the embedded conformers before backend refinement
* ``conformers/selected/`` for the final selected conformers

See also
--------

* :doc:`quickstart`
* :doc:`usage`
* :doc:`workflows`
