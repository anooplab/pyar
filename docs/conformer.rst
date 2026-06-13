Conformer Search
================

Use conformer search when you want PyAR to build and rank a set of plausible
starting geometries for a single molecule before optional backend refinement.
The workflow is RDKit-based and is aimed at flexible molecules where a single
local minimum is not enough to describe the search space.

The default settings are tuned for folded, compact conformers: multiple RDKit
seeds, random-coordinate embedding, chemistry-guided evolutionary torsion
search, tighter pruning, a compactness-biased backend pool, and a larger
backend-refinement window than the final output set.

Basic commands
--------------

.. code-block:: bash

   pyar-conformer CCO
   pyar-cli conformer CCO
   pyar-conformer input.sdf --software xtb

What it does
------------

* embeds multiple RDKit conformers with ETKDGv3
* perturbs rotatable-bond torsions around embedded conformers
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
* ``--compactness-fraction`` reserves part of the backend pool for compact,
  contact-rich conformers before the diversity fill
* ``--use-random-coords`` starts RDKit from random coordinates instead of the
  default distance-geometry eigenvector start
* ``--torsion-kicks`` or ``--no-torsion-kicks`` controls the local
  torsion-perturbation stage
* ``--torsion-mode`` chooses the torsion sampler; ``evolve`` is the default,
  ``mc`` keeps the tabu-style walk behavior, ``grid`` keeps the deterministic
  scan behavior, and ``random`` keeps the stochastic kick behavior
* ``--torsion-rounds`` controls how many successive torsion-kick rounds are run
* ``--torsion-kicks-per-conformer`` controls how many torsion trials are made
  around each embedded conformer
* ``--torsion-dedup-rms`` removes near-duplicate trial conformers before
  backend selection
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
