Aggregation and Cluster Search
==============================

Use aggregation when you want PyAR to build and screen low-energy structures
from one or more molecular fragments. This is the main task for noncovalent
complexes, molecular clusters, weakly bound assemblies, and formula-based
structure generation.

Typical chemistry questions
---------------------------

Aggregation is useful when you want to ask questions such as:

* What are plausible low-energy structures of a molecular dimer or cluster?
* How can two or more fragments pack through hydrogen bonding, dispersion, or
  ion-pairing interactions?
* Which structures should be selected for a later xTB, ORCA, Gaussian, or ML
  refinement?
* Can I generate candidate structures from a formula before doing expensive
  calculations?

Basic commands
--------------

Aggregate two one-atom fragments:

.. code-block:: bash

   pyar-cli aggregate C H -as 1 4 -N 8
   pyar-cli -a C H -as 1 4 -N 8

Generate trial structures directly from a formula:

.. code-block:: bash

   pyar-cli --aggregate --formula C5H4 -N 8

Run a fragment-cluster search with a backend optimizer:

.. code-block:: bash

   pyar-cli -s water.xyz water.xyz --software xtb -ss 10 -N 16 -c 0 0 -m 1 1

How to think about the output
-----------------------------

For most chemistry users, the important outputs are the selected structures
and their energies. PyAR removes near-duplicates and keeps a smaller set of
candidate geometries for inspection or higher-level refinement.

Disconnected covalent graphs are expected for noncovalent aggregates and
solvation-like assemblies, so they are not treated as invalid here. A
connectivity preference is mainly useful for atomic or formula-driven growth
where the goal is to keep a single-component cluster together.
Open-shell or radical inputs also tend to keep that preference on, because a
simple disconnectedness test is not enough to describe their chemistry.

A typical aggregation run creates a directory structure like:

.. code-block:: text

   aggregates/
     state.json
     ag_.../
       selected/
     selected/
       stoichiometry_.../

Useful files to inspect:

* ``aggregates/state.json`` records the request, restart state, and provenance.
* ``selected/`` contains the selected candidate structures.
* Energy-table output helps rank structures by relative energy.

Restart behaviour
-----------------

Aggregation restart state is stored as readable JSON. Re-running an
interrupted aggregation with the same request resumes unfinished pathways while
reusing existing step outputs. Older ``pyar.log`` pathway markers are imported
once into JSON state when a legacy ``aggregates/`` calculation is resumed.

Next steps
----------

After an aggregation run, common follow-up steps are:

* run ``pyar-energy-table`` on selected structures
* cluster or deduplicate the final candidates
* refine selected structures with a higher-level backend
* use selected aggregates as input for reaction, solvation, or external DFT
  calculations
