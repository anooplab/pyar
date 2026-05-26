Workflows
=========

This page describes the main execution paths exposed by ``pyar-cli``.

Aggregation
-----------

Aggregation searches for low-energy packings of one or more fragments.

Examples:

.. code-block:: bash

   pyar-cli -a C H -as 1 4 -N 8
   pyar-cli --aggregate --formula C5H4 -N 8

Aggregation restart state is stored as readable JSON:

.. code-block:: text

   aggregates/
     state.json
     ag_.../
       selected/
     selected/
       stoichiometry_.../

``state.json`` records the input geometry and calculation settings, selected
pathway order, completed pathways, pathway-level selected results, and final
selected results. Re-running an interrupted aggregation with the same request
resumes only unfinished pathways while reusing their existing step outputs.
Legacy ``pyar.log`` pathway markers are imported once into JSON state when an
older ``aggregates/`` calculation is resumed.

Reaction
--------

Reaction searches operate on exactly two input structures.

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

For ``xtb`` and ``aimnet_2`` reaction runs, PyAR uses geomeTRIC/TRIC to
optimize the backend energy/force objective plus the AFIR bias. Ordinary
aggregation and standalone optimization continue to use each backend's native
optimizer. A bonded reaction candidate is relaxed again with ``gamma=0.0``
before product identity is assessed.

Reaction restart state is stored as readable JSON and XYZ snapshots:

.. code-block:: text

   reaction/
     state.json
     state/
       geometries/
     gamma_0100/
     products/

``state.json`` records the numeric gamma schedule, current cycle, pending and
retained geometries, completed jobs, discovered products, and the calculation
settings used for restart validation. Re-running the same command in an
interrupted calculation directory resumes compatible pending work. Completed
state is retained as a run record. An existing ``reaction/`` directory without
a compatible state record is never overwritten automatically; start from a new
directory or remove archived output deliberately.

Legacy ``jobs.pkl`` reaction checkpoints are imported once when their gamma
schedule is unambiguous. A legacy checkpoint whose formatted keys have lost
distinct fractional gamma values exits with a clear error instead of resuming
an uncertain calculation.

Solvation
---------

Solvation takes a solute and a solvent fragment and explores solvent placement.

.. code-block:: bash

   pyar-cli -s solute.xyz solvent.xyz --software xtb -ss 10 -N 16

Bond scan
---------

Bond scanning evaluates a distance scan between two fragments.

.. code-block:: bash

   pyar-cli --scan-bond 1 2 A.xyz B.xyz -N 8
