Workflows
=========

This page describes the main chemistry workflows exposed by ``pyar-cli``.
Use it to choose the right command for the kind of structure search you want
to run.

Aggregation
-----------

Aggregation searches for low-energy packings of one or more fragments. This
is the workflow to use for molecular clusters, noncovalent complexes, and
small aggregate models.

Examples:

.. code-block:: bash

   pyar-cli aggregate C H -as 1 4 -N 8
   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16
   pyar-cli scan-bond 1 2 A.xyz B.xyz -N 8
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

For a chemistry researcher, the main outputs to inspect are:

* ``aggregates/state.json`` for restart and provenance
* ``selected/`` for the chosen low-energy candidates
* the energy table output for quick ranking of the structures

Reaction
--------

Reaction searches operate on exactly two input structures and are meant for
reaction discovery, bond formation, or close-contact pathway exploration.

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

For ``xtb`` and ``aimnet_2`` reaction runs, PyAR uses geomeTRIC/TRIC to
optimize the backend energy/force objective plus the AFIR bias. Ordinary
aggregation and standalone optimization continue to use each backend's native
optimizer. A bonded reaction candidate is relaxed again with ``gamma=0.0``
before product identity is assessed.

For a chemist, the reaction workflow is useful when you want to:

* search for candidate products without hand-building every starting guess
* compare multiple orientations of the same reactants
* inspect whether a close-contact structure relaxes back to starting material
  or becomes a new product
* review the trace summary for energetic trends and bond changes

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

Useful files to inspect are:

* ``reaction/state.json`` for restart state
* ``reaction_trace/trace.jsonl`` for the raw evaluation history
* ``path_summary.csv`` for a compact energy/bond summary
* ``candidate_ts/`` for candidate geometries of interest

Legacy ``jobs.pkl`` reaction checkpoints are imported once when their gamma
schedule is unambiguous. A legacy checkpoint whose formatted keys have lost
distinct fractional gamma values exits with a clear error instead of resuming
an uncertain calculation.

Solvation
---------

``solvate`` is the command name, but the workflow is broader than solvent
placement. It explores how a central core grows when units are added around
it. Microsolvation is a common use case, and so is adding ligands to a
transition metal center to build an organometallic complex.

.. code-block:: bash

   pyar-cli -s solute.xyz solvent.xyz --software xtb -ss 10 -N 16
   pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16

Solvation restart state is stored as readable JSON:

.. code-block:: text

   solvation/
     state.json
     state/
       geometries/
   aggregate_002/
   aggregate_003/

``state.json`` records the input seeds, added fragment, calculation settings,
next cycle, completed cycles, and the current seeds to continue from.
Re-running an interrupted solvation with the same request resumes from the
last completed cycle and reuses the stored seed geometries.

For a chemistry researcher, the main outputs to inspect are:

* ``solvation/state.json`` for restart and cycle progress
* ``solvation/state/geometries/`` for saved seed structures
* the selected seed geometries from the final cycle

Bond scan
---------

Bond scanning evaluates a distance scan between two fragments. It is a simple
way to probe whether a bond-forming or bond-breaking coordinate behaves as
expected before committing to a more expensive reaction search.

.. code-block:: bash

   pyar-cli --scan-bond 1 2 A.xyz B.xyz -N 8
   pyar-cli scan-bond 1 2 A.xyz B.xyz -N 8
