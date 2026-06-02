Workflow Internals
==================

This page is a technical map of the main PyAR workflow implementations. If you
are using PyAR for chemistry, start with the task pages instead:
:doc:`aggregate`, :doc:`react`, and :doc:`solvate`.

The word "workflow" is used here in the developer sense: a coordinated set of
sampling, optimisation, selection, restart, and reporting steps.

Aggregation internals
---------------------

Aggregation searches for low-energy packings of one or more fragments. This is
the internal route behind molecular clusters, noncovalent complexes, and small
aggregate models.

Connectivity filtering is policy-controlled here: it is useful for atomic or
formula-driven growth, but noncovalent aggregate and solvation workflows keep
disconnected covalent graphs by default.

The CLI exposes this choice as ``--connectivity-policy {auto,off,prefer,strict}``.
Use ``auto`` for the workflow default, ``off`` for noncovalent complexes and
solvation, and ``prefer`` or ``strict`` for atomic/formula growth.

Examples:

.. code-block:: bash

   pyar-cli aggregate C H -as 1 4 -N 8
   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16
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

Reaction internals
------------------

Reaction searches operate on exactly two input structures and are meant for
reaction discovery, bond formation, or close-contact pathway exploration.

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

The geomeTRIC/TRIC reaction route is used for registered energy-gradient
providers. At present, this route is wired for ``xtb``, ``aimnet_2``,
``orca``, and ``gaussian``. In practice, ``xtb`` and ``aimnet_2`` are the
easier immediately usable options; ``orca`` and ``gaussian`` require the
corresponding executable and should be validated on the target installation.
Ordinary aggregation and standalone optimization continue to use each
backend's native optimizer. A bonded reaction candidate is relaxed again with
``gamma=0.0`` before product identity is assessed.
Reaction product validity is determined by molecular identity and bond-change
logic, not by single-component covalent connectivity.

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

.. code-block:: text

   reaction/
     state.json
     gamma_.../
       orientation_.../
         reaction_trace/
           trace.jsonl
           steps/step_*.xyz
         path_summary.csv
         candidate_ts/
           highest_backend_energy.xyz
           highest_total_energy.xyz
           pre_product_geometry.xyz
           max_bond_change.xyz
           metadata.json
         trace_plots/
           reaction_profile.png

``backend_energy_hartree`` is the physical backend energy without the AFIR
bias. ``total_energy_hartree`` is the optimization objective, including AFIR,
that geomeTRIC follows. The candidate file
``candidate_ts/highest_backend_energy.xyz`` is usually the first structure to
inspect for later NEB, string, dimer, or TS attempts. The file
``candidate_ts/pre_product_geometry.xyz`` is based on the first persistent
connectivity change and should not be treated as a confirmed transition state.

Legacy ``jobs.pkl`` reaction checkpoints are imported once when their gamma
schedule is unambiguous. A legacy checkpoint whose formatted keys have lost
distinct fractional gamma values exits with a clear error instead of resuming
an uncertain calculation.

Solvation internals
-------------------

``solvate`` is the command name, but the internal route is broader than solvent
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

Bond-scan internals
-------------------

Bond scanning evaluates a distance scan between two fragments. It is a simple
way to probe whether a bond-forming or bond-breaking coordinate behaves as
expected before committing to a more expensive reaction search.

.. code-block:: bash
