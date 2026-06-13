Usage
=====

The main command-line entry point is ``pyar-cli``.

Examples:

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

   pyar-cli -a C H -as 1 4 -N 8

   pyar-cli --aggregate --formula C5H4 -N 8

   pyar-conformer input.sdf --software xtb --num-seeds 3 --backend-top-n 50

The repository README contains additional examples for clustering,
aggregation, solvation, and reaction searches.

Energy tables
-------------

To print a relative-energy table from raw ``.xyz`` files, use:

.. code-block:: bash

   pyar-energy-table *.xyz

The command prints absolute energies, relative energies in kcal/mol, and the
global minimum directly to the terminal.

Reaction trace analysis
-----------------------

To summarize a recorded AFIR trace and optionally generate PNG plots:

.. code-block:: bash

   pyar-reaction-trace .
   pyar-reaction-trace . --plot
   pyar-reaction-trace . --plot-only
   pyar-cli trace . --plot

The command writes ``path_summary.csv`` and ``candidate_ts/`` in the job
directory unless ``--plot-only`` is used, and places plots in
``trace_plots/`` unless ``--plot-directory`` is set.

The summary distinguishes the physical backend energy from the AFIR-biased
optimization objective:

* ``backend_energy_hartree``: backend electronic, ML, or xTB energy without AFIR
* ``afir_energy_hartree``: artificial AFIR contribution
* ``total_energy_hartree``: optimization objective used by geomeTRIC
* ``backend_relative_kcalmol``: backend-energy change relative to the first
  recorded trace frame

The candidate file ``candidate_ts/highest_backend_energy.xyz`` is usually the
first structure to inspect for future NEB, string, dimer, or TS workflows. It
is not a confirmed transition state.

Utilities
---------

Several smaller helper commands are available for inspection and benchmarking:

.. code-block:: bash

   pyar-energy-table *.xyz
   pyar-clustering *.xyz -a maxmin -n 8
   pyar-similarity -f "*.xyz" -t 0.005
   pyar-descriptor *.xyz
   pyar-conformer input.sdf --num-conformers 200 --num-seeds 4 --use-random-coords
   pyar-trial-generation -N 8 -i seed.xyz monomer.xyz --plot
   pyar-optimiser structure.xyz -c 0 -m 1 --software xtb

``pyar-energy-table`` prints relative energies, ``pyar-clustering`` selects a
diverse subset, ``pyar-conformer`` generates and optionally refines RDKit
conformers, ``pyar-similarity`` reports near-duplicate structures,
``pyar-descriptor`` writes compact cluster descriptors,
``pyar-trial-generation`` builds candidate orientations, and
``pyar-optimiser`` runs the standalone geometry optimizer.
