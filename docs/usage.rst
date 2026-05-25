Usage
=====

The main command-line entry point is ``pyar-cli``.

Examples:

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

   pyar-cli -a C H -as 1 4 -N 8

   pyar-cli --aggregate --formula C5H4 -N 8

The repository README contains additional examples for clustering,
aggregation, solvation, and reaction searches.

Energy tables
-------------

To print a relative-energy table from raw ``.xyz`` files, use:

.. code-block:: bash

   pyar-energy-table *.xyz

The command prints absolute energies, relative energies in kcal/mol, and the
global minimum directly to the terminal.
