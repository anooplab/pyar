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

Reaction
--------

Reaction searches operate on exactly two input structures.

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software orca

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
