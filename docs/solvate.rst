Solvate
========

Use solvation when you want PyAR to add one or more copies of a fragment
around a central structure. Microsolvation is the classic use case, but the
same workflow also applies to ligand addition and local growth around a core.

Basic commands
--------------

.. code-block:: bash

   pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16
   pyar-cli -s solute.xyz solvent.xyz --software xtb -ss 10 -N 16

What it does
------------

* starts from a solute or central core
* generates trial orientations of the added fragment
* selects diverse low-energy structures
* uses selected structures as seeds for the next growth cycle

Useful outputs
--------------

* ``solvation/state.json`` for restart and provenance
* ``solvation/state/geometries/`` for saved seed structures
* the selected structures from the final cycle

See also
--------

* :doc:`quickstart`
* :doc:`workflows`
