Quickstart
==========

PyAR helps you search chemically meaningful structures:

* low-energy aggregates and clusters from fragments or formulas
* AFIR-style reaction candidates from two reactants
* microsolvation shells around a solute

Install
-------

.. code-block:: bash

   python -m pip install pyar-chem

Check the CLI
-------------

.. code-block:: bash

   pyar-cli --help

Choose a chemistry task
-----------------------

Aggregate two fragments or start from a formula:

.. code-block:: bash

   pyar-cli -a C H -as 1 4 -N 8
   pyar-cli --aggregate --formula C5H4 -N 8

Run a reaction search between two reactants:

.. code-block:: bash

   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

Run a solvation search:

.. code-block:: bash

   pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16
   pyar-cli -s solute.xyz solvent.xyz --software xtb -ss 10 -N 16

For aggregate and solvation workflows, ``--connectivity-policy`` controls
whether PyAR prefers covalent connectivity when choosing survivors. Use
``auto`` for the default chemistry-aware choice, ``off`` to never filter by
connectivity, and ``prefer`` or ``strict`` when you are building atomic or
formula-driven clusters.

See also
--------

* :doc:`aggregate`
* :doc:`react`
* :doc:`solvate`
* :doc:`installation`
* :doc:`external_programs`
* :doc:`publications`
