Quickstart
==========

PyAR helps you search chemically meaningful structures:

* low-energy aggregates and clusters from fragments or formulas
* AFIR-style reaction candidates from two reactants
* microsolvation shells around a solute
* RDKit conformer ensembles before backend refinement

Install
-------

.. code-block:: bash

   python -m pip install pyar-chem

Check the CLI
-------------

.. code-block:: bash

   pyar --help

Choose a chemistry task
-----------------------

Aggregate two fragments or start from a formula:

.. code-block:: bash

   pyar aggregate C H --size 1 4 --backend xtb
   pyar aggregate --formula C5H4 --backend xtb --check

Run a reaction search between two reactants:

.. code-block:: bash

   pyar react A.xyz B.xyz --backend xtb --bias-max 100 --through all

Run a solvation search:

.. code-block:: bash

   pyar microsolvate solute.xyz solvent.xyz --count 10 --backend xtb

Run a conformer search:

.. code-block:: bash

   pyar conformer input.sdf --backend xtb

``pyar-cli`` and ``pyar-*`` entry points remain available for existing scripts
and legacy workflows. In particular, ``pyar-cli solvate`` retains its legacy
algorithm and restart state; modern ``pyar solvate`` aliases microsolvation.

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
* :doc:`conformer`
* :doc:`installation`
* :doc:`external_programs`
* :doc:`publications`
