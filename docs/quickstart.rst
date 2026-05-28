Quickstart
==========

PyAR helps you search chemically meaningful structures:

* low-energy aggregates and clusters from fragments or formulas
* AFIR-style reaction candidates from two reactants
* microsolvation shells around a solute
* bond scans between two fragments

If you want to see published chemistry problems that used these tasks, see
:doc:`publications`.

.. seealso::

   :doc:`publications` for a short list of papers that used PyAR in
   aggregation, reaction discovery, solvation-style growth, catalyst
   formation, and chemical-space exploration.

Install
-------

.. code-block:: bash

   python -m pip install .

For development:

.. code-block:: bash

   python -m pip install -e .

Check the CLI
-------------

.. code-block:: bash

   pyar-cli --help

Choose a chemistry task
-----------------------

Aggregate two one-atom fragments when you want to find low-energy dimer or
cluster geometries:

.. code-block:: bash

   pyar-cli aggregate C H -as 1 4 -N 8
   pyar-cli -a C H -as 1 4 -N 8

Generate a starting aggregate directly from a formula when you already know
the composition:

.. code-block:: bash

   pyar-cli --aggregate --formula C5H4 -N 8

Run a reaction search when you want to explore bond formation or close-contact
reaction pathways between two reactants:

.. code-block:: bash

   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

Scan a bond when you want a simple reaction-coordinate probe:

.. code-block:: bash

   pyar-cli scan-bond 1 2 A.xyz B.xyz -N 8
   pyar-cli --scan-bond 1 2 A.xyz B.xyz -N 8

Run a solvation search when you want to add units around a central core.
Microsolvation is the classic case, but the same task is also useful for
building coordination complexes, such as adding ligands to a transition metal
center:

.. code-block:: bash

   pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16
   pyar-cli -s solute.xyz solvent.xyz --software xtb -ss 10 -N 16

Notes
-----

* In aggregate mode, bare element symbols such as ``C`` and ``H`` are accepted.
* If ``--software`` is omitted for aggregate mode, PyAR generates trial
  geometries only and skips quantum-chemistry optimization.
* ``--formula`` is only valid together with ``--aggregate``.
* Reaction runs with registered energy-gradient providers such as ``xtb`` and
  ``aimnet_2`` use geomeTRIC/TRIC for AFIR-biased optimization, then use
  ``gamma=0.0`` for product relaxation.
