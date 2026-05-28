MOPAC Backend
=============

MOPAC is a semiempirical quantum-chemistry backend available in PyAR through
legacy optimisation routes.

When to use MOPAC
-----------------

Use ``mopac`` when:

* MOPAC is installed and available on your system
* you specifically want a MOPAC semiempirical method
* you are working with legacy PyAR inputs or comparisons

Example command
---------------

.. code-block:: bash

   pyar-optimiser structure.xyz -c 0 -m 1 --software mopac

Practical notes
---------------

* MOPAC must be installed separately.
* Check charge, multiplicity, and method support on a small molecule before
  running batch calculations.
* For new broad exploratory work, ``xtb`` is usually the simpler first choice.

AFIR/geomeTRIC status
---------------------

``mopac`` is not currently registered as an energy-gradient provider for the
geomeTRIC-backed AFIR reaction route.
