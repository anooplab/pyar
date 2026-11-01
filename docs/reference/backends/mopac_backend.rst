MOPAC Backend
=============

MOPAC is a semiempirical quantum-chemistry backend available in PyAR through
native optimization routes.

When to use MOPAC
-----------------

Use ``mopac`` when:

* MOPAC is installed and available on your system
* you specifically want a MOPAC semiempirical method
* your scientific protocol specifies a supported MOPAC method

Example command
---------------

.. code-block:: bash

   pyar optimize structure.xyz -c 0 -m 1 --backend mopac

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
