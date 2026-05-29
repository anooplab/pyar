xTB Guide
=========

This page summarizes how PyAR uses xTB at runtime.

Homepage and installation
-------------------------

xTB is free/open-source software from the Grimme group.

* Homepage/source: https://github.com/grimme-lab/xtb
* Documentation and setup: https://xtb-docs.readthedocs.io/en/latest/setup.html

After installing xTB, check that it is visible to your shell:

.. code-block:: bash

   xtb --version

When to use xTB
---------------

xTB is usually the best first backend for new PyAR users. It is fast enough for
exploration and useful for aggregation, solvation, standalone optimisation, and
AFIR/geomeTRIC reaction searches.

Wrappers
--------

PyAR uses xTB in several places:

* ``pyar.backends.xtb`` for standalone xTB optimization
* ``pyar.backends.xtb_aiqm1`` for xTB followed by AIQM1 refinement
* ``pyar.backends.xtb_aimnet2`` for xTB followed by AIMNet2 refinement
* ``pyar.backends.xtb_turbo`` for the retained Turbomole/AFIR compatibility
  workflow that calls xTB for gradients
* ``pyar.backends.geometric`` for geomeTRIC/TRIC reaction optimization using
  xTB energy and force evaluations plus the AFIR bias

The ``xtb_turbo`` path requires Turbomole for coordinate updates and
convergence checking even though xTB supplies the physical gradient. Reaction
runs selected with ``--software xtb`` now use the geomeTRIC channel over a
composed xTB plus AFIR objective. See :doc:`reaction_optimization`.

Example commands
----------------

Standalone optimisation:

.. code-block:: bash

   pyar-optimiser structure.xyz -c 0 -m 1 --software xtb

Small solvation or cluster-growth job:

.. code-block:: bash

   pyar-cli solvate water.xyz water.xyz --software xtb -ss 1 -N 4 -c 0 0 -m 1 1

Reaction search:

.. code-block:: bash

   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

Shared behaviour
----------------

* ``--nprocs`` is passed through to xTB as ``--parallel``.
* Charge, multiplicity, and spin handling are normalized in
  ``pyar.backends.xtb_utils``.
* The wrappers fail fast if ``xtb`` is not available on ``PATH``.

See the full API documentation in :doc:`api`.
