xTB Guide
=========

This page summarizes how PyAR uses xTB at runtime.

Wrappers
--------

PyAR uses xTB in four places:

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
composed xTB plus AFIR objective. See
:doc:`reaction_optimization`.

Shared behavior
---------------

* ``--nprocs`` is passed through to xTB as ``--parallel``.
* Charge, multiplicity, and spin handling are normalized in
  ``pyar.backends.xtb_utils``.
* The wrappers fail fast if ``xtb`` is not available on ``PATH``.

See the full API documentation in :doc:`api`.
