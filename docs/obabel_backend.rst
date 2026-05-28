OpenBabel Backend
=================

OpenBabel is available in PyAR as a lightweight molecular preparation and
optimisation route through commands such as ``obabel``, ``babel``,
``obminimize``, and ``obenergy``.

When to use OpenBabel
---------------------

Use ``obabel`` when:

* you need a lightweight structure preparation or force-field-style route
* you want a quick geometry cleanup before more serious calculations
* charge and multiplicity handling are not central to the calculation

Example command
---------------

.. code-block:: bash

   pyar-optimiser structure.xyz --software obabel

Practical notes
---------------

* OpenBabel tools must be installed separately and available on ``PATH``.
* This route is best treated as lightweight preparation, not final chemistry.
* For reaction or mechanistic work, prefer xTB, AIMNet2, ORCA, Gaussian, or a
  validated higher-level backend.

AFIR/geomeTRIC status
---------------------

``obabel`` is not currently registered as an energy-gradient provider for the
geomeTRIC-backed AFIR reaction route.
