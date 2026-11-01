Supported Python interfaces
===========================

Use canonical scientific modules rather than retired compatibility imports:

.. code-block:: python

   from pyar.core.molecule import Molecule
   from pyar.workflows.aggregate import aggregate
   from pyar.workflows.grow import grow
   from pyar.workflows.microsolvation import microsolvate
   from pyar.workflows.scan_bond import run_scan_bond
   from pyar.neb import run_neb

Consult the signatures/docstrings below before constructing requests. Python
calls must supply the same scientific state/provenance as CLI requests. A
low-level optimizer call is not a substitute for a workflow's request/preflight,
selection or scientific acceptance. Private functions may change and are not a
persistence compatibility promise.

.. automodule:: pyar.core.molecule
   :members: Molecule

.. automodule:: pyar.optimiser
   :members: optimise, bulk_optimize

See :doc:`/generated_api` for implementation details. Backend APIs and stored
protocol versions must be tested independently; import success is not a real
backend calculation.
