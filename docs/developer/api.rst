Supported Python interfaces
===========================

Use canonical scientific modules rather than retired compatibility imports:

.. code-block:: python

   from pyar.core.molecule import Molecule
   from pyar.conformer.request import ConformerRequest
   from pyar.workflows.conformer import conformer_search, run_conformer_search
   from pyar.workflows.aggregate import aggregate, run_aggregate
   from pyar.workflows.grow import grow, run_grow
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

Request-taking conformer execution
----------------------------------

``run_conformer_search(request, root_directory=...)`` accepts an already resolved
``ConformerRequest`` and executes exactly the same fresh-run workflow as
``conformer_search(...)``. The latter retains its supported keyword signature
and performs request conversion once. CLI execution reuses its resolved request.
Both interfaces retain stereo checks, seeded generation, physical backend energy
ranking, immutable artifacts, and rejection of existing nonresumable runs.

.. autofunction:: pyar.workflows.conformer.run_conformer_search

Request-taking aggregation and growth
-------------------------------------

``run_aggregate(request, recover_stale=False)`` accepts an ``AggregateRequest``.
``run_grow(request, output="grow", recover_stale=False)`` accepts a ``GrowRequest``.
The existing ``aggregate``, ``aggregate_from_formulas`` and ``grow`` interfaces
remain supported. CLI request resolution is reused directly; pathway ordering,
cumulative electronic states, survivor budgets, immutable step commits and
restart authorization remain owned by the workflow.

.. autofunction:: pyar.workflows.aggregate.run_aggregate

.. autofunction:: pyar.workflows.grow.run_grow
