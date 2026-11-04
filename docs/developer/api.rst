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
   from pyar.workflows.scan_bond import run_scan_bond, run_scan_bond_request
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

Request-taking bond scans
-------------------------

``run_scan_bond_request(settings, recover_stale=False)`` accepts a mapping of
the existing ``run_scan_bond`` keyword fields: ``input_a``, ``input_b`` (``None``
for an intramolecular scan), ``atoms``, ``orientations``, ``qc_params`` and
``output_dir``, plus the same optional grid and continuation controls. The
mapping is copied and bound against the established signature; the CLI reuses
one invocation for preflight and execution. ``recover_stale`` is an execution
control passed separately. Backend electronic state must be resolved as for the
CLI; this API does not infer different scientific defaults.

This is an invocation mapping, not a stored signed request or a manifest.
``ScanLayout`` still creates/verifies request signatures, immutable scan and
relaxation evidence, artifact hashes, product membership and restart authority.
``ScanRun`` and per-orientation contexts identify the owners of that evidence.

.. autofunction:: pyar.workflows.scan_bond.run_scan_bond_request

Resolved reaction discovery
---------------------------

``run_reaction(request, directory=..., recover_stale=False)`` accepts the
``ReactionRequest`` returned by ``resolve_reaction_request``. The command line
uses this same entry point. The existing ``react(...)`` argument API remains
supported. Discovery owns restart verification; the request's descriptive
``restart_request`` field cannot authorize reuse of persisted artifacts.

.. autofunction:: pyar.workflows.reaction.run_reaction

Legacy graph conversion
-----------------------

``pyar.checkvalidity.guess_correct_molecules`` reads zero-based atom-indexed
weighted edge lists and writes the existing MOL2 schema. Atom labels retain
their weighted-degree element inference regardless of edge insertion order.
This historical heuristic is not a chemical-validity or valence certification.
Rejected empty, disconnected, out-of-range or excessive-degree graphs leave existing output
untouched; unexpected parser or graph implementation failures propagate.

Graph conversion and valence JSON loading do not require the historical
``pyar.sascore`` provider. ``calculate_property`` still requires that provider;
no alternative SA scoring algorithm is selected automatically.
