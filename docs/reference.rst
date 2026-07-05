Reference
=========

This page collects the main modules that make up PyAR.

Core types
----------

Canonical core, I/O, and sampling APIs are documented in :doc:`api`.

Main workflows
--------------

.. automodule:: pyar.scripts.explore
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.optimiser
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.react
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.reaction_trace
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.energy_table
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.trial_generation
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.conformer
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.clustering
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.similarity
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.descriptor
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.benchmark_clustering
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.conformer_benchmark
   :members:
   :undoc-members:

.. automodule:: pyar.benchmarks.conformer
   :members:
   :undoc-members:

.. automodule:: pyar.scripts.benchmark_orientations
   :members:
   :undoc-members:

.. automodule:: pyar.workflow_results
   :members:
   :undoc-members:

.. automodule:: pyar.backends
   :members:
   :undoc-members:

.. automodule:: pyar.backends.orca
   :members:
   :undoc-members:

.. automodule:: pyar.backends.orca_aiqm1
   :members:
   :undoc-members:

.. automodule:: pyar.backends.psi4
   :members:
   :undoc-members:

.. automodule:: pyar.backends.gaussian
   :members:
   :undoc-members:

.. automodule:: pyar.backends.babel
   :members:
   :undoc-members:

.. automodule:: pyar.backends.mopac
   :members:
   :undoc-members:

.. automodule:: pyar.backends.ani
   :members:
   :undoc-members:

.. automodule:: pyar.backends.aimnet_2
   :members:
   :undoc-members:

.. automodule:: pyar.backends.aiqm1_mlatom
   :members:
   :undoc-members:

.. automodule:: pyar.backends.mlatom_aiqm1
   :members:
   :undoc-members:

.. automodule:: pyar.backends.mlopt
   :members:
   :undoc-members:

.. automodule:: pyar.state.aggregate
   :members:
   :undoc-members:

.. automodule:: pyar.state.reaction
   :members:
   :undoc-members:

.. automodule:: pyar.state.solvation
   :members:
   :undoc-members:

Feature request packages
------------------------

These packages own validated request models and the persisted restart-contract
payloads for their workflows. The workflow modules below remain the public
orchestration entry points.

.. automodule:: pyar.aggregation
   :no-index:
   :members:
   :undoc-members:

.. automodule:: pyar.aggregation.request
   :no-index:
   :members:
   :undoc-members:

.. automodule:: pyar.solvation
   :no-index:
   :members:
   :undoc-members:

.. automodule:: pyar.solvation.request
   :no-index:
   :members:
   :undoc-members:

.. automodule:: pyar.conformer
   :no-index:
   :members:
   :undoc-members:

.. automodule:: pyar.conformer.request
   :no-index:
   :members:
   :undoc-members:

Workflow orchestration
----------------------

.. automodule:: pyar.workflows.aggregate
   :members:
   :undoc-members:

.. automodule:: pyar.workflows.conformer
   :members:
   :undoc-members:

.. automodule:: pyar.reaction_analysis
   :members:
   :undoc-members:

.. automodule:: pyar.reaction_trace
   :members:
   :undoc-members:

.. automodule:: pyar.workflows.reaction
   :members:
   :undoc-members:

.. automodule:: pyar.workflows.solvation
   :members:
   :undoc-members:

.. automodule:: pyar.solvation_state
   :members:
   :undoc-members:

.. automodule:: pyar.optimiser
   :members:
   :undoc-members:

.. automodule:: pyar.scan
   :members:
   :undoc-members:

Supporting utilities
--------------------

.. automodule:: pyar.property
   :members:
   :undoc-members:

Selection services
------------------

.. automodule:: pyar.selection.basin_memory
   :members:
   :undoc-members:

.. automodule:: pyar.selection.deduplication
   :members:
   :undoc-members:

.. automodule:: pyar.selection.diversity
   :members:
   :undoc-members:

.. automodule:: pyar.selection.clustering
   :members:
   :undoc-members:

For xTB-specific APIs, see :doc:`api`.
