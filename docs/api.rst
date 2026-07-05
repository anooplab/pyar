API Reference
=============

Core package
------------

.. automodule:: pyar.core.molecule
   :members:
   :undoc-members:

.. automodule:: pyar.io.xyz
   :members:
   :undoc-members:

.. automodule:: pyar.core.geometry
   :members:
   :undoc-members:

.. automodule:: pyar.molecule_merge
   :members:
   :undoc-members:

Command-line entry point
------------------------

.. automodule:: pyar.cli
   :members: argument_parse, main
   :undoc-members:

Feature request APIs
--------------------

These modules validate workflow inputs and define the persisted request payloads
used for restart compatibility. The workflow modules still run calculations.

.. automodule:: pyar.aggregation.request
   :members:
   :undoc-members:

.. automodule:: pyar.solvation.request
   :members:
   :undoc-members:

.. automodule:: pyar.conformer.request
   :members:
   :undoc-members:

xTB backends
-------------

.. automodule:: pyar.backends.xtb_utils
   :members:
   :undoc-members:

.. automodule:: pyar.backends.xtb
   :members:
   :undoc-members:

.. automodule:: pyar.backends.xtb_aiqm1
   :members:
   :undoc-members:

.. automodule:: pyar.backends.xtb_aimnet2
   :members:
   :undoc-members:

.. automodule:: pyar.backends.xtb_turbo
   :members:
   :undoc-members:

.. automodule:: pyar.backends.geometric
   :members:
   :undoc-members:

Backend capabilities and providers
----------------------------------

.. automodule:: pyar.backend_capabilities
   :members:
   :undoc-members:

.. automodule:: pyar.energy_gradient_providers
   :members:
   :undoc-members:

Sampling API
------------

.. automodule:: pyar.sampling.sphere
   :members:
   :undoc-members:

.. automodule:: pyar.sampling.rotation
   :members:
   :undoc-members:

.. automodule:: pyar.sampling.metrics
   :members:
   :undoc-members:

.. automodule:: pyar.sampling.trial_generator
   :members:
   :undoc-members:

Selection API
-------------

Selection services are documented in :doc:`generated_api` to avoid duplicating
the same object references in multiple API pages.

Scripts package
---------------

.. automodule:: pyar.scripts
   :members:
   :undoc-members:
