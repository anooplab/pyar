Save and inspect computational settings
=======================================

A configuration captures a resolved request; a profile reuses methodology.
Inspect values and their sources before execution, particularly backend,
electronic state, selection controls and budgets.

.. code-block:: bash

   pyar energy water.xyz --backend xtb --write-config water-energy.toml --check
   pyar config validate water-energy.toml
   pyar run water-energy.toml --check
   pyar run water-energy.toml
   pyar profile --help

``--check`` validates without calculation or workflow output. ``--dry-run`` also
prints a plan for supported calculation commands. Neither is native scientific
validation. Explicit CLI overrides and profile/default values are recorded so
that an apparently omitted setting is not lost from provenance.

An executed workflow writes ``pyar-run.toml`` with request, settings and observed
status. Validate saved requests with the same PyAR series; unsupported formats
are rejected. Preserve the original environment for older records rather than
rewriting a schema number by hand. Never store credentials in shared profiles.
