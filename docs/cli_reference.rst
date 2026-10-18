Command-line reference
======================

PyAR installs the console scripts below. ``pyar`` is the task-oriented
interface for new work. ``pyar-cli`` preserves the legacy interface. Other
scripts remain available for established workflows and specialist tools.
Use each command's ``--help`` output for its current arguments.

Core and workflow commands
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - Purpose
   * - ``pyar``
     - Modern task-oriented interface and diagnostics.
   * - ``pyar-cli``
     - Legacy aggregate, reaction, solvation, scan, and utility interface.
   * - ``pyar-react``
     - Legacy reaction-search command.
   * - ``pyar-neb``
     - Legacy staged reaction-path command.
   * - ``pyar-conformer``
     - Conformer generation and optional backend refinement.
   * - ``pyar-optimiser``
     - Standalone geometry optimization.
   * - ``pyar-trial-generation``
     - Generate trial encounter geometries.
   * - ``pyar-explore``
     - Legacy structure-exploration interface.
   * - ``pyar-aimnet2-ase-opt``
     - AIMNet2/ASE optimization helper; model assets must be installed separately.
   * - ``pyar-mlopt``
     - MLatom optimization helper; install the ``ml`` extra and backend requirements.

Analysis and reports
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - Purpose
   * - ``pyar-energy-table``
     - Print absolute and relative energies from XYZ files.
   * - ``pyar-reaction-trace``
     - Summarize a reaction trace and optionally plot it.
   * - ``pyar-basin-memory``
     - Inspect or validate saved basin-memory records.
   * - ``pyar-clustering``
     - Cluster, compare, or select XYZ structures.
   * - ``pyar-similarity``
     - Find similar structures in XYZ files.
   * - ``pyar-descriptor``
     - Calculate legacy structure descriptors; see the ``analysis`` extra.

Conformer and benchmark tools
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - Purpose
   * - ``pyar-conformer-benchmark``
     - Run the conformer benchmark harness against a supplied case file.
   * - ``pyar-benchmark-clustering``
     - Run clustering benchmark cases supplied by the user.
   * - ``pyar-benchmark-orientations``
     - Inspect encounter-orientation sampling on supplied structures.
   * - ``pyar-benchmark-ts``
     - Run the transition-state optimizer benchmark harness.

Benchmark command availability or a successful synthetic harness run does not
establish scientific validation. Use the approved benchmark cases and
acceptance criteria for any scientific claim.

For the complete command syntax, run ``pyar --help``, ``pyar-cli --help`` or
``COMMAND --help``. Installation options and external executable requirements
are described in :doc:`installation` and :doc:`external_programs`.
