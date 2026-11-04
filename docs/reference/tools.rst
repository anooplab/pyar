Specialist and developer tools
==============================

These tools expose distinct functionality through the same executable. They do
not all perform quantum-chemical calculations, and they have task-specific
outputs. Use a fresh analysis directory when a tool writes files.

.. list-table::
   :header-rows: 1

   * - Command
     - Purpose and interpretation
   * - ``pyar explore``
     - Formula-directed composite geometry generation; inspect atom counts and placement
   * - ``pyar trial``
     - Trial directions, orientations and geometric constructions, not optimized results
   * - ``pyar basin-memory``
     - Read-only geometry-backed basin registry inspection; retain source provenance
   * - ``pyar similarity``
     - Grigoryan–Springborg pool similarity and duplicate/unique outputs
   * - ``pyar descriptor``
     - Cluster descriptors; descriptor similarity is not independent chemical identity
   * - ``pyar mlopt``
     - Optional MLatom AIQM1 driver; native calculations untested; requires MNDO, DFT-D4 and Gaussian
   * - ``pyar aimnet2-ase-opt``
     - Optimization with an explicit AIMNet2 model and its model/runtime requirements
   * - ``pyar benchmark conformer``
     - Developer conformer benchmark inputs, budgets and reports
   * - ``pyar benchmark clustering``
     - Clustering benchmark campaigns, not new reference chemistry
   * - ``pyar benchmark orientations``
     - Geometric orientation sampling benchmark
   * - ``pyar benchmark ts``
     - Transition-state optimizer benchmark execution/collection

.. code-block:: bash

   pyar explore --help
   pyar trial --help
   pyar basin-memory --help
   pyar similarity --help
   pyar descriptor --help
   pyar mlopt --help
   pyar aimnet2-ase-opt --help
   pyar benchmark conformer --help
   pyar benchmark clustering --help
   pyar benchmark orientations --help
   pyar benchmark ts --help

The reference in :doc:`cli` records available arguments. Preserve threshold,
feature/model, selection and source settings with specialist outputs. Do not
interpret a descriptor ranking as an electronic energy ranking, or a benchmark
completion as a broad scientific validation. These help paths were checked; ML
backend calculations and complete benchmark campaigns were not newly run for
this documentation change.

See :doc:`backends/optional_unvalidated` for the current MLatom/AIQM1, Gaussian
and Turbomole installation requirements and native-validation limits.

Work with a pool of geometries
-------------------------------

Run these in an analysis directory containing copies of the XYZs you intend
to compare. Keep the quoted glob intact so the similarity tool sees the pool:

.. code-block:: bash

   pyar similarity --files "selected_*.xyz" --threshold 0.005
   pyar descriptor "selected_*.xyz"

The threshold demonstrates syntax; choose and report a threshold suitable for
your system. ``similarity`` writes ``Info_Duplicates.txt``,
``02Duplicates_coords.xyz`` and ``01Clean_Duplicates_coords.xyz``. It uses a
pair-distance similarity measure; this is not a graph-identity certificate.
Use pools with consistent atom composition/order assumptions and inspect the
reported pairs. It can use multiple processes; begin with a small pool.

``descriptor`` writes per-input ``.mb`` values, ``unique_cluster_properties.csv``,
``unique_files.xyz`` and ``duplicate_files.xyz``. These are geometric descriptor
outputs, not newly calculated electronic energies. Its optional Python
requirements are checked on execution. Preserve the original XYZ provenance.

Generate geometrical candidates
---------------------------------

.. code-block:: bash

   pyar trial -N 4 -i seed.xyz ligand.xyz --make-configurations
   pyar explore seed.xyz ligand.xyz --formula C2H2 --pop 2

Provide chemically appropriate seed/addend files. ``trial`` writes trial
placement XYZs such as ``mol_000.xyz``. ``explore`` generates formula-directed
composites such as ``geometry_000.xyz``. These commands construct geometries;
view atom counts, composition and contacts before passing them to an optimizer.
Their small budgets illustrate use and do not establish representative sampling.

Inspect a basin registry
--------------------------

.. code-block:: bash

   pyar basin-memory selected/basin_registry.json

Use the actual registry path emitted by your selection workflow. The command
prints JSON counts and schema information without changing the file. Unsupported
old registries are rejected; no migration is attempted. A stored representative
is an input to selection memory, not proof of an independently validated minimum.

Benchmark responsibly
-----------------------

Read ``pyar benchmark TASK --help`` before running a campaign. Prepare the task's
manifest/dataset, methods, budgets and reference provenance first. Benchmark
execution can be expensive and require external programs. Collection can report
incomplete or failed cases; do not remove them to improve aggregate success.
Keep raw results distinct from a summary table and compare only matched cases.
Model-specific optimization tools require their explicitly selected model/assets
and may have different capabilities from ``pyar optimize``. Their help and the
backend guide specify the input/output argument order; test a small system
before a batch.
