Sequential Seed Growth
======================

Modern interface
----------------

.. code-block:: bash

   pyar grow metal.xyz ligand.xyz --count 4 --backend xtb
   pyar grow cluster.xyz H.xyz --count 12
   pyar grow seed.xyz monomer.xyz --count 4 --check

The modern command uses ``--backend`` for optional optimization. Omission runs
bounded geometry-only growth without invented energies. Orientations default to
8 and the survivor budget to 12. Output defaults to ``grow/``; ``--output DIR``
selects another directory. Every selected intermediate stage is retained.

Charge/multiplicity accept one value for both fragments or one per fragment;
omitted multiplicity follows electron parity. The existing stagewise spin rule
and site-index translation remain unchanged. ``--check`` validates the complete
request, backend dependencies when selected, and existing state/snapshots without
creating directories, changing restart state or generating structures.

Grow preserves a fixed seed and repeats one addend, whereas aggregate searches a
final composition with alternative pathways. Grow is chemistry-neutral.
Legacy ``pyar-cli grow`` retains its ``--software`` interface.
Legacy ``solvate`` remains on its existing restart-compatible path and emits
a deprecation warning; use :doc:`microsolvate` for solute-centred first-shell
construction.

``grow`` starts from one specific seed structure and repeatedly adds one chosen
species. Candidate structures are generated and selected after every addition;
all intermediate selected stages are retained. This is a chemistry-neutral
workflow for ligands, adsorbates, clusters, atoms, or solvent molecules.

In contrast, :doc:`aggregate` searches a requested final composition. No
component is permanently privileged as the seed, and PyAR may explore
alternative build pathways. Its primary result is the final composition.

Examples
--------

.. code-block:: bash

   # Search a final 2A + 3B composition
   pyar-cli aggregate A.xyz B.xyz --aggregate-size 2 3 -N 16 --software xtb

   # Add four ligands to a specified complex
   pyar-cli grow metal.xyz ligand.xyz --count 4 -N 16 \
       --maximum-number-of-seeds 5 --software xtb --xtb-model gfn2

   # Bounded geometry generation, without electronic-structure calculations
   pyar-cli grow cluster.xyz atom.xyz --count 6 -N 32 \
       --maximum-number-of-seeds 5 --connectivity-policy prefer

Inputs must be XYZ files. ``--count`` is required; orientations default to 8
and the survivor budget defaults to 12. The shared placement engine retains
its established single-orientation shortcut for monatomic seeds. Incoming
atoms have no rotational degrees of freedom. Sampling uses Fibonacci approach
directions and Halton quaternion rotations with the existing contact placement.

``--charge`` and ``--multiplicity`` accept one value for both inputs or one per
input. Each stage uses the combined fragment charge and PyAR's existing
fragment spin-combination rule, validated before work starts. Charge defaults to zero; omitted multiplicity follows electron parity.
``--site I J`` uses 0-based indices local to the original seed and addend;
the addend index is translated into the growing geometry at each step.

Optimization and selection
--------------------------

With ``--software``, the shared addition engine uses existing optimization,
cycle-limited candidate handling, staged refinement where supported, clustering,
structural comparison, and basin-memory policy. Method/model controls are
backend-specific; xTB receives no DFT method or basis defaults.

Without a backend, candidates undergo graph-first duplicate removal,
connectivity filtering, and deterministic max-min diversity selection using
existing pool descriptors and distance services. No energies are invented.
The selection algorithm override applies to energy-based clustering; the
geometry-only path always uses max-min diversity after duplicate pruning.

``--connectivity-policy auto|off|prefer|strict`` retains existing meanings.
Ordinary molecular cluster growth does not require covalent attachment;
coordinate-derived connectivity is not formal bond-order perception.

Output and restart
------------------

The default run directory is ``grow/``; override it with ``--output DIR``.
``request.json`` records the resolved request, and atomic ``state.json`` records
completed additions. ``step_000/selected/`` retains the initial seed. Each
``step_NNN/`` retains trial/jobs and a selected pool with selection diagnostics.
``final/`` holds the final pool and structured ``summary.json``.

A matching invocation resumes from the last completed addition, reusing
existing job/trial artifacts for an interrupted addition. Completed runs return
the saved result. Geometry, count, backend, sampling, site, survivor budget, or
selection changes reject reuse; use a fresh output directory. Intermediate
structures remain available. Fragment charge or spin is not inferred from
coordinate connectivity.

Python API
----------

.. code-block:: python

   from pyar.growth.request import GrowRequest
   from pyar.workflows.grow import grow

   result = grow(GrowRequest(seed, monomer, count=4,
                             number_of_orientations=16,
                             maximum_number_of_seeds=5))

``GrowResult`` reports status, state and selected paths, completed additions,
per-stage survivor counts, backend, connectivity/selection policy, and sampling.
Existing ``solvate`` semantics are unchanged.
