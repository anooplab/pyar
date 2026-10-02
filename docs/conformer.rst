Conformer Search
================

Use conformer search when you want PyAR to build and rank a set of plausible
starting geometries for a single molecule before optional backend refinement.
The workflow is RDKit-based and is aimed at flexible molecules where a single
local minimum is not enough to describe the search space.

The default settings are tuned for basin coverage: multiple RDKit seeds,
larger RDKit conformer pools, random-coordinate embedding, stratified random
torsion kicks,
tighter pruning, a larger backend-refinement window than the final output set,
and a balanced pool that protects low-energy, geometrically diverse,
contact-rich folded, open, and outlier conformer families.

Basic commands
--------------

.. code-block:: bash

   pyar-conformer CCO
   pyar-cli conformer CCO
   pyar-conformer input.sdf --software xtb

What it does
------------

* embeds multiple RDKit conformers with ETKDGv3
* perturbs rotatable-bond torsions around embedded conformers
* optionally starts from multiple seeds for broader coverage
* can keep a larger backend-refinement pool than the final output set
* can use random coordinates or a tighter prune threshold for more breadth
* optionally refines the selected pool with an existing PyAR backend

Useful options
--------------

* ``--num-conformers`` controls how many conformers are embedded per seed
* ``--num-seeds`` repeats embedding with deterministic seed offsets
* ``--backend-top-n`` widens the set sent to backend refinement
* ``--diversity-fraction`` protects part of that backend pool for graph-mapped
  RMSD diversity instead of pure RDKit energy ranking
* ``--compactness-fraction`` protects a contact-rich folded quota and a matching
  open-conformer quota, so folded structures are sampled without crowding out
  extended conformer families
* ``--use-random-coords`` starts RDKit from random coordinates instead of the
  default distance-geometry eigenvector start
* ``--torsion-kicks`` or ``--no-torsion-kicks`` controls the local
  torsion-perturbation stage
* ``--torsion-mode`` is ``random``. The sampler cycles through rotatable bonds
  and coarse angle bins while mixing one-bond, two/three-bond, and broader
  torsion kicks, then locally minimizes each trial with the selected RDKit
  force field.
* ``--torsion-rounds`` controls how many successive torsion-kick rounds are run
* ``--torsion-kicks-per-conformer`` controls how many torsion trials are made
  around each embedded conformer
* ``--torsion-dedup-rms`` removes near-duplicate trial conformers before
  backend selection
* ``--dedup-atom-mode heavy`` (default) scores heavy atoms; ``all`` includes
  hydrogens. This choice applies to generation-pool deduplication, diversity
  selection, and final deduplication. Hydrogen-only systems use all atoms.
* ``--rms-threshold`` or ``--prune-rms-threshold`` sets RDKit's greedy
  pruning threshold during embedding

Comparison policy
-----------------

Conformer comparison uses the shared graph-first deduplication policy.
Element-labelled coordinate graphs constrain atom permutations, and RMSD is
calculated after translation and proper rotation. Full composition and
connectivity are checked in both atom modes, including hydrogen attachments.
Terminal hydrogen permutations are contracted in heavy-atom mode to avoid
enumerating mappings that cannot change the scored distance. Known charge or
multiplicity differences prevent removal.

If graph mapping exceeds its limit, the shared bidirectional iRMSD check can
confirm a duplicate only after verifying its returned atom correspondence.
In heavy-atom mode, that correspondence is scored on heavy atoms after the
full mapping has passed validation. Backend errors, diagnostics, timeouts,
and missing or invalid geometries retain the candidate: **in doubt, keep**.
An uncertain diversity distance is not reported as a measured RMSD.

The existing thresholds remain in Angstrom: generation collapse uses
``max(torsion_dedup_rms, rms_threshold, 0.5)``; final selection uses
``torsion_dedup_rms`` without backend refinement and
``max(torsion_dedup_rms, 0.75)`` after refinement. A confirmed duplicate must
fall strictly below the threshold. Generation collapse keeps the lower RDKit
energy representative; final selection preserves the supplied energy ranking.
The output may contain fewer than ``top_n`` unique conformers.

RDKit's embedding pruning remains a separate generation efficiency setting;
it does not establish final identity. Set ``--rms-threshold 0`` to disable
embedding pruning when hydrogen geometry coverage matters in ``all`` mode.
Coordinate graphs encode adjacency, not bond order or chemical identity, so
these comparisons do not replace a chemical identity analysis.

Useful outputs
--------------

* ``conformers/state.json`` for restart/provenance-style metadata
* ``conformers/summary.csv`` for the ranked conformer table, including the
  basin-retention reason, nonbonded contact count, and heavy-atom radius of
  gyration when available
* ``conformers/state.json`` also records the stronger generation-collapse
  threshold used before backend refinement as ``generation_dedup_rms``
* the request records ``dedup_atom_mode`` and ``final_dedup_rms``; the
  ``comparison_policy`` section records the method, atom mode, uncertainty
  policy, and comparison/fallback counts for each stage (state schema 3)
* ``conformers/rdkit/`` for the embedded conformers before backend refinement
* ``conformers/selected/`` for the final selected conformers

Developer API note
------------------

``pyar.conformer.ConformerRequest`` owns validation and the persisted
``conformers/state.json`` request payload. ``pyar.workflows.conformer`` remains
the orchestration layer for RDKit loading, conformer generation, torsion-kick
sampling, selection, optional backend refinement, and ``ConformerResult``
creation.

Benchmarking conformer search
-----------------------------

PyAR ships a JSON benchmark runner for diagnosing why a known or reference
low-energy conformer was missed. The benchmark is intentionally diagnostic:
it does not add new search algorithms. It runs the existing conformer workflow,
compares the reference geometry against generated and selected conformers, and
classifies the result so the next algorithmic change can be chosen from data.

.. code-block:: bash

   pyar-conformer-benchmark benchmarks/conformer/small.json
   pyar-conformer-benchmark benchmark.json --num-conformers 200 --num-seeds 5 --top-n 20

The diagnosis classes are:

* ``found``: the reference geometry is present in the final selected set
* ``generated_lost_selection``: the reference was generated but not selected
* ``generated_lost_backend_or_dedup``: the reference was generated but lost
  after backend refinement or deduplication
* ``never_generated``: no generated conformer is reference-like
* ``wrong_ranking``: a reference-like conformer was generated but ranked
  outside the chosen energy window
* ``input_chemistry_issue``: the input or reference geometry prevents a
  meaningful comparison
* ``uncertain``: the available evidence is insufficient for a conservative
  diagnosis

For RDKit-only runs, ``--energy-window`` is interpreted in the native RDKit
force-field energy units used for ranking generated conformers. When backend
refinement is requested, backend energies are reported separately in Hartree
fields so RDKit ranking and backend ranking are not mixed. The per-case
``diagnosis.json`` also records the reference-like conformer's RDKit rank,
backend rank when available, refinement reason, and whether it appears to have
been removed before the final selected outputs.

The command writes ``benchmark_summary.csv``, ``benchmark_summary.json``, and
per-case ``diagnosis.json`` files under the requested output directory.

See also
--------

* :doc:`quickstart`
* :doc:`usage`
* :doc:`workflows`
