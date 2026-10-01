Usage
=====

The main command-line entry point is ``pyar-cli``.

Examples:

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

   pyar-cli -a C H -as 1 4 -N 8

   pyar-cli --aggregate --formula C5H4 -N 8

   pyar-conformer input.sdf --software xtb --num-seeds 3 --backend-top-n 50

The repository README contains additional examples for clustering,
aggregation, solvation, and reaction searches.

Energy tables
-------------

To print a relative-energy table from raw ``.xyz`` files, use:

.. code-block:: bash

   pyar-energy-table *.xyz

The command prints absolute energies, relative energies in kcal/mol, and the
global minimum directly to the terminal.

Reaction trace analysis
-----------------------

To summarize a recorded AFIR trace and optionally generate PNG plots:

.. code-block:: bash

   pyar-reaction-trace .
   pyar-reaction-trace . --plot
   pyar-reaction-trace . --plot-only
   pyar-cli trace . --plot

The command writes ``path_summary.csv`` and ``candidate_ts/`` in the job
directory unless ``--plot-only`` is used, and places plots in
``trace_plots/`` unless ``--plot-directory`` is set.

The summary distinguishes the physical backend energy from the AFIR-biased
optimization objective:

* ``backend_energy_hartree``: backend electronic, ML, or xTB energy without AFIR
* ``afir_energy_hartree``: artificial AFIR contribution
* ``total_energy_hartree``: optimization objective used by geomeTRIC
* ``backend_relative_kcalmol``: backend-energy change relative to the first
  recorded trace frame

The candidate file ``candidate_ts/highest_backend_energy.xyz`` is usually the
first structure to inspect for future NEB, string, dimer, or TS workflows. It
is not a confirmed transition state.

Utilities
---------

Several smaller helper commands are available for inspection and benchmarking:

.. code-block:: bash

   pyar-energy-table *.xyz
   pyar-clustering *.xyz -a maxmin -n 8
   pyar-clustering *.xyz --feature soap --distance cosine -a agglomerative -n 8
   pyar-clustering *.xyz --mode labels --feature distance-histogram --distance euclidean --labels-output labels.csv --report-output clustering.json
   pyar-clustering *.xyz --mode analyze --coordinate-model covalent-radii --bond-scale 1.15 --structure-report structure.json
   pyar-clustering *.xyz --mode analyze --coordinate-model distance-cutoff --bond-cutoff 3.0
   pyar-similarity -f "*.xyz" -t 0.005
   pyar-descriptor *.xyz
   pyar-conformer input.sdf --num-conformers 200 --num-seeds 4 --use-random-coords
   pyar-conformer-benchmark benchmarks/conformer/small.json
   pyar-trial-generation -N 8 -i seed.xyz monomer.xyz --plot
   pyar-optimiser structure.xyz -c 0 -m 1 --software xtb

``pyar-energy-table`` prints relative energies, ``pyar-clustering`` selects a
diverse subset, ``pyar-conformer`` generates and optionally refines RDKit
conformers, ``pyar-similarity`` reports near-duplicate structures,
``pyar-descriptor`` writes compact cluster descriptors,
``pyar-conformer-benchmark`` diagnoses why reference conformers are missed,
``pyar-trial-generation`` builds candidate orientations, and
``pyar-optimiser`` runs the standalone geometry optimizer.

Standalone structure clustering
-------------------------------

``pyar-clustering`` supports ``auto``, ``mbtr``, ``soap``, and
``distance-histogram`` features. ``auto`` uses MBTR for molecular conformers,
isomers, and uncertain XYZ-only pools, and SOAP for recognized atomic clusters
and molecular aggregates. MBTR includes pair-distance and angular terms; SOAP
gives an averaged local-environment representation; the pair-distance
histogram is a NumPy-only fallback that preserves element-pair distance
distributions. Every feature is computed with one species vocabulary for the
whole input pool. Descriptor failures are recorded and trigger the next
available feature. The inferred system class and the selected descriptor are
included in JSON reports. Use ``--system-type`` (or aggregate
``--selection-system-type``) when XYZ geometry alone does not identify the
system reliably.

Use ``pyar-clustering --mode analyze`` to inspect XYZ-only geometric
connectivity. The default ``covalent-radii`` model connects atom pairs within
``--bond-scale`` times the sum of their covalent radii; ``--coordinate-model
none`` reports atoms without assigning adjacency. For systems such as atomic
clusters, ``distance-cutoff`` accepts an explicit ``--bond-cutoff`` in
Angstrom. This operation assigns no bond orders, charges, spins, or chemical
identity. Its connected components and edge changes are geometric evidence and
must not be treated as definitive molecular or reaction labels. In aggregate
workflows, PyAR records input summaries in
``aggregates/structural_analysis/input.json`` and writes coordinate-only
input-to-output comparisons before similarity selection.

The cluster-label choices are ``auto``, ``hybrid``, ``hdbscan``,
``agglomerative``, ``dbscan``, and ``optics``. ``auto`` and ``hybrid`` try
HDBSCAN first, then average-linkage agglomerative if HDBSCAN is unavailable,
fails, or labels every item as noise. ``--distance`` selects Euclidean,
Manhattan, or cosine distance on standardized features, or ``graph-rmsd``,
``fragment-rmsd``, or ``soap-rematch`` on structures. Euclidean remains the
default. DBSCAN estimates epsilon from k-neighbour distances in the metric used;
``--eps`` can override it. ``maxmin`` requests automatic clustering followed
by max-min trimming of cluster minima when the set exceeds the seed budget.
If the descriptor has no variation, or if the requested clusterer and
average-linkage fallback both fail, PyAR assigns
each structure its own label and records that last-resort fallback. This keeps
the input pool available for the existing downstream budget-trimming rule.
Labels mode writes one row per input structure, including noise label ``-1``.
JSON reports record the actual feature and algorithm and every fallback.
Cluster mode also accepts ``--report-output``; its report includes the actual
post-filter candidate paths, labels when clustering ran, selected names, and
the selection stage. It records ``algorithm_used: not-run`` when the filtered
pool already fits the budget. All density algorithms use the same explicit
cosine distance convention for zero vectors; HDBSCAN uses a precomputed matrix
for this metric.

Structural distance matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``graph-rmsd`` uses the common graph/symmetry RMSD comparison and defaults to
heavy atoms. ``--distance-atom-mode all`` includes hydrogens; use this when
hydrogen orientation matters in molecular aggregates. Element-labelled
coordinate graphs constrain the correspondence. Incomplete mapping may use
the shared verified bidirectional iRMSD check.

``fragment-rmsd`` defaults to all atoms. It identifies molecular fragments
independently from each structure's current coordinate graph, matches
isomorphic fragments, and scores their full placement under one global proper
rotation. It retains fragment orientation and interfragment separation. It
does not independently align each monomer and discard the packing information.
Matching and rotation hypotheses are bounded; the symmetric result is a
conservative mapped RMSD upper bound, not a guaranteed global minimum.
Single-atom fragments or changed constituent graphs cannot use this backend.

``soap-rematch`` computes normalized local SOAP environments with one species
vocabulary and matches them using entropy-regularized transport. This differs
from Euclidean distance on the existing averaged SOAP vector. Its normalized
kernel dissimilarity is dimensionless. ``--soap-cutoff`` (default 5 Angstrom)
and ``--rematch-alpha`` (default 1) control the local environments and entropy
regularization. Transport has an iteration limit and uses log-domain updates.
SOAP needs DScribe and ASE. SOAP cannot distinguish reflections, and packing
changes entirely outside its cutoff can be invisible.

Clustering and max-min trimming of cluster minima use the same complete
distance matrix. Structural matrices are supplied to HDBSCAN, DBSCAN, OPTICS,
and average-linkage clustering as precomputed dissimilarities. A failing
fragment matrix tries graph RMSD, then SOAP/REMatch; a failing graph matrix
tries SOAP/REMatch. The final fallback is Euclidean distance on the chosen
standardized descriptor. The whole pool changes backend together: missing
pairs are never replaced by zero, infinity, or values in another unit system.
After a backend fallback, an explicit DBSCAN ``--eps`` is discarded and
re-estimated in the successful backend's scale. JSON reports record this,
the requested and actual distances, units, parameters, and fallback reasons.

If a descriptor distance has no resolved variation, singleton labels preserve
candidates. Magnitude-aware feature scaling avoids amplifying roundoff in
near-constant descriptor columns. Exact zero graph/fragment distances remain
valid for the selected atom scope; heavy-atom mode deliberately ignores
hydrogen geometry. These clustering distances do not define deduplication or
chemical identity.

Examples:

.. code-block:: bash

   pyar-clustering pool/*.xyz --mode labels --distance fragment-rmsd -a agglomerative --report-output packing.json
   pyar-clustering pool/*.xyz --distance soap-rematch --soap-cutoff 7 --rematch-alpha 0.5 -n 8
   pyar-clustering conformers/*.xyz --distance graph-rmsd --distance-atom-mode heavy -n 8
   pyar-cli -a monomer.xyz -as 2 --selection-distance fragment-rmsd

The aggregate workflow accepts the distance names with their documented
defaults. Custom structural parameters are available in the standalone CLI
and in ``cluster_molecules(..., distance_options={...})`` or
``choose_geometries(..., distance_options={...})``. ``ClusteringResult`` exposes
the structural ``distance_matrix``; ``pyar.selection.compute_distance_matrix``
provides the matrix directly. Descriptor features are computed only if needed
for a structural-distance fallback.

Basin memory
~~~~~~~~~~~~

Selection remembers prior basin representatives in
``selected/stoichiometry_<formula>/basin_registry.json`` (or the flat
``selected/basin_registry.json`` used by pathway selection). Schema 2 stores
the representative XYZ geometry and an explicitly versioned element-pair
distance histogram. On later runs, memory novelty is recomputed from those
geometries using the current requested feature and distance policy; old
descriptor numbers are never compared across descriptor versions. If a full
comparison matrix cannot be built, basin memory leaves the candidate pool
unchanged.

Schema-1 registries contain normalized Coulomb-spectrum fingerprints, but the
historical fingerprint routine could return a radial-coordinate fallback
instead. Migration therefore preserves each vector as an opaque
``pyar-fingerprint-v1`` legacy descriptor and does not use it to prune
candidates. It cannot reconstruct geometry that was not archived. Migrate
explicitly, first inspecting the dry-run report:

.. code-block:: bash

   pyar-basin-memory --dry-run selected/stoichiometry_C2H6/basin_registry.json
   pyar-basin-memory selected/stoichiometry_C2H6/basin_registry.json

Migration is atomic and retains legacy vectors. Future-version or malformed
registries are left untouched; selection disables memory pruning and writing
for that registry rather than risking data loss.
Workflow selection writes ``selection_diagnostics.json`` alongside selected
geometries with the memory schema, archived geometry count, ignored opaque
legacy count, requested and actual feature/distance backends, and fallbacks.

The constructed validation corpus is in ``benchmarks/clustering_distances``:
84 geometries, seven packing/torsion-family pools, and 42 independently verified
rigid-transform/permutation witnesses. Run it with:

.. code-block:: bash

   python -m pyar.scripts.benchmark_distances benchmarks/clustering_distances/dataset.json --audit-only
   python -m pyar.scripts.benchmark_distances benchmarks/clustering_distances/dataset.json --output distances.json

These are geometric regression fixtures, not potential-energy basin labels.
They do not establish a universal distance choice or production thresholds.
See the `DScribe local-kernel tutorial <https://singroup.github.io/dscribe/latest/tutorials/similarity_analysis/kernels.html>`_
and `De et al. (2016) <https://doi.org/10.1039/C6CP00415F>`_ for SOAP/REMatch.

Feature choice depends on the structures and the question. For conformers of
one molecule, a topology-aware geometry metric is preferable when distinguishing
basins; current MBTR/SOAP clustering is a structural approximation. For atomic
clusters, SOAP is a useful candidate because it compares local environments.
For molecular aggregates, averaged SOAP and pair distributions can blur
fragment arrangements, so inspect results and prefer a fragment-aware geometry
representation when available. For constitutional isomers or different
molecules with the same formula, first group by chemical identity/connectivity;
formula equality alone is not a suitable clustering feature. The standalone
module reports its representation and fallback path so these choices can be
benchmarked without changing workflow seed policy.

These recommendations follow the descriptor scope: MBTR encodes distributions
of k-body terms, while SOAP describes local atomic environments. REMatch-SOAP
was developed to compare whole structures using pairwise environment
similarities, which is more expressive than the averaged SOAP vector currently
used here. See `DScribe MBTR documentation
<https://singroup.github.io/dscribe/latest/tutorials/descriptors/mbtr.html>`_,
`DScribe SOAP documentation
<https://singroup.github.io/dscribe/latest/tutorials/descriptors/soap.html>`_,
and `De et al., Comparing molecules and solids across structural and
alchemical space <https://doi.org/10.1063/1.4940029>`_. HDBSCAN is useful when
clusters have different densities and outliers should be explicit; it does
not choose a scientifically optimal number of seed geometries. See `McInnes,
Healy, and Astels (2017) <https://doi.org/10.21105/joss.00205>`_.

The current feature-to-system mapping is a provisional engineering policy, not
a benchmark-established ranking. The existing clustering benchmark lacks
expert-labelled conformational basins across the relevant system classes.
Use the class override and JSON reports to audit a run; a dedicated labelled
benchmark is required before treating these defaults as scientifically
validated. Isomer pools with inferred graph-distinct structures are clustered
within separate element-labelled topology groups. This graph inference uses
coordinates only and does not infer bond orders, charges, or intended
fragments.
Automatic molecular-aggregate classification also partitions differing
constituent graphs, so a changed mixture is not assigned a shared cluster
solely through its descriptor.

Seed selection remains cluster-first: retain the lowest-energy geometry from
each cluster, then use max-min only if the number of cluster minima exceeds the
seed budget. If there are fewer cluster minima than requested, the selector
returns fewer seeds and does not top up from non-minimum members.

Examples::

   pyar-clustering pool/*.xyz --mode labels --feature mbtr -a hybrid -n 12 --report-output report.json
   pyar-clustering pool/*.xyz --mode labels --feature soap -a agglomerative --labels-output labels.csv
   pyar-clustering pool/*.xyz -n 8 --report-output selection.json
   pyar-cli -a monomer.xyz -as 2 --features soap --selection-algorithm agglomerative --selection-distance cosine

Aggregate workflows accept the same feature, clusterer, and distance settings
through ``--features``, ``--selection-algorithm``, ``--selection-distance``,
and ``--selection-system-type``. The chosen values are stored in aggregation
state and must match when resuming. Selection policy version 2 records the
updated distance scaling and aggregate topology partitioning. States with an
older policy version require the PyAR version that created them or a fresh
calculation directory.
The selection policy version is recorded too. Calculations with a different
or unversioned policy must be resumed with the PyAR version that created their
state, or started in a new directory; PyAR will not silently mix selection
semantics across an interrupted calculation.
