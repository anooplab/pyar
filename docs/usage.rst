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

``pyar-clustering`` supports ``mbtr``, ``soap``, and ``distance-histogram``
features. MBTR includes pair-distance and angular terms and is the default;
SOAP gives an averaged local-environment representation; the pair-distance
histogram is a NumPy-only fallback that preserves element-pair distance
distributions. Every feature is computed with one species vocabulary for the
whole input pool. Descriptor failures are recorded and trigger the next
available feature.

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

The cluster-label algorithms are ``hybrid`` (HDBSCAN, then average-linkage
agglomerative if the dependency is missing or every item is noise), ``hdbscan``,
``agglomerative``, ``dbscan``, and ``optics``. ``--distance`` selects Euclidean,
Manhattan, or cosine distance on standardized features (Euclidean is the
default). DBSCAN estimates epsilon from k-neighbour distances in that metric;
``--eps`` can override it. ``maxmin`` is a fixed-budget selector, not a
clustering algorithm.
If the requested clusterer and average-linkage fallback both fail, PyAR assigns
each structure its own label and records that last-resort fallback. This keeps
the input pool available for the existing downstream budget-trimming rule.
Labels mode writes one row per input structure, including noise label ``-1``.
JSON reports record the actual feature and algorithm and every fallback.

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

Examples::

   pyar-clustering pool/*.xyz --mode labels --feature mbtr -a hybrid -n 12 --report-output report.json
   pyar-clustering pool/*.xyz --mode labels --feature soap -a agglomerative --labels-output labels.csv
   pyar-cli -a monomer.xyz -as 2 --features soap --selection-algorithm agglomerative --selection-distance cosine

Aggregate workflows accept the same feature, clusterer, and distance settings
through ``--features``, ``--selection-algorithm``, and ``--selection-distance``.
The chosen values are stored in aggregation state and must match when resuming.
