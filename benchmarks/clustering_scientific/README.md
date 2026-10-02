# Scientific clustering benchmark (pilot)

This is a reproducible pilot for PyAR's cluster-first seed selection. It
benchmarks algorithm and descriptor choices against independently optimized
conformer basins, then checks geometric coverage and basin recovery after a
perturbed downstream optimization. It is not a universal policy decision.

## Data and attribution

The checked-in XYZ files are three complete conformer ensembles from
[MPCONF196GEN](https://github.com/rowansci/MPCONF196GEN-benchmark): CAMVES_I,
FGG55, and WG01. Upstream describes its data and code as CC BY 4.0. The local
copy includes the upstream license and README in `data/mpconf196gen/`; cite the
benchmark as described there. SHA-256 hashes are recorded in each run manifest.

These are CREST-generated macrocycle conformer pools, not broad coverage of
atomic clusters, constitutional isomers, or non-covalent aggregates. The three
pools contain 10, 127, and 72 geometries, respectively. We preserve every
source geometry and do not treat the upstream frame order as a basin label.

The upstream README says its conformer energies are Egret-1 single-point
energies in kcal/mol. In these three ensembles, bare numeric XYZ comments
instead match independently recomputed GFN2-xTB final energies in Hartree to
within `8e-7 Eh`, including rounding and geometry relaxation. That conflicts
with the upstream description. Seed representatives are ranked by fresh xTB
optimization energy; both values and their differences are retained for audit.

## Reference basin protocol

Each source geometry is optimized independently with GFN2-xTB, neutral charge,
and singlet multiplicity. The harness requires xTB's normal-termination line,
a final energy, and gradient norm at or below `1e-3 Eh/alpha`. Reference basin
labels are then assigned by complete element-labelled coordinate-graph
isomorphism followed by heavy-atom Kabsch RMSD below `0.35 Å`. Terminal
hydrogen connectivity and counts constrain graph matching, but hydrogen
coordinates do not contribute to RMSD; hydrogen-only rearrangements are
outside this pilot's basin definition. Incomplete graph mappings are never
merged. The label artifact also reports sensitivity at `0.25`,
`0.35`, and `0.50 Å`; labels remain a reviewed operational reference, not an
experimental truth. Bond connectivity is inferred from coordinates, so graph
inference and the RMSD threshold remain sources of uncertainty.

Energy-window recall is measured relative to the lowest freshly optimized
GFN2-xTB basin present in each selected ensemble. It does not estimate coverage
of the full potential-energy surface.

Each condition calls PyAR's `cluster_molecules`, retains the lowest-xTB-energy
member of each non-noise cluster, keeps noise points as independent candidates,
then applies the production max-min budget trim when needed. The separate
pre-clustering deduplication stage is excluded to isolate clustering and
feature quality; these results therefore do not measure the complete workflow
selector. The budget is 12 seeds. The primary outcome
is the number and fraction of independently labelled basins represented among
the selected seeds. Coverage is the nearest RMSD from each reference basin to
any selected optimized seed using source-order-corresponding
heavy-atom Kabsch RMSD (the data preserve a common atom order). This coverage
distance does not enumerate graph automorphisms; graph-constrained basin labels
and basin retention are reported separately. For downstream robustness,
each selected geometry receives deterministic Gaussian coordinate noise
(`sigma=0.10 Å`) and is optimized again with GFN2-xTB; the report measures
successful terminations separately from verified recovery of each seed's
own unperturbed basin. A complete graph comparison above the cutoff is
reported as a basin change; a graph mapping that hits its limit is uncertain.
Optimizer convergence alone does not count as basin recovery.

The tested factorial is algorithms `auto` (HDBSCAN first with agglomerative
fallback), `agglomerative`, `dbscan`,
`optics`, and `maxmin` crossed with `mbtr`, `soap`, and `distance-histogram`.
All feature and algorithm fallback diagnostics are retained. Requested method
names must be interpreted alongside the actual method and fallback reported
by PyAR.

## Reproduce

From the repository root, with the PyAR environment active and `xtb` on PATH:

```bash
bash benchmarks/clustering_scientific/run_benchmark.sh
```

This stores optimization inputs, outputs, and logs, reference labels, per-
condition JSON/CSV, manifests, and the aggregate analysis under
`benchmarks/clustering_scientific/runs/`. No run data is written to `/tmp`.
Use the Python module directly for a smaller diagnostic run:

```bash
python -m pyar.scripts.scientific_clustering_benchmark \
  benchmarks/clustering_scientific/data/mpconf196gen/FGG55_crest_conformers.xyz \
  --output benchmarks/clustering_scientific/runs/FGG55 \
  --algorithms agglomerative --features mbtr --max-seeds 12
```

`--limit` is only for harness debugging and must not be used for reported
benchmark results. The run shell script pins no software binaries; manifests
capture the xTB version string and the PyAR source revision should be recorded
by the caller when archiving/publishing results.

## Interpreting the pilot

The macro-average over three macrocycles can compare candidate settings for
this conformer-selection use case only. A broader scientific policy still
needs independently labelled pools for small-molecule conformers/isomers,
atomic clusters, and molecular aggregates, with a separate reference protocol
for each structural class. This pilot should not be presented as validating a
single feature or algorithm for all those classes.

## Atomic-cluster extension

The LJ13 and Au13 pilot uses published minima and cluster entries without
claiming verified basin labels. It records geometric pair-spectrum coverage,
source energies, and provenance. See
[`atomic_clusters/README.md`](atomic_clusters/README.md) for interpretation,
reproduction, and limitations, and [`data/SOURCES.md`](data/SOURCES.md) for
source attribution and checksums.

For a fixed-composition aggregate similarity study, see the small W6 water
cluster pilot in [`water_similarity/README.md`](water_similarity/README.md)
and its [results](water_similarity/report.md).

## Cross-domain meta-analysis

The exploratory cross-domain synthesis is in
[`meta_analysis/report.md`](meta_analysis/report.md), with machine-readable
tables and JSON alongside it. It reports conformer pairwise ROC/AUC and
cluster-assignment confusion matrices, water-cluster proxy ROC/AUC and
high-specificity operating points, and atomic-cluster coverage metrics. The
three domains are kept separate because their reference labels differ and the
atomic-cluster pilot has no validated same/different labels. Rebuild it with:

```bash
python -m pyar.scripts.meta_analyze_clustering_benchmarks \
  --output benchmarks/clustering_scientific/meta_analysis
```
