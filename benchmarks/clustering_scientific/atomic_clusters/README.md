# Atomic cluster clustering benchmark

This pilot applies the same cluster-first, lowest-energy-per-cluster, then
max-min-budget-trim selector to two distinct atomic-cluster datasets:

* **LJ13:** 1,510 Cambridge Cluster Database local minima for the
  Lennard-Jones 13-atom landscape (reduced sigma and epsilon units).
* **Au13:** 30 gold-cluster entries from the Quantum Cluster Database's
  published NOMAD archive (coordinates in angstrom and total energies in eV).

Dataset provenance, source terms, unit conversions, and checksums are in
[`../data/SOURCES.md`](../data/SOURCES.md). Download sources with
`bash benchmarks/clustering_scientific/download_atomic_cluster_data.sh`.

## Run

From the repository root with PyAR dependencies installed:

```bash
python -m pyar.scripts.atomic_cluster_benchmark \
  --lj-archive benchmarks/clustering_scientific/data/cambridge_cluster_database/LJ13.tar.bz2 \
  --au13-archive benchmarks/clustering_scientific/data/qcd_au13/nomad_archive.zip \
  --output benchmarks/clustering_scientific/runs/atomic_clusters
```

All source structures and their provenance are written to each dataset's
`structures.csv` and `structures.xyz`. Each requested algorithm/feature
condition records the actual fallbacks, cluster-minimum count, selected
source IDs and energy interval. `manifest.json` records the source hashes.

## Interpretation and limits

The CCD and QCD files supply structures and energies, but they do not supply a
common, independently validated set of basin labels for this comparison. Each
source minimum/entry is retained as a distinct record; uncertain matches are
not collapsed. The reported coverage is the nearest selected **sorted
all-pairs distance spectrum RMS**, normalized by the input spectrum scale.
This proxy is invariant to rigid motion and atom ordering, but distinct
structures can share a pair-distance spectrum. It is not a structural RMSD,
basin recall, or proof that the selected geometry can recover a physical
minimum.

LJ13 is a useful controlled landscape but covers only one size and one model
potential. Au13 is a small real-metal set with its own mixed provenance and
source-computation limitations. This benchmark can characterize selector
behavior on these sets; it cannot establish a universal algorithm/feature
policy for metal clusters or atomic clusters generally. No downstream
optimization or energy re-evaluation is performed.
