# Water-cluster structural similarity pilot

This small benchmark is sampled from the W6 within-5-kcal/mol minima file of
the University of Washington water-cluster database. It preserves the fixed
(H2O)6 composition and picks 24 evenly spaced energy ranks from the 80 source
structures, including both ends of the source energy interval. The full
80-structure source text file and its attribution/hash are documented in
[`../data/SOURCES.md`](../data/SOURCES.md).

## Build and reproduce

```bash
bash benchmarks/clustering_scientific/download_water_cluster_sample.sh
python -m pyar.scripts.water_cluster_similarity_benchmark \
  benchmarks/clustering_scientific/data/water_clusters/W6_geoms_5.0_KCal-1hgztfv.txt \
  --output benchmarks/clustering_scientific/runs/water_similarity
```

The build creates 24 XYZ structures, source IDs and energies, a complete
pairwise structural-distance table, a square reference-distance matrix, and
Euclidean matrices for MBTR, SOAP, and the distance histogram. It compares
feature-distance rankings with the reference using Spearman correlation and
nearest-pair overlap. Pair comparisons are checkpointed, so an interrupted
run can resume.

## Similarity reference and limitations

The reference is PyAR's fragment-matched global Kabsch RMSD upper bound, with
all atoms included. Coordinate graph inference must find six disconnected
H2O fragments in every record; the comparator allows water-fragment
permutations and equivalent hydrogen mappings. It reports the mapped RMSD
upper bound because its bounded mapping procedure does not guarantee the
global minimum over every possible mapping. We retain continuous pairwise
distances and do not invent similarity classes or basin labels. The source
records are distinct TTM2.1-F minima; being in the same energy window does
not imply that two structures belong to the same basin.

This first slice tests distance ranking within one water-cluster size and one
potential. It does not validate an all-chemistry aggregate policy. Water
clusters are a useful controlled start because their identical monomers allow
explicit fragment permutation matching.
