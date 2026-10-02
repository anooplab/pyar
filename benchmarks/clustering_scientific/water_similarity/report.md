# W6 structural similarity pilot

The source slice contains **80 distinct TTM2.1-F (H2O)6 minima** in the
database's 5 kcal/mol window. The reproducible sample contains **24** source
geometries evenly spaced over the energy-sorted source list, spanning source
indices 0–79 and energies −46.532 to −41.593 kcal/mol. The full-resolution
sample artifact and all 276 pair comparisons are saved under
[`../runs/water_similarity/`](../runs/water_similarity/).

## Similarity reference

Each pair was compared with PyAR's fragment-matched global Kabsch RMSD upper
bound, including all atoms and allowing all six water fragments and the two
hydrogens within each water to permute. All 276 comparisons completed. The
distances range from 0.305 to 1.711 Å; the median is 1.138 Å. The closest pair
is source records 27 and 31 at 0.305 Å, followed by records 48 and 52 at
0.393 Å. These are continuous distances, not manually assigned same/different
labels. Source energy proximity is not treated as a structural-similarity
label.

## Feature ranking against the reference

| PyAR feature | Spearman rank correlation | Overlap among nearest 10% of pairs | Precision at 10% |
|---|---:|---:|---:|
| Distance histogram | 0.723 | 13 / 28 | 0.464 |
| MBTR | 0.672 | 18 / 28 | 0.643 |
| SOAP | 0.563 | 9 / 28 | 0.321 |

For this sample, the distance histogram ranks pairwise similarity best
overall by Spearman correlation, while MBTR retrieves more of the reference's
nearest 10% of pairs. SOAP is weaker on both measures in this pilot. The
closest reference pair (27, 31) is also among the closest pairs for all three
features. This is an exploratory result from one size and one potential; pair
distances are dependent, and these statistics do not establish a general
feature policy.

## Scope

The all-atom fragment comparison preserves water orientation and allows
identical water molecules to permute. The comparator reports an upper bound
on the global minimum RMSD under mappings; a complete mapping search does not
make that bound an exact global minimum. The benchmark tests within-composition
distance ranking. It does not test composition classification, basin recovery,
or the full clustering-and-seed-selection workflow. Provenance, source URL,
terms caveat, and checksum are documented in
[`../data/SOURCES.md`](../data/SOURCES.md).

Reproduce with:

```bash
bash benchmarks/clustering_scientific/download_water_cluster_sample.sh
python -m pyar.scripts.water_cluster_similarity_benchmark \
  benchmarks/clustering_scientific/data/water_clusters/W6_geoms_5.0_KCal-1hgztfv.txt \
  --output benchmarks/clustering_scientific/runs/water_similarity
```
