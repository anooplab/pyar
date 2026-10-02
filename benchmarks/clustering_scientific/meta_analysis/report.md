# Cross-domain clustering meta-analysis

This analysis compares the molecular conformer, atomic-cluster, and water-cluster pilots. It keeps each domain's reference label and outcome separate; a pooled AUC across these unlike targets would be misleading.

## Reference evidence by domain

| Domain | Available reference | What can be measured | Limitation |
|---|---|---|---|
| Molecular conformers | Independent GFN2-xTB optimizations; graph identity plus heavy-atom Kabsch RMSD below 0.35 Å | Basin-pair ROC/AUC and confusion matrices for each clusterer-feature setting | Three macrocycles; pair counts are dependent and class balance differs |
| LJ13 and Au13 atomic clusters | Distinct source minima/entries, energies, and pair-spectrum coverage | Source-energy retention and normalized geometric coverage | No independently validated same/different similarity labels; no sensitivity/specificity/AUC |
| W6 water clusters | All-atom fragment-matched Kabsch RMSD upper-bound distances | Feature ROC/AUC and confusion matrices at declared distance cutoffs | One composition; proxy labels, not verified basin labels; reference distances are upper bounds |

## Conformer feature ROC/AUC

The positive class is a pair assigned to the same operational reference basin. Lower feature distances predict a positive pair. Values are ROC AUC by source pool and the unweighted macro-average.

| Feature | CAMVES_I | FGG55 | WG01 | Macro AUC |
|---|---:|---:|---:|---:|
| mbtr | 0.871 | 0.997 | 0.986 | 0.951 |
| soap | 0.905 | 0.998 | 0.998 | 0.967 |
| distance-histogram | 0.903 | 0.971 | 0.942 | 0.939 |

soap has the highest observed same-basin pair AUC in this three-pool comparison. AUC ranks pair distances; it does not by itself establish a safe merge threshold.

## Cluster-assignment confusion summary

For each clusterer-feature condition, a pair is predicted positive if it has the same non-noise cluster label. Noise points are treated as separate retained candidates. The table below macro-averages metrics across the three conformer pools; undefined metrics in a pool are omitted from that metric's mean.

| Algorithm | Feature | Sensitivity | Specificity | Precision | False-positive rate | False-discovery rate | Miss rate | Basin budget efficiency |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| optics | distance-histogram | 0.464 | 0.974 | 0.547 | 0.026 | 0.453 | 0.536 | 0.917 |
| optics | mbtr | 0.645 | 0.972 | 0.559 | 0.028 | 0.441 | 0.355 | 0.889 |
| optics | soap | 0.677 | 0.972 | 0.553 | 0.028 | 0.447 | 0.323 | 0.917 |
| dbscan | soap | 0.695 | 0.969 | 0.476 | 0.031 | 0.524 | 0.305 | 0.917 |
| dbscan | mbtr | 0.659 | 0.965 | 0.426 | 0.035 | 0.574 | 0.341 | 0.972 |
| auto | soap | 0.730 | 0.957 | 0.403 | 0.043 | 0.597 | 0.270 | 0.917 |
| maxmin | soap | 0.730 | 0.957 | 0.403 | 0.043 | 0.597 | 0.270 | 0.917 |
| agglomerative | soap | 0.667 | 0.910 | 0.089 | 0.090 | 0.911 | 0.333 | 1.000 |
| agglomerative | mbtr | 0.658 | 0.871 | 0.048 | 0.129 | 0.952 | 0.342 | 0.972 |
| auto | mbtr | 0.749 | 0.830 | 0.321 | 0.170 | 0.679 | 0.251 | 0.778 |
| maxmin | mbtr | 0.749 | 0.830 | 0.321 | 0.170 | 0.679 | 0.251 | 0.778 |
| dbscan | distance-histogram | 0.687 | 0.720 | 0.347 | 0.280 | 0.653 | 0.313 | 0.972 |
| auto | distance-histogram | 0.714 | 0.580 | 0.233 | 0.420 | 0.767 | 0.286 | 0.917 |
| maxmin | distance-histogram | 0.714 | 0.580 | 0.233 | 0.420 | 0.767 | 0.286 | 0.917 |
| agglomerative | distance-histogram | 0.658 | 0.554 | 0.015 | 0.446 | 0.985 | 0.342 | 1.000 |

The class imbalance matters: FGG55 and WG01 contain many more different-basin pairs than same-basin pairs, so specificity can look high while precision remains low. Use the per-pool 2×2 counts in `conformer_pair_confusion.csv` alongside these rates. The seed-budget basin efficiency divides recovered basins by the maximum distinct basins the fixed 12-seed budget could retain; it avoids calling 12/98 a 12% failure when the budget itself is 12.

## Water-cluster proxy ROC/AUC and conservative operating points

Reference-positive pairs are those with fragment-RMSD upper bound at or below the listed cutoff. For the 0.75 Å proxy cutoff, feature thresholds are selected on this same 24-structure sample to maximize sensitivity subject to the listed minimum specificity. These are apparent, in-sample operating points, not validated deployment thresholds.

| Feature | Reference cutoff (Å) | Positive pairs | AUC | Specificity target | Actual specificity | Sensitivity | TP | FP | FN | TN |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| mbtr | 0.75 | 15 | 0.956 | — | — | — | — | — | — | — |
| mbtr | 0.75 | 15 | 0.956 | 0.95 | 0.969 | 0.667 | 10 | 8 | 5 | 253 |
| mbtr | 0.75 | 15 | 0.956 | 0.99 | 0.996 | 0.333 | 5 | 1 | 10 | 260 |
| soap | 0.75 | 15 | 0.851 | — | — | — | — | — | — | — |
| soap | 0.75 | 15 | 0.851 | 0.95 | 0.950 | 0.400 | 6 | 13 | 9 | 248 |
| soap | 0.75 | 15 | 0.851 | 0.99 | 0.992 | 0.200 | 3 | 2 | 12 | 259 |
| distance-histogram | 0.75 | 15 | 0.885 | — | — | — | — | — | — | — |
| distance-histogram | 0.75 | 15 | 0.885 | 0.95 | 0.958 | 0.600 | 9 | 11 | 6 | 250 |
| distance-histogram | 0.75 | 15 | 0.885 | 0.99 | 1.000 | 0.333 | 5 | 0 | 10 | 261 |

At the 0.75 Å proxy cutoff there are 15 positive pairs and 261 proxy-negative pairs. The high-specificity points deliberately accept false negatives: with a keep-in-doubt rule, an uncertain pair is retained rather than merged. Because the reference distance is an upper bound, proxy-negative pairs are not proven dissimilar; the resulting specificity is against this operational proxy only.

## Atomic-cluster outcomes

All 15 tested conditions on each atomic dataset retained its lowest-source-energy entry. Best normalized pair-spectrum coverage differs by dataset:

| Dataset | Best setting by pair-spectrum coverage | Selected | Normalized mean nearest-spectrum RMS |
|---|---|---:|---:|
| lj13 | dbscan / mbtr | 12 | 0.05097 |
| qcd_au13 | auto / soap | 12 | 0.02050 |

Sorted pair-distance spectra are not unique structural identifiers, so these coverage results cannot be converted into confusion matrices without new independent labels.

## Provisional policy

1. **Deduplication:** retain the conservative `in doubt, keep` rule. Use feature distances only to find candidate neighbors; merge only after a complete, chemistry-aware structural comparison confirms equivalence. Incomplete graph mappings or threshold-borderline matches remain separate.
2. **Conformer feature:** SOAP is the strongest tested basin-pair ranking feature across the three macrocycle pools. Use it as the leading candidate for further held-out testing, not as a universal default yet.
3. **Clustering algorithm:** no algorithm-feature pair dominates the cluster confusion, specificity, and seed-budget metrics. The HDBSCAN-first `auto` policy is not independently validated against DBSCAN or OPTICS by this analysis; `maxmin` shares the automatic/HDBSCAN cluster labels and is a cluster-then-trim selection mode, not a separate partitioner.
4. **Atomic and aggregate systems:** keep feature/algorithm choices system-aware. Two atomic datasets disagree on their best proxy setting, and the W6 benchmark shows descriptor ranking depends on whether the objective is global rank correlation or retrieval of the nearest pairs.
5. **Operating point:** the high-specificity water settings miss most proxy-positive pairs in this sample. Do not adopt these in-sample thresholds; evaluate retrieval recall and false-discovery rate on held-out cluster sizes and a chemically diverse aggregate set before choosing an operating point.
6. **Cutoff interpretation:** the 0.35 Å conformer reference label is not the production duplicate-removal cutoff. Current generation and backend-final reduction floors are 0.50 and 0.75 Å, respectively; this benchmark does not establish those broader removals as safe. Keep the distinction explicit until direct deletion-safety labels are available.

## Reproducibility and caution

The pairwise entries are not independent observations, and this is a small exploratory set. Conformer reference basins are operational labels produced by the saved graph/Kabsch protocol. Water thresholds are RMSD-proxy categories. Atomic-cluster similarity labels remain unavailable. No single pooled accuracy number is reported across these different targets.

Rebuild this analysis with `python -m pyar.scripts.meta_analyze_clustering_benchmarks --output benchmarks/clustering_scientific/meta_analysis`.
