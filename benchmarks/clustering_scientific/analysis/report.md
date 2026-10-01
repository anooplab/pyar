# Scientific clustering benchmark results

This is a macrocycle conformer-selection pilot. Reference labels come from independently GFN2-xTB-optimized geometries grouped by complete coordinate-graph isomorphism and heavy-atom Kabsch RMSD. See the benchmark README for attribution, methods, and limits.

The run covers 209 starting geometries, 145 reference basins at 0.35 Å, 45 algorithm/feature-setting runs, and 89 unique perturbed seed optimizations.

The highest weighted basin retention was **agglomerative + soap**: 28/145 basins (19.3%). It retained 6/10 basins within 1 kcal/mol and 18/58 within 3 kcal/mol of each pool's lowest reference energy.

Across all 450 selected-seed/settings occurrences, xTB optimization success was 100.0% and verified recovery of each seed's own basin was 100.0%. These are robustness checks at the tested 0.10 Å perturbation level.

## Test systems

| System | Geometries | Reference basins |
|---|---:|---:|
| CAMVES_I | 10 | 4 |
| FGG55 | 127 | 98 |
| WG01 | 72 | 43 |

## Reference-label sensitivity and source audit

| System | Basins at 0.25 Å | Basins at 0.35 Å | Basins at 0.50 Å | Inferred topology groups | Max comment/xTB energy difference (Eh) |
|---|---:|---:|---:|---:|---:|
| CAMVES_I | 5 | 4 | 3 | 1 | 3.22e-07 |
| FGG55 | 103 | 98 | 92 | 1 | 7.55e-07 |
| WG01 | 46 | 43 | 37 | 1 | 7.96e-07 |

All source comment energies match fresh xTB final energies numerically within `8e-7 Eh`. This conflicts with the upstream README's Egret-1/kcal-per-mole description; fresh GFN2-xTB Hartree energies are used for cluster minima, and source values are retained for audit.

## Settings (ranked)

| Rank | Algorithm | Feature | Basin recall (weighted) | System-macro recall | Basin recall (≤1 kcal/mol, weighted) | Basin recall (≤3 kcal/mol, weighted) | Mean nearest indexed heavy-atom RMSD (Å) | Perturbed optimization success | Verified basin recovery | Fallback systems (alg/feat/dist) | Mean selection time (s) |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | agglomerative | soap | 0.193 | 0.467 | 0.600 | 0.310 | 0.592 | 1.000 | 1.000 | 0/0/0 | 0.077 |
| 2 | agglomerative | distance-histogram | 0.193 | 0.467 | 0.500 | 0.241 | 0.759 | 1.000 | 1.000 | 0/0/0 | 0.110 |
| 3 | dbscan | mbtr | 0.186 | 0.459 | 0.500 | 0.121 | 0.609 | 1.000 | 1.000 | 0/0/0 | 2.443 |
| 4 | agglomerative | mbtr | 0.186 | 0.459 | 0.500 | 0.190 | 0.624 | 1.000 | 1.000 | 0/0/0 | 2.417 |
| 5 | dbscan | distance-histogram | 0.186 | 0.459 | 0.500 | 0.224 | 0.774 | 1.000 | 1.000 | 0/0/0 | 0.110 |
| 6 | optics | soap | 0.186 | 0.384 | 0.400 | 0.207 | 0.627 | 1.000 | 1.000 | 0/0/0 | 0.190 |
| 7 | dbscan | soap | 0.186 | 0.384 | 0.400 | 0.155 | 0.629 | 1.000 | 1.000 | 0/0/0 | 0.080 |
| 8 | hybrid | soap | 0.186 | 0.384 | 0.400 | 0.190 | 0.635 | 1.000 | 1.000 | 0/0/0 | 0.099 |
| 9 | maxmin | soap | 0.186 | 0.384 | 0.400 | 0.190 | 0.635 | 1.000 | 1.000 | 0/0/0 | 0.104 |
| 10 | optics | distance-histogram | 0.186 | 0.384 | 0.400 | 0.207 | 0.775 | 1.000 | 1.000 | 0/0/0 | 0.129 |
| 11 | hybrid | distance-histogram | 0.186 | 0.384 | 0.400 | 0.207 | 0.782 | 1.000 | 1.000 | 0/0/0 | 0.111 |
| 12 | maxmin | distance-histogram | 0.186 | 0.384 | 0.400 | 0.207 | 0.782 | 1.000 | 1.000 | 0/0/0 | 0.116 |
| 13 | optics | mbtr | 0.179 | 0.376 | 0.400 | 0.103 | 0.621 | 1.000 | 1.000 | 0/0/0 | 3.358 |
| 14 | maxmin | mbtr | 0.152 | 0.345 | 0.400 | 0.103 | 0.683 | 1.000 | 1.000 | 0/0/0 | 2.458 |
| 15 | hybrid | mbtr | 0.152 | 0.345 | 0.400 | 0.103 | 0.683 | 1.000 | 1.000 | 0/0/0 | 2.976 |

## Scope and caveats

- Only macrocycle conformer pools are represented.
- The independent clustering benchmark excludes the workflow's pre-clustering deduplication stage.
- Reference basins depend on coordinate-inferred connectivity and the stated RMSD cutoff.
- Heavy-atom RMSD ignores hydrogen-coordinate-only rearrangements, while retaining hydrogen connectivity/count constraints.
- Downstream optimization success is measured after deterministic 0.10 Angstrom Gaussian coordinate perturbations.
- Successful optimization is distinct from recovery of the selected seed's basin; complete-comparison basin changes and incomplete graph mappings are reported separately.
- Results do not establish policy for isomer pools, atomic clusters, or non-covalent aggregates.

The ordering is descriptive for these three test pools, not a universal feature/algorithm policy. Requested methods can fall back; per-system actual methods are recorded in `summary.json`.
