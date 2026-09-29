# Structure-comparison policy validation — 2026-09-29

## Dataset audit

126 geometries and 189 pairs passed the independent audit with zero errors.
Regenerating the data in the recorded environment reproduced both XYZ and CSV
byte for byte. The audit checks identifiers, references, coordinates, element
counts, SMILES-derived identity, stored MMFF energies, convergence metadata,
water geometry, and recorded atom mappings with proper rigid transforms.
Recomputed MMFF energies differed by at most 4.92e-9 kcal/mol; the largest
gradient component was 0.000733 kcal/mol/Angstrom.

Only 99 pairs have binary ground truth: 26 exact duplicates, 72 same-formula
connectivity isomers, and one different-composition pair. The other 90 pairs
are diagnostics. Six perturbations have contacts below 0.5 Angstrom and remain
flagged stress cases. They are excluded from scored accuracy.

The corpus corrects earlier fixture problems: torsion bins no longer stand in
for distinct conformers, coordinate noise no longer claims to preserve
connectivity, and the water dimer no longer contains the 0.89 Angstrom H–H
contact. The generator records exact atom mappings and excludes unconverged
MMFF optimizations.

## Comparator results

The original Coulomb eigenvalue energy rule, the Coulomb eigenvalue prefilter
plus Kabsch RMSD comparator, graph RMSD, iRMSD, and the active combined policy
were checked in both argument orders on each pair and on every self-pair:
2,520 isolated calls.
No process timeout, exception, invalid distance, or input mutation occurred.

| Method | Exact duplicates recognized | Scored negatives rejected | Abstentions |
|---|---:|---:|---:|
| Coulomb eigenvalue energy rule | 24/26 | 72/72 | 2 positives and 1 negative without energy |
| Coulomb eigenvalue prefilter plus Kabsch RMSD | 5/26 | 73/73 | 0 |
| Element-labeled graph RMSD alone | 25/26 | 73/73 | 1 positive; 2 additional diagnostic pairs incomplete |
| iRMSD 0.1.2 | 26/26 | 73/73 | 2 calls used an internal fallback |
| Graph-first plus gated bidirectional iRMSD | 26/26 | 73/73 | 0 |

For the RMSD methods, the counts were unchanged at strict thresholds 0.05,
0.1, 0.2, and 0.5 Angstrom. The original Coulomb eigenvalue energy rule
requires `|dE| < 1e-5` and sorted eigenvalue distance `< 1.0`; the benchmark
uses MMFF energies where available as a proxy. This is not a test with the
production xTB energy. Two known duplicates and one different-composition pair
have no energy metadata, so the rule abstains on them. The earlier claim that the original Coulomb method was
the worst performer was incorrect: that result came from comparing against
the newer Coulomb-prefilter plus Kabsch method.

The 72 isomer negatives are correlated across eight families and share
reference geometries. These results validate constructed invariances and
regression behavior; they are not population accuracy or threshold calibration.

### Coulomb eigenvalue energy rule

All 504 forward, reverse, and self calls completed. On scored fixtures it
recognized 24 exact duplicates and rejected all 72 scored negatives. It could
not score two positives and one negative because energy was unavailable. This
method remains a useful descriptor rule available at installation time; these fixtures do not establish
that it is universally accurate.

### Coulomb eigenvalue prefilter plus Kabsch RMSD

All calls completed, but 21 of 26 known exact duplicates were split at every
tested RMSD threshold. For example, `pair_0001` returned 1.3134 Angstrom for
a recorded rigidly transformed and permuted copy. This is a different method
from the original Coulomb eigenvalue energy rule.

### Element-labeled graph RMSD

492 pair/self calls completed and 12 were incomplete when counting both
argument orders and self checks. The scored forward cases recognized 25 of 26
exact duplicates and rejected all 73 known negatives. The known positive
abstention is Au13 pair `pair_0151`. `pair_0152` (Au13 motifs) and `pair_0153`
(perturbed Au13) are diagnostic incompletes. Six self-comparisons also hit the
10,000-isomorphism limit. Completed comparisons were direction-consistent;
all completed self-checks passed.

### iRMSD

Two calls for `pair_0169` emitted a native warning that atom topologies differed
and iRMSD fell back to quaternion RMSD without reordering atoms. Those values
are not trusted. Seven isomer pairs had direction-dependent distances differing
by 0.0326–0.0629 Angstrom, though their classifications did not change at the
tested thresholds. All self-checks passed.

## Active deduplication policy

`remove_similar` uses element-labeled graph RMSD first. If graphs match but
graph mapping enumeration reaches its cap, it runs iRMSD in both directions in
an isolated process. Native output, errors, timeouts, invalid distances, or
incomplete output cause an abstention. It uses the larger of the two distances
and deletes only if both are below the existing adaptive threshold. Graph
mismatches never invoke iRMSD. This follows “in doubt, keep.”

The combined policy completed all 504 calls; its 12 secondary iRMSD checks
returned clean output, with no directional disagreements. It recognized all 26 known duplicates and rejected
all 73 scored negatives at each tested threshold. It recovered the exact Au13
duplicate that graph RMSD alone could not finish, while retaining the Au13
icosahedron/cuboctahedron and distorted-cluster pairs. All 12 graph-incomplete
calls reached the iRMSD second check; the exact duplicate scored zero in both
orders. This evidence supports the combined policy on these fixtures, not a
claim of population-level accuracy.

## Verification and limitations

- Full test suite: **539 passed, 183 subtests passed, 553 warnings**.
- Dataset audit: **126 geometries, 189 pairs, zero audit errors**.
- Systematic comparator run: **2,520 isolated calls**.
- `git diff --check`: passed.

The corpus lacks expert-reviewed conformer/basin labels, optimized metal
clusters, and representative production reaction products. Validate the
combined policy on representative production structures before broader claims.
Results and per-call evidence are in `audit.json`, `pairwise_results.jsonl`,
and `summary.json`.
