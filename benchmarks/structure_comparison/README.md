# Structure-comparison fixtures v0.2

The dataset contains **126 geometries and 189 pairs**. It covers small organic
molecules, idealized Au/Ar atomic clusters, rigid water dimers and trimers,
and coordinate perturbations. Coordinates are in Angstrom. This is a
construction-based regression corpus; it does not establish potential-energy
basin identity or calibrate production equivalence thresholds.

## Files and reproduction

- `structures.xyz`: multi-frame XYZ, identified by `structure_id` in each comment.
- `pairs.csv`: references, relation labels, confidence, and provenance.
- `generate_dataset.py`: fixed-seed geometry generation; requires RDKit.
- `systematic_test.py`: independent data audit and isolated comparator runs.
- `results/audit.json`: label evidence, geometry diagnostics, package versions,
  and hashes of data and comparison code.
- `results/pairwise_results.jsonl`: individual measurements, runtime, input
  mutation checks, native backend diagnostics, exceptions, and timeouts.
- `results/summary.json`: threshold sweeps, abstentions, directional consistency,
  and self-comparison results.
- `results/REPORT.md`: interpretation of the recorded run.

From the repository root, using the environment with RDKit and iRMSD:

```sh
./.venv/bin/python benchmarks/structure_comparison/generate_dataset.py
./.venv/bin/python benchmarks/structure_comparison/systematic_test.py --audit-only
./.venv/bin/python benchmarks/structure_comparison/systematic_test.py --jobs 4
./.venv/bin/python -m pytest -q tests/test_structure_comparison.py tests/test_structure_comparison_dataset.py tests/test_deduplication_policy.py
```

Regeneration depends on the RDKit and NumPy versions. The stored coordinates
are the authoritative inputs; hashes and installed versions accompany results.

## Label evidence

| Relation | Pairs | Evidence and interpretation |
|---|---:|---|
| `same_structure` | 26 | Recorded atom correspondence and a proper rigid transform, checked independently of the tested comparators. Includes Au13 and fragment-permuted water. |
| `different_connectivity` | 72 | Matching formula and different SMILES-derived covalent connectivity. These test chemical false merges, not an exact geometric distance. |
| `different_composition` | 1 | Identical coordinates with Au versus Ar elements. |
| `unreviewed_conformer_pair` | 42 | Different indexed torsion bins after converged MMFF optimization; no equivalence/basin judgment. |
| `coordinate_perturbation` | 32 | Eight molecules at Gaussian displacement scales 0.05, 0.15, 0.35, and 0.65 A; connectivity and basin membership are unknown. |
| `ambiguous` | 12 | Additional unoptimized coordinate perturbations, without equivalence labels. |
| `different_cluster_motif` | 1 | Idealized icosahedral versus cuboctahedral Au13; pair-distance spectra independently prove noncongruence. |
| `different_cluster_arrangement` | 2 | Constructed water dimer/trimer arrangements with different pair-distance spectra. |
| `distorted_cluster` | 1 | Perturbed Au13, without basin or equivalence assertion. |

Only the first three categories (99 pairs) receive binary classification
scores. The remaining 90 pairs are diagnostics. The 72 molecular negatives
come from eight isomer families and share reference geometries; they are
correlated fixtures, not 72 independent chemical systems. No population-level
accuracy or confidence intervals are claimed.

## Corrections from v0.1

The previous `different_conformer` labels incorrectly treated indexed torsion
bins as proof of distinct conformers. Methyl permutations and nearby bin
boundaries can create different signatures for equivalent geometries. Those
42 pairs are now explicitly unreviewed. Comparator outputs were not used to
relabel them as either duplicates or distinct conformers.

The previous `same_connectivity_distorted` label incorrectly asserted that
noise preserves physical connectivity. Those pairs now describe only the
coordinate perturbation. Six have pair separations below 0.5 A and are flagged
as pathological stress cases in `audit.json`, not physical molecular references.
Consequently, graph rejection of a noisy geometry is not automatically an error.

The original water dimer placed opposing hydrogens about 0.89 A apart. The
acceptor monomer now faces away from the donor. The audit checks monomer bond
lengths and rejects intermolecular contacts below 1.2 A. These are constructed
arrangements, not optimized hydrogen-bond minima.

The generator now rejects unavailable MMFF parameters and excludes unconverged
optimizations. MMFF energies are recomputed from the serialized geometry during
the audit. Convergence alone does not prove a Hessian-confirmed local minimum.
A conformer cross-product bug was also corrected: each side now uses its own
conformer count.

## Benchmark interpretation

The original Coulomb eigenvalue energy rule, the Coulomb eigenvalue prefilter
plus Kabsch implementation, graph RMSD, iRMSD, and the active graph-first
combined policy are evaluated in both directions for every pair and on every
self-pair: 2,520 isolated calls.
Child-process isolation captures native Fortran output and enforces a
30-second call timeout.

Parameters: graph bond scale 1.15 and maximum 10,000 isomorphisms; iRMSD inversion
off (`2`); strict RMSD thresholds 0.05, 0.1, 0.2, and 0.5 A. Threshold sweeps reuse
measured distances. Timing excludes process startup.
`coulomb_eigenvalue_prefilter_rmsd` uses sorted Coulomb eigenvalues as a
prefilter before permutation-aware Kabsch RMSD. The
`coulomb_eigenvalue_energy_rule` reconstructs the original selection rule,
using MMFF energy as a proxy where available (`|dE| < 1e-5 kcal/mol`) and
sorted Coulomb-matrix eigenvalue distance `< 1.0`; this proxy applies only to
constructed organic fixtures and is not an xTB-energy benchmark.

The active deduplication policy uses graph RMSD first. If connectivity matches
but graph mapping enumeration is incomplete, it runs iRMSD in both directions
in an isolated process. It removes a candidate only when both calls complete
without diagnostics and both distances are below the adaptive threshold; it
uses the larger distance. Graph mismatches, diagnostics, errors, or remaining
incomplete results keep both candidates. This preserves “in doubt, keep.”

Native fallback warnings, errors, incomplete graph searches, and input mutation
are recorded explicitly and excluded from successful classification counts.
Abstentions on known positives and negatives are reported separately. Direction
agreement uses a 1e-6 A tolerance; self-distance should be below 1e-6 A.

The test set still lacks expert basin labels, optimized metal-cluster
benchmarks, and representative reaction-product sets. XYZ-derived graphs
should not be assumed reliable for every perturbation or atomic cluster. The
active policy treats uncertainty conservatively by retaining candidates when
graph comparison is incomplete. A native iRMSD fallback is not a validated
iRMSD measurement.

## External sources considered

No third-party coordinates are redistributed. Organic coordinates were generated
with RDKit ETKDGv3/MMFF; cluster coordinates are analytical constructions.

[3DCS](https://huggingface.co/datasets/EscheWang/3dcs) provides conformer, chirality,
and trajectory data under CC BY-SA 4.0. [DockRMSD](https://pmc.ncbi.nlm.nih.gov/articles/PMC6556049/)
provides a complementary ligand-pose mapping benchmark. Both are candidates for
separate external-validation adapters.
