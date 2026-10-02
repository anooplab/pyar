# Review of selection, deduplication, and validation

Review date: 2 October 2026. Window: 1 October 03:23:44 to 2 October
11:23:44 Asia/Kolkata, plus the working tree at review time. Assessment:
**Needs revision before using these results to choose scientific defaults.**
The pilots remain useful exploratory evidence. This is a targeted review, not
a complete audit of PyAR or verification of all external source calculations.
Production files and original benchmark results were preserved; this directory
contains the review, runnable checks, and isolated reproduction fixtures.

## Changes covered

| Commit | Change | Review emphasis |
|---|---|---|
| `65755e2` | Structural clustering distances | Graph/fragment RMSD, SOAP/REMatch, automatic features, topology grouping, fallback behavior |
| `4db6577` | Geometry-backed basin memory | Archive insertion, novelty filtering, conformer comparison integration |
| `f48a6b6` | Scientific clustering pilot | Optimized reference labels, algorithm/feature comparison, downstream checks |
| Uncommitted | Atomic/water pilots, meta-analysis, algorithm naming | Reference quality, metrics, reuse safety, `auto`/HDBSCAN terminology |

The older deduplication components were inspected where these changes call
them. An inherited policy issue is identified below as inherited, rather than
being attributed to these three commits.

## Implemented policy versus the evidence

| Stage | Current behavior | Evidence gap |
|---|---|---|
| Aggregation/standalone deduplication | All-atom graph RMSD; validated bidirectional iRMSD after incomplete graph mapping; adaptive cutoff 0.05–0.15 Å, default 0.10 Å without samples | The new real-data clustering pilots exclude this removal stage |
| Conformer generation | Heavy atoms by default; generation cutoff at least 0.50 Å | Larger than the 0.35 Å reference-group cutoff; generation floor survives the comparison refactor |
| Backend-refined conformers | Heavy atoms by default; final cutoff at least 0.75 Å | Not calibrated by the reported clustering pilot |
| Clustering | `auto` attempts HDBSCAN, then average-linkage agglomerative; singleton preservation if no usable labels | Operational robustness is established more clearly than scientific superiority |
| Feature `auto` | MBTR for molecular/isomer/unknown pools; SOAP for atomic and molecular aggregates | A provisional mapping, not a learned or validated universal choice |
| Distances | Standardized-feature Euclidean by default; graph RMSD, fragment RMSD and SOAP/REMatch optional | Real-data pilots compare Euclidean descriptors; the structural-distance comparison is predominantly synthetic |
| Selection | Deduplication, optional archive novelty reduction, then clustering when over budget; minima and noise candidates; max-min trims excess | Full workflow, memory effects, and budget-dependent bypasses are not measured by the isolated clustering harness |

The Coulomb-spectrum fingerprint still helps select the sample used to estimate
the aggregation RMSD cutoff; that sample is scored with permutation-aware
Kabsch RMSD. Final duplicate removal uses graph/iRMSD evidence. Therefore the
cutoff estimator and final comparator are not yet the same method.

## Priority findings

### 1. Repeated equivalent geometries can crowd distinct entries out of basin memory

**P1, demonstrated regression.** Archive insertion hashes the raw coordinates
and atom ordering in
[`basin_memory.py`](../../../pyar/selection/basin_memory.py), lines 297–339.
A translation, rotation, or reordering gets a new storage ID; capacity trimming
then keeps the last entries. The previous archive suppressed identical invariant
fingerprints, although that approach had its own collision problems.

The saved probe used a capacity of three. After inserting one distinct geometry
and repeated translations of another, only the translated copies remained.
Graph RMSD independently confirmed the copies were equivalent. This can waste
the normal 200-entry capacity and change subsequent novelty selection.

**Proposed fix:** retain raw hashes for provenance, but recognize verified
equivalent representatives before insertion, update their metadata, and retain
uncertain cases. Use an invariant descriptor only to propose comparisons;
never make descriptor collisions sufficient for removal. Define archive
eviction in terms of representative diversity rather than duplicate snapshots.
Evidence: [policy_probe.json](policy_probe.json) and
[comparison_probe.json](comparison_probe.json).

### 2. Benchmark resume can reuse results for different source geometries

**P1, demonstrated reuse defect.**
[`scientific_clustering_benchmark.py`](../../../pyar/scripts/scientific_clustering_benchmark.py),
lines 106–107, accepts a completed `frame_N` optimization without checking its
input coordinates or calculation identity.
[`water_cluster_similarity_benchmark.py`](../../../pyar/scripts/water_cluster_similarity_benchmark.py),
lines 164–169, checks pair indices and completion but not the source content or
comparator configuration before reuse; it subsequently writes a fresh source hash.

In an isolated fixture, a new 4.0 Å H–H input returned a cached 0.74 Å geometry
without launching an optimizer. In another fixture, changing water coordinates
preserved the old pair distance while the manifest acquired the new source hash.
These fixtures deliberately mock expensive computation to test reuse behavior.

**Proposed fix:** version cache records and check normalized input identities,
coordinates, source checksum, calculation/comparator settings, and implementation
or executable identity before accepting results. Reject or rebuild mismatches;
do not relabel stale calculations with new provenance.

This is a reproducible risk, not evidence that the saved conformer runs are
already stale: all 209 saved optimization inputs match the current source
coordinates, and their reference source hashes match. Evidence:
[benchmark_probe_results.json](benchmark_probe_results.json).

### 3. Production conformer cutoffs and reference labels answer different questions

**P1 for a deduplication policy decision; inherited cutoffs.**
[`conformer/request.py`](../../../pyar/conformer/request.py), lines 130 and 179,
enforces generation and backend-final floors of 0.50 and 0.75 Å. These floors
predate the reviewed comparison refactor. Setting the exposed torsion cutoff
below them does not lower those effective thresholds.

Applying the current heavy-atom graph comparator in energy order to the full
saved optimized pools gives:

| Pool | Reference groups at 0.35 Å | Representatives at 0.50 Å | Representatives at 0.75 Å |
|---|---:|---:|---:|
| CAMVES_I | 4 | 3 | 3 |
| FGG55 | 98 | 92 | 68 |
| WG01 | 43 | 37 | 30 |

This is a cutoff replay, not a rerun of the full conformer-generation workflow.
It shows the policy mismatch; it does not establish that every collapsed group
is a distinct physical minimum. Even the 0.35 Å reference groups contain
optimized energy spans up to 2.31 kcal/mol in FGG55 and 2.27 kcal/mol in WG01.
Heavy-atom neighborhoods can include differing hydrogen arrangements, and the
greedy representative rule is not an independent proof of basin identity.

**Proposed fix:** expose effective generation/final duplicate cutoffs directly,
distinguish duplicate identity from deliberately coarse conformer reduction,
and evaluate removal safety against reviewed geometry/energy evidence. Keep
uncertain pairs. Do not replace 0.75 Å with 0.35 Å merely because the benchmark
currently uses the latter. Evidence: [comparison_probe.py](comparison_probe.py).

### 4. Specificity and the current “false-merge rate” do not measure deletion safety

**P1 interpretation issue.**
[`meta_analyze_clustering_benchmarks.py`](../../../pyar/scripts/meta_analyze_clustering_benchmarks.py),
line 52, names `FP/(FP+TN)` a false-merge rate. This is the false-positive rate.
The fraction of proposed matches that are wrong is `FP/(TP+FP)`, or one minus
precision. Both denominators are useful, but answer different questions.

For the water SOAP operating point at the 0.75 Å proxy cutoff and 0.99
specificity target, the independently recomputed counts are TP=3, FP=2, FN=12,
TN=259. Specificity is 99.23%, yet 2/5 proposed matches disagree with the proxy.
The reference is an RMSD upper bound, so those disagreements are not verified
chemical false deletions. They still demonstrate why high specificity alone
does not establish a safe operating point. Across the three descriptors, these
0.99-specificity operating points miss 67–80% of proxy-positive pairs; that is
not a demonstrated small miss rate.

In conformer clustering, same-cluster pairs are predictions of co-clustering,
not actual deduplication actions. A cluster can legitimately contain distinct
geometries. Feature ranking, cluster coverage, and duplicate deletion must have
separate outcomes.

**Proposed fix:** report FPR, precision/false-discovery fraction, recall, and
actual counts; evaluate actual removals separately. For candidate-neighbor
screening, prioritize retrieval recall so possible duplicates reach the
confirming comparator. For final deletion, prioritize confirmed false-deletion
risk and abstentions. The definitions agree with the
[scikit-learn metrics guide](https://scikit-learn.org/stable/modules/model_evaluation.html).
Independent calculations are in [policy_probe.py](policy_probe.py).

### 5. The benchmark comparison does not yet identify the best production policy

The source data contain 209 conformers, 1,510 LJ13 minima, 30 Au13 entries, and
24 sampled W6 geometries. The main limitation is the targets and coverage of
methods, not just the number of structures.

- The conformer and atomic harness fixes descriptor Euclidean distance. It does
  not compare the newly implemented graph RMSD, fragment RMSD, and SOAP/REMatch
  on these real-data selection tasks. Water measures descriptor ranking only.
- `maxmin` and `auto` run the same clustering and trimming policy. They are
  duplicate conditions, not an independent max-min-only control.
- Production bypasses clustering when the filtered pool fits the budget
  (`selection/clustering.py:209`). The harness always clusters. CAMVES_I has
  ten geometries and a twelve-seed budget, so its clusterer ranking is not a
  behavior production would execute for that pool.
- HDBSCAN/SOAP selects four LJ13 representatives while several alternatives
  select twelve. Coverage rankings combine representation quality with the
  different numbers selected. Report coverage versus actual cost/budget and
  underfilling explicitly; do not top up production minima contrary to policy.
- Retaining the global energy minimum is enforced by cluster-minimum selection
  and the first max-min anchor. It is a useful regression invariant, not evidence
  favoring an algorithm. Au13 energy rankings also need a verified common
  electronic-structure protocol before serving as a scientific oracle.
- Explicit algorithms still have operational fallbacks. A real DBSCAN
  all-noise probe returned agglomerative labels. For scientific comparisons,
  evaluate named methods with substitutions disabled or report the method as
  unavailable, and evaluate the automatic fallback policy separately.
- Meta-analysis allows descriptor fallback, groups by requested names, and has
  hard-coded winning-feature statements. Current saved descriptors did not
  fall back, but regeneration in another environment can make the narrative
  untrue. Derive claims from the actual method and record hashes/versions.

Add energy-only, repeated random, and genuine max-min-only experimental
controls. Include both isolated stages and the production pipeline with
deduplication and archive memory. These controls are benchmark arms; they do
not require changing the user-approved production sequence.

### 6. Aggregate-distance scaling needs a deliberate next experiment

The fragment comparator has a fixed default limit of 720 fragment assignments.
Seven identical water fragments require 7! assignments, so even an identical
seven-water geometry immediately reaches the fragment-permutation limit.
This is safe abstention and fallback behavior, but the successful W6 study is
exactly at the largest fully permuted identical-fragment count under that limit.
It does not validate performance on larger clusters. Increasing
`max_mappings` alone does not change the separate permutation limit.

Mixed-topology pools request a full structural matrix before topology grouping.
A simple ethanol/dimethyl-ether pool therefore falls back from graph RMSD to
SOAP/REMatch even though the two topologies are subsequently kept separate.
Consider distances within each topology group, with a separately defined
between-group selection policy.

The water harness evaluates one fragment-comparison direction and mirrors it;
production uses the larger of both directions. Three pairs nearest the 0.75 Å
cutoff agreed to numerical precision in both directions during this review.
This sampled check found no changed classifications, but the harness should
use the same declared protocol before being treated as production validation.
Evidence: [water_direction_probe.json](water_direction_probe.json).

## Advice for the next development cycle

1. **Repair reliability and reporting first:** archive duplicate saturation,
   cache identity, metric names/denominators, actual-method reporting, and
   generated conclusions. Keep the `auto → HDBSCAN → agglomerative` operational
   default and the named alternatives while those repairs are made.
2. **Define three evaluation targets:** exact/near duplicate removal, broader
   geometry neighborhoods, and resource-limited representative selection.
   Use separate cutoffs, labels, and reports. Include known rigid/permutation
   duplicate witnesses, verified distinct structures, and an explicit uncertain
   category for every domain. The existing synthetic distance corpus is useful
   for invariance tests, but it cannot establish optimized basin identity.
3. **Use suitable method candidates:** for molecular conformers compare SOAP
   and MBTR with graph/symmetry RMSD; for molecular aggregates compare all-atom
   fragment matching, local SOAP/REMatch, MBTR, and averaged SOAP; for atomic
   clusters compare SOAP/MBTR and validated permutation-aware structural
   comparisons, treating inferred metallic adjacency cautiously. Compare
   scaling/normalization choices along with the descriptor, since the current
   reported result is for column-standardized Euclidean distances.
4. **Compare algorithms at explicit budgets:** HDBSCAN, average-linkage
   agglomerative, DBSCAN and OPTICS, plus the automatic policy and experimental
   controls. Measure unique-reference coverage, low-energy coverage, selected
   count, cost, fallback rate and uncertainty. Tune HDBSCAN's density parameters
   as part of development, not on the final holdout; their influence is described
   in the [HDBSCAN parameter guide](https://hdbscan.readthedocs.io/en/latest/parameter_selection.html).
5. **Expand by independent systems:** additional molecular families and sizes,
   LJ sizes beyond LJ13, metal sets with consistent calculation settings, water
   sizes below and above six, and mixed/flexible noncovalent aggregates. Hold out
   entire molecules, cluster sizes or source families; splitting dependent pairs
   from the same structure across development/test sets would overstate transfer.
   [Grouped cross-validation](https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation-iterators-for-grouped-data)
   describes the relevant separation principle. Report uncertainty by independent
   system, not by pretending every pair is an independent experiment.

SOAP's conformer macro AUC of 0.967 versus MBTR's 0.951 makes it a candidate
for the next conformer comparison. It is not evidence for changing every
system's feature default. The water and atomic observations favor different
settings under different objectives. A stronger benchmark design should precede
a universal feature or threshold policy.

## Review coverage and limits

The counts below refer to the four analytical reports (conformer, atomic,
water, meta-analysis), or the three benchmark source/reuse paths where stated.
They record observed problems in that scoped inventory, not a percentage of
scientific correctness or proof that every underlying result was audited.

| Artifact quality category | Observed defects | Assessment |
|---|---|---|
| Usefulness and completeness | 0 / 4 | Pilots disclose their limited domains; full policy choice remains an unmet follow-up question |
| Analytical clarity | 1 / 4 | Meta-analysis uses an ambiguous false-merge label and overemphasizes specificity |
| Visual/interaction consistency | N/A | Static Markdown and CSV/JSON outputs; no interactive dashboard was reviewed |

| Analytical category | Observed defects | Assessment |
|---|---|---|
| Source authority/confidence | 0 / 3 | Source snapshots/provenance inspected; external calculations and Au13 protocol uniformity not independently revalidated |
| Value accuracy | 0 / 4 | All 45 conformer confusion counts and water AUC/operating points independently agree; conformer AUC rerun agrees; atomic ranking inspected, not every coverage cell recomputed |
| Within-chart agreement | N/A | No charts in these reports |
| Complete source details | N/A | No dashboard source-detail surfaces; local source/artifact links were inspected |
| Cross-artifact consistency | 1 / 4 | Meta policy says incomplete mappings remain separate, but verified iRMSD can confirm a match after incomplete graph enumeration |
| Data-quality controls | 2 / 3 | Conformer and water reuse paths fail input-identity checks; no evidence the current conformer runs were stale |
| Conclusion support | 1 / 4 | Specificity-only operating-point guidance cannot establish safe deduplication; the domain-specific ranking caveats remain appropriate |

At review time, focused verification passed **255 tests and 38 subtests** across
clustering, distances, deduplication, conformers, basin memory, benchmark
scripts, and workflow state. The follow-up below adds regression coverage for
the defects found. Passing tests do not validate the scientific policy.

## Implementation follow-up

The follow-up changes in this branch address the demonstrated reliability and
reporting defects:

| Finding | Follow-up | Evidence |
|---|---|---|
| Basin-memory duplicate saturation | Index candidates by composition and invariant pair distances; discard only after complete graph-first comparison at 1e-5 Å. Compact verified duplicates before capacity trimming; uncertain cases remain. | `tests/test_basin_memory.py::test_memory_compacts_verified_rigid_and_permuted_copies_before_capacity` |
| Stale conformer optimization cache | Require matching ordered symbols, coordinates, xTB executable bytes/path, calculation arguments, normal termination, and optimized-geometry checksum. | `tests/test_scientific_clustering_benchmark.py::test_optimization_cache_is_invalidated_by_changed_geometry` |
| Stale water pair cache | Bind resumed pairs to source and sampled-geometry hashes, implementation hash, comparator settings, valid indices and complete finite rows. | `tests/test_water_cluster_similarity_benchmark.py::test_water_pair_cache_checks_source_and_comparator_identity`; full 276-pair rebuild |
| Misleading pair metric name | Report false-positive rate and false-discovery rate as separate quantities; update generated CSV/JSON/report. | `tests/test_meta_analyze_clustering_benchmarks.py::test_false_positive_rate_and_false_discovery_rate_use_distinct_denominators` |
| Environment-dependent feature fallbacks/report prose | Feature AUC comparison disables substitution, records unavailable features, and derives its top-feature statement from computed results. | `meta_analysis/report.md` and regenerated machine-readable tables |
| Machine-specific paths in saved runs | Run manifests use source basenames and checksums instead of this checkout's absolute paths. | Inspect saved run manifests |

The scientific limits remain: current conformer cutoff floors were not changed
without stronger removal-safety labels; benchmarked cluster co-assignment is
not equivalent to deletion; independent max-min/energy/random controls,
production-pipeline validation, explicit named-method fallback controls,
larger aggregates, and mixed-topology distance handling remain follow-up
experiments. The `auto` default remains HDBSCAN with agglomerative fallback.

## Reproduce this review

From the repository root:

```bash
.venv/bin/python benchmarks/clustering_scientific/review_20261002/benchmark_probe.py
.venv/bin/python benchmarks/clustering_scientific/review_20261002/comparison_probe.py
.venv/bin/python benchmarks/clustering_scientific/review_20261002/policy_probe.py
.venv/bin/python benchmarks/clustering_scientific/review_20261002/water_direction_probe.py
```

The first script includes deliberate mocked cache fixtures; the second and
third include synthetic coordinate fixtures for storage and fallback checks.
They are reproduction evidence, not additions to the scientific benchmark
population. Saved source results are read only. Scripts and outputs remain
here rather than in temporary storage.
