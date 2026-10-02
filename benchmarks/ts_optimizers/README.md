# Transition state optimizer benchmark

This directory documents the paired benchmark runner. It does not contain
third-party dataset files. The first comparison uses the same starting
geometry, charge, multiplicity, electronic structure model, force threshold,
and step ceiling for geomeTRIC and Sella.

## RGD1-TSopt-GFN2 pilot

The adapter supports the external [RGD1-TSopt-GFN2 dataset](https://doi.org/10.5281/zenodo.20489312)
(version 1.0). The Zenodo record describes 1,500 distinct neutral, closed-shell
C/H/N/O transition states and 4,500 deterministic starting geometries across
easy, med, and hard tiers. Each reference saddle, start, and reaction endpoint
is stored in XYZ; `manifest.tsv` gives atom counts and displacement metadata.
The reference implementation is tblite 0.6.0 at SCF accuracy 0.01. PyAR runs
the normal `xtb` executable with `--gfn 2`; matching the model family does not
claim exact numerical identity between tblite and `xtb`.

The Zenodo record states copyright and does not give an open redistribution
license. The adapter reads an extracted local copy and writes only a manifest
with paths, hashes, and source metadata. Do not commit the dataset archive or
structures into this repository.

Extract the dataset locally, then make the default 15-reaction pilot (45
starting geometries, with each reaction represented at all three difficulty
tiers):

```bash
pyar-benchmark-ts prepare-rgd1 \
  --dataset-dir /path/to/RGD1-TSopt-GFN2 \
  --output /path/to/RGD1-TSopt-GFN2/pyar_pilot.json
```

Selection is deterministic and stratified over four atom-count ranges. The
manifest records the random seed, selected reaction IDs, source dataset
manifest/settings hashes, displacement metadata, and the reference energy.
Choose 7–16 reactions for a 21–48 input pilot, or set `--tiers` to include a
subset of tiers. `--ts-fmax` and `--ts-max-cycles` are shared by both optimizers;
the adapter defaults are 0.02 eV/angstrom and 200 steps. A first paired smoke
case was nonstationary for Sella at 0.03 eV/angstrom but passed independent
frequency/IRC/endpoint checks at 0.02, so 0.02 is the provisional pilot target.
The full pilot will check it across the sampled set; this does not select a
production optimizer default.

Run the paired jobs in a deterministic shuffled order. Each optimizer/case is
launched in a separate process, and Sella gets a fresh temporary JAX compilation
cache for each run so warm-cache effects do not depend on job order:

```bash
python benchmarks/ts_optimizers/prepare/run_pilot.py \
  /path/to/RGD1-TSopt-GFN2/pyar_pilot.json \
  --output /path/to/rgd1_pilot_run \
  --schedule-seed 20261003

pyar-benchmark-ts collect /path/to/rgd1_pilot_run
```

The execution order and cache policy are recorded in `pilot_execution.json`.

## Manifest

Create a JSON manifest with:

- `name`;
- `source`: dataset name, version, license, and DOI or URL;
- `qc_model`: `software`, plus model-specific settings such as `xtb_model`;
- `settings`: required shared `ts_fmax` and `ts_max_cycles`, with optional
  frequency, IRC, and endpoint validation controls;
- `cases`: unique IDs, a `ts_guess`, reference TS, charge, multiplicity,
  preserved difficulty label, and optional reactant/product endpoints.

Geometry paths are relative to the manifest unless absolute. All geometries in
a case must use the same ordered atoms. Reactant and product must be supplied
together. xTB manifests must name `xtb_model` explicitly; use `gfn2` when the
reference data are defined on GFN2-xTB. The original manifest is copied to
`benchmark_manifest.json`. Runs record the manifest hash, input hashes, PyAR
and optimizer versions, platform, model, and available backend version
information for each execution.

Both optimizers receive the same explicit `ts_fmax` and `ts_max_cycles`. The
runner does not substitute an optimizer after a failure. It preserves every
stage artifact in the case/optimizer directory. When endpoints are supplied,
it also runs the existing relaxation, frequency, IRC, and endpoint validation
stages. Those costs are reported separately from TS optimization.

## Run one paired job

Run each case/optimizer independently; this works with shell parallelism and
scheduler arrays without embedding scheduler code:

```bash
pyar-benchmark-ts run benchmark.json \
  --case case_0001_medium \
  --optimizer geometric \
  --output benchmark_run

pyar-benchmark-ts run benchmark.json \
  --case case_0001_medium \
  --optimizer sella \
  --output benchmark_run
```

Each job writes `cases/<case-id>/<optimizer>/benchmark_result.json` alongside
the unmodified PyAR stage summaries, trajectories, and copied input geometries.
Existing or interrupted job directories are never overwritten; use a fresh
output directory to repeat a job.

## Collect

```bash
pyar-benchmark-ts collect benchmark_run
```

Collection writes `runs.csv`, `runs.json`, `paired.csv`, and `summary.json`.
Missing jobs appear as `incomplete_run`; they are not treated as failures or
successful optimizer outcomes. Paired differences are emitted only when both
runs use the manifest's matching starting-geometry hash. The summary is
descriptive; it does not select a default optimizer.

## Outcome definitions

- `optimizer_not_converged`: optimizer returned without convergence.
- `converged_not_stationary`: frequency validation found a nonstationary result.
- `stationary_not_first_order_saddle`: stationary point did not have exactly
  one significant imaginary frequency.
- `validated_first_order_saddle`: frequency confirmed a first-order saddle and
  endpoint information was unavailable.
- `first_order_saddle_wrong_connection`: saddle validation passed but supplied
  endpoint validation did not confirm the requested reaction connection.
- `reaction_connected_success`: saddle, IRC, and endpoint checks passed.
- `*_exception`: the named calculation/validation stage raised an exception.

Optimizer success is not inferred from a return code alone. When endpoints are
available, the primary success is the complete reaction connection. Optimizer
cost and frequency/IRC/endpoint validation cost remain separate.
