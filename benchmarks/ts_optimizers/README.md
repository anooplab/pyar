# Transition state optimizer benchmark

This directory documents the paired benchmark runner. It does not contain a
dataset or benchmark results yet. The first comparison uses the same starting
geometry, charge, multiplicity, electronic structure model, force threshold,
and step ceiling for geomeTRIC and Sella.

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
