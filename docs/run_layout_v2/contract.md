# Proposed v2 contract freeze for architecture review

This is a concrete proposal, not implemented behavior. Approval freezes these
names and ownership rules; a later milestone must record and seek review for a
conflicting live-code finding. Existing v1 directories are never moved or rewritten.

## Common roots, optional files and mappings

Retain current public roots and output-selection semantics. `layout.json` is
the only v2 marker. A `state.json` version alone never identifies layout v2.
No `state.json` is invented for scan-bond; conformer retains completion/provenance
state without gaining resume. All roots may have immutable `inputs/` snapshots,
derived `summary.json`/`summary.csv` and `results/manifest.json` when applicable.
CLI audit records keep `pyar-run.toml` and append-history suffixes.
Reports and manifests never establish completion or restart validity.

`NNN` is zero-padded, at least three digits, not a maximum index. Display indexes
are persisted and never recomputed from energy rank. The stable identity is an
opaque UUID stored in metadata; directory/basename spelling is presentation only.
Stage and job labels remain their existing names. Scratch files are retained
under a native execution tree, never copied into `results/`.

### Grow pilot

```text
grow/
  layout.json
  request.json
  state.json
  pyar-run*.toml
  inputs/{seed.xyz,monomer.xyz,metadata.json}
  steps/
    step_000/selected/selected_NNN.xyz
    step_NNN/
      selection.json
      selected/selected_NNN.xyz
      attempts/attempt_NNN/execution/
        <unchanged add_one trial files, job_*, selected/, basin registry, diagnostics>
  summary.json
  results/{structure_NNN.xyz,manifest.json}
```

Seed is step 000; additions start at 001. A committed `selection.json` contains
the complete selected-reference pool, decision metadata and scientific request
identity. An incomplete step retries its existing uncommitted attempt so native
jobs can be reused; a retry after a committed attempt requires a new attempt and
cannot overwrite it. The workflow controls this decision, not a generic engine.
`steps/step_NNN/selected/` is published once on step completion and thereafter
immutable. Final result copies never become restart seeds.

| Baseline | Fresh v2 | Ownership |
|---|---|---|
| `step_000/selected` | `steps/step_000/selected` | immutable seed reference |
| `step_NNN` native execution | `steps/step_NNN/attempts/attempt_NNN/execution` | existing executor cwd/names |
| `step_NNN/selected/selected_NNN.xyz` snapshots | `steps/step_NNN/selected/selected_NNN.xyz` | committed survivor pool |
| `final/selected_NNN.xyz` | `results/structure_NNN.xyz` | immutable final copy |
| `final/summary.json` | `summary.json` | rebuildable report |

This explicitly amends the source plan's literal `{candidates,jobs}` split:
splitting shared growth execution without a workdir API would change cache
discovery. Artifact roles can distinguish candidates/jobs within `execution/`
without renaming the native files. Review this amendment before PR 1/2.

### Microsolvation

```text
microsolvation/
  layout.json
  request.json
  state.json
  pyar-run*.toml
  inputs/{solute.xyz,solvent.xyz,metadata.json}
  surface/initial_surface.json
  steps/
    step_000/selected/selected_NNN.xyz
    step_NNN/
      coverage.json
      selection.json
      selected/selected_NNN.xyz
      attempts/attempt_NNN/execution/candidates/
        <candidate XYZ, job_<candidate>/result_<candidate>.xyz, native scratch>
  summary.json
  results/{structure_NNN.xyz,manifest.json}
```

Map old `step_NNN` into `steps/step_NNN`, `final` to final result copies, and
keep the native candidates/job relationship inside execution. Original solute
indices, solvent fragment identities, surface sampling parameters, shell filter
and selected coverage are scientifically authoritative committed metadata.
`surface/initial_surface.json` is written once, not regenerated during resume.
Legacy solvate is excluded.

### Aggregate

```text
aggregates/
  layout.json
  request.json
  state.json
  pyar-run*.toml
  inputs/{component_NNN.xyz,metadata.json}
  pathways/pathway_NNN/
    pathway.json
    steps/step_NNN/
      selection.json
      selected/selected_NNN.xyz
      attempts/attempt_NNN/execution/<native ag_*, trial, job_*, selected trees>
  final_selection.json
  summary.json
  results/{aggregate_NNN.xyz,manifest.json}
```

`pathway.json` records stable pathway identity and exact addition order/composition;
not an energy-ranked label. The baseline `ag_*` trees are retained within native
execution, not parsed for v2 scientific identity. Final cross-path clustering
consumes explicitly committed per-path survivor references rather than a broad
filesystem glob; retain old glob/import logic for v1. `final_selection.json`
commits the existing final selection decisions and survivor origins. No redundant
stoichiometry nesting for a single target. Existing pathway sampling and selection
remain unchanged. Old `selected/<stoichiometry>/result_*.xyz` corresponds to new
final copies, not a newly defined selection procedure.

### Scan-bond

```text
scan_bond/
  layout.json
  request.json
  pyar-run*.toml
  inputs/{fragment_a.xyz,fragment_b.xyz,metadata.json}
  orientation_NNN/
    scan_commit.json
    start.xyz
    attempts/attempt_NNN/execution/
      scan/<native scan files, profile, trajectory, final frame>
      ts_candidates/<native candidate files and metadata>
      relax/<native relaxation attempts>
      result_relaxed.xyz
    reaction_path/
      inputs/{reactant.xyz,product.xyz,ts_guess.xyz,metadata.json}
      attempts/attempt_NNN/execution/<native cumulative stage output>
      stages/<stage>/attempts/attempt_NNN/commit.json
    result.json
  summary.{json,csv}
  results/manifest.json
```

Keep `orientation_NNN` directly under the root: the optional `orientations/`
level from the source plan is rejected here to limit consumer churn. Scan profiles
and native candidate files move only as an entire executor subtree for fresh v2;
the record supplies exact refs to consumers. Do not create another directory of
copied candidates. Manifest role `scan-maximum-candidate` explicitly disclaims
TS validation. The native directory spelling `ts_candidates` remains available
to the native execution/recovery code.

The cumulative `run_neb` executor consumes sibling stage summaries and native
artifact paths. PR 5 therefore keeps each full continuation attempt together;
per-stage immutable commits index that native execution instead of relocating
its dependent files. Changed continuation settings create a different attempt.
This follows the same executor-subtree rule as grow and avoids changing native
stage hashes/reuse. Standalone stage CLIs remain unchanged.

`scan_commit.json` commits scan identity, artifacts/profile and validation.
`through` and downstream convergence options are separate continuation controls;
extending `scan` to `all` preserves the scan commit and its bytes. Free relaxation
belongs to its own committed attempt and must not change scan identity. Raw
scan trajectories are indexed, not placed in results. A validated TS or endpoint
can be copied to results only with its successful source-stage role and gates.

### Reaction discovery and characterization

```text
reaction/
  layout.json
  request.json
  state.json
  pyar-run*.toml
  inputs/{reactant_a.xyz,reactant_b.xyz,metadata.json}
  checkpoints/geometries/<pending pool snapshots>
  discovery/
    cycles/cycle_NNN/cycle.json
    cycles/cycle_NNN/attempts/attempt_NNN/execution/
      <native gamma_*/orientation*/job_*, reaction_trace/, candidate_ts/, release evidence>
    products/product_NNN/{structure.xyz,product.json}
    commit.json
  characterization/
    pathways.json
    products/product_NNN/pathway_NNN/
      pathway.json
      inputs/{reactant.xyz,product.xyz,ts_guess.xyz,metadata.json}
      stages/<stage>/attempts/attempt_NNN/{commit.json,execution/<native stage output>}
  summary.json
  results/{product_NNN.xyz,manifest.json,<validated path structures when available>}
```

Old `gamma_*` and `trial_geometries` become native discovery execution contents;
old restart geometry snapshots correspond to `checkpoints/geometries`. Old
`products/` correspond to committed discovery-product records; old
`pathways/product_NNN/route_001` correspond to characterization path records.
No existing bytes or downstream hashes are rebased.

Cycle metadata records actual gamma/controller decisions; directory numbers do
not encode their floating-point values. Products have stable IDs independent of
discovery/rank order. Preserve current one-selected-route-per-product behavior;
`pathway_NNN` is extensible naming, not authorization to search more routes.
`discovery/commit.json` is published once after discovery completes and includes
physical settings, accepted product IDs, source trace hashes, unbiased release
evidence and relevant immutable scientific records. Characterization hashes
this commit instead of mutable progress state. Characterization never modifies it.
Future report regeneration cannot invalidate the discovery commit.

### Conformer

```text
conformers/
  layout.json
  state.json
  pyar-run*.toml
  inputs/<native converted inputs plus metadata.json>
  generated/{seed_NNNNNN_conf_NNNN.xyz,records.json}
  refinement/execution/{job_<name>/result_<name>.xyz,<native scratch>}
  completion.json
  summary.{json,csv}
  results/{conf_NNNN.xyz,manifest.json}
```

Map old `input/` to `inputs/`, `rdkit/` to `generated/`, `backend/` to
`refinement/execution/`, `selected/` to `results/`. Preserve seeded conformer
identity, torsion source identity, selection reason, RDKit force-field energy
and physical backend energy separately. `conf_NNNN.xyz` is a display rank;
stable artifact ID is not derived from it. `completion.json` records the existing
completion/selection outcome; `state.json` remains inspection metadata, not a new
resume API. Benchmark selected-file discovery must dispatch on layout version
and consume the manifest, rather than accidentally finding no structures.

### Legacy solvation and standalone stage CLIs

`solvation/state.json`, `solvation/state/geometries/`, caller-root sibling
`aggregate_NNN/` and their selected/native job subtrees are unchanged. No v2 marker
is written for legacy solvate. `neb_run/`, `ts_run/`, explicit IRC run directories
and other standalone utilities retain their present format in these milestones;
a future separately approved extension may opt them into layout v2.

## Concrete metadata proposal

These are minimum required fields. Workflow-specific scientific payloads remain
owned by their workflows; no generic state-machine class is proposed.

| File/type | Required fields | Mutability and hash ownership |
|---|---|---|
| `layout.json` | `format: "pyar-run-layout"`, `layout_version: 2`, UUID `run_id`, canonical `workflow` | immutable; independent of workflow/state schemas |
| `request.json` (supported workflows) | `request_schema_version: 1`, `run_id`, `workflow`, `scientific_request`, `inputs` (artifact refs) | immutable once created; continuation-only controls excluded |
| ArtifactRef | `artifact_id` UUID, `role`, `path` run-relative POSIX string, `sha256` lowercase 64 hex digits, optional `source_record_id` UUID | committed reference immutable; no absolute path, `..`, symlink escape or cwd rebasing |
| State | `workflow_state_schema_version`, `run_id`, `workflow`, exact resolved request identity, progress/status, committed references, uncommitted active attempt | atomic mutable progress; each workflow owns validation |
| Scientific completion record | `stage_result_schema_version: 1`, UUID `record_id`, `run_id`, `stage`, UUID `attempt_id`, `physical_method`, `scientific_settings`, input/output ArtifactRefs, dependency record IDs and hashes, `status`, `validation` | successful committed record immutable; failures remain separate attempts; gate result explicitly recorded |
| `results/manifest.json` | `format: "pyar-results-manifest"`, `manifest_schema_version: 1`, `run_id`, `workflow`, `artifacts` | derived index; each entry has ArtifactRef, presentation rank if meaningful, successful source record ID and role; rebuildable |
| `inputs/metadata.json` | `input_schema_version: 1`, `run_id`, source label, saved ArtifactRefs, resolved charge/multiplicity/scftype where relevant | immutable snapshots, no machine-specific source path required |
| JSON/CSV summary | own `report_schema_version: 1` for JSON; workflow payload and source record IDs | derived; never hash report bytes as v2 scientific dependency |

Proposed new workflow-state schema versions: grow 2, microsolvation 2, aggregate
2, reaction 3, conformer 4. Scan-bond has scan request/commit records, no generic
state. Legacy schemas and their `version` fields remain exactly as today.
Final schemas are separately dispatched through the explicit layout marker;
equal numeric versions across kinds carry no implication of compatibility.

Canonical JSON for hashed v2 records: UTF-8, sorted string keys, compact
separators `(',', ':')`, finite JSON numbers (`allow_nan=False`), `ensure_ascii=False`,
one trailing LF. Hash actual published bytes using SHA-256. No hash of a pretty
report or a reconstructed object substitutes for the recorded file hash.
Energy records include value, unit and convention, e.g. physical potential
`{"value": -1.0, "unit": "hartree", "convention": "physical_potential"}`;
RDKit MMFF/UFF values retain their force-field label/native kcal/mol convention.
Geometry-only energies are null, never zero. Biased objective values are separate
and never named barriers. Preserve current numerical precision/XYZ writer when
publishing structures; changing scientific serialization requires separate review.

Artifact IDs are assigned once when an artifact is created and persisted before
reuse; source molecule/conformer identity is also recorded explicitly. Manifests
may change display rank while retaining IDs. Rebuilding a manifest uses committed
records and verifies hashes; it cannot manufacture successful scientific gates.

## Publication, detection and locking

Fresh construction is side-effect-free. `--check`/`--dry-run` read/validate only:
no lock, directory, input copy, marker or state writes. Execution acquires a
single-writer OS advisory lock on a sibling `.<root-name>.pyar-v2.lock` before
checking/creating a fresh root. Keeping the lock outside the root avoids changing
empty-directory detection and serializes first creation. Restrict support to
filesystems on which the chosen OS lock has tested semantics; never unlink the
lock file to recover a held lock. Process death releases the OS lock; stale owner
metadata may be replaced only after acquiring it, with an explicit recovery
option and recorded recovery event. No active lock can be force-broken by the CLI.
No cross-host/distributed locking guarantee or thread-safe `os.chdir` guarantee
is implied. PR 1 must validate the narrowly scoped lock before the pilot.

Execution first snapshots inputs and writes valid metadata. Missing marker plus
nonempty recognized v1 records dispatches to the unchanged v1 owner; no marker is
added. Recognition is workflow-specific: scan-bond uses its request signature
and native result/recovery records, without requiring `state.json`. Other
workflows use their established state/request contracts. Directory spelling or
an arbitrary JSON file alone does not establish a recognized run.
Nonempty unrecognized output is rejected. Unknown version, malformed
marker, marker/workflow/run-ID mismatch or mixed v1/v2 authority fails closed.
Interrupted initialization with a valid marker but incomplete required metadata
fails with an actionable explicit recovery/fresh-output instruction; it is not
silently interpreted as a fresh run.

For new scientific stages: write temporary artifacts in the same filesystem,
flush/fsync needed bytes, atomically publish, verify hashes, then publish the
immutable completion record, and finally atomically update mutable progress.
Flush parent directories for durable publication where supported. Crash before
completion leaves an uncommitted attempt; crash after completion but before
progress update may recover only by validating the complete record and artifacts.
Completion/gate validation belongs to the workflow, never mere file presence.
Reports/manifests are published last and may be rebuilt. Never replace a
committed artifact path during retry. Native backends can update files only in
their active uncommitted execution tree.

The current NEB engine writes cumulative mutable stage summaries into one output
directory. Adapting it to isolated v2 attempts and immutable scientific records
is a **PR 5/6 requirement**, not a change to be bundled into primitives/grow.
Dependencies must carry forward exact validated upstream structures/settings
and gates; do not copy arbitrary old summaries and rewrite their hash keys.

## CLI/API compatibility acceptance

V1 result values and inspect record selection stay unchanged. Fresh v2 results
explicitly return paths to the new final structures; preserve the caller-visible
absolute/relative convention of each owning API, using the documented run-root
base internally. This path change must appear in workflow docs and its PR's
before/after tree, not be hidden by basename-symmetric renaming. `state_path`
points to state only where it already exists; scan-bond retains its dict result
shape/`output_dir`. `pyar inspect` must add layout-aware reading without treating
the derived manifest as authoritative or changing latest CLI-record precedence.

Every owning milestone must update hard-coded CLI, summary and benchmark
consumers alongside its v2 writer, with v1/v2 fixtures for both. Existing roots,
backend physics, preflight dependency assumptions, legacy aliases, CLI record
fallbacks, and standalone utilities remain compatible. No implementation
milestone is authorized by PR 0 until architecture review accepts this contract.
