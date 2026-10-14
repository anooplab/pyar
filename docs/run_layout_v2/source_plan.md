# PyAR Run Layout v2: Architecture Contract and Implementation Plan

**Status:** PROPOSED — PR 0 audit and explicit approval required before implementation
**Prepared:** 2026-10-13
**Repository:** `anooplab/pyar`
**Scope:** filesystem layout, artifact references, persistence integrity, and backward compatibility
**Baseline:** capture the exact commit SHA and current tests in PR 0; do not trust this document as a substitute for live code inspection.

## 1. Objective

Make every new PyAR calculation directory understandable to a computational chemist and navigable by tools, **without changing chemistry, computational methods, selection, restart validity, or cache semantics**. Distinguish (1) the user-facing run layout, (2) backend-native execution directories, (3) scientific stage states/provenance, and (4) report indexes.

**Architectural rule:** A directory is an artifact location, not a source of scientific truth. Stable identities, requested physics, validated stage dependencies, and selection decisions belong in explicit metadata. Existing backend-specific paths may remain internal when changing them would compromise caching or execution.

## 2. Non-negotiable invariants

1. No change in candidate geometry generation, random-number sampling, spin/charge, energy/gradient conventions, AFIR/softmin biasing, optimization, deduplication, ranking, connectivity, or scientific gate logic.
2. Never automatically relocate, rewrite, or convert an existing v1 run, whether complete or partial. Retain its established restart/inspection path. Old run states that already refer to files must continue to resolve them.
3. Do not rewrite an upstream state/stage summary when its **byte hash** is used by completed downstream stages. Provide a new v2 stable provenance contract for new stages; preserve v1 logic unchanged.
4. A stage is reusable only when its recorded physical method, scientific settings, upstream dependencies, artifact hashes, and scientific validation gates pass. `scan-bond --through scan` followed by `--through all` must reuse verified scan work.
5. `react` discovery followed by NEB, TS optimization, TS frequency, IRC, endpoint relaxation, and endpoint frequency remains possible without repeating verified work. A discovered product is not a confirmed pathway or activation barrier.
6. No mutable process-level cwd assumptions may silently change. Replace global `os.chdir` separately, behind explicit executor workdir contracts; do not combine whole-backend rewrites with this migration.
7. Do not promise new resume semantics for workflows that currently lack them (notably `conformer`).
8. Never change public CLI defaults and returned paths just for symmetric spelling. Distinguish modern `pyar solvate` (alias for `microsolvate`) from deprecated `pyar-cli solvate` (legacy compatibility workflow).
9. Respect optional dependencies; layout tests must run without ORCA, xTB, AIMNet2, or RDKit unless a test specifically requires them. Use mocked providers and optional smoke tests.
10. `--check` and `--dry-run` must not create or mutate calculation output. Do not use destructive migration of user directories.

## 3. Current-code audit targets (verify in PR 0)

- `AGENTS.md`, `pyproject.toml`, `pyar/run_config.py`, `pyar/modern_cli.py`, `pyar/workflow_results.py`: coding rules, version, execution records, preflight, result paths.
- `pyar/workflows/{aggregate,grow,microsolvation,solvation,reaction,reaction_characterization,scan_bond,scan_path,conformer}.py`.
- `pyar/state/{aggregate,grow,microsolvation,solvation,reaction}.py`, plus any other workflow state owners.
- `pyar/growth/{service,reporting}.py`; `pyar/optimiser.py`: `job_<name>`, `result_<name>.xyz`, `selected/`, and cache-sensitive discovery.
- `pyar/neb.py`: stage summaries, physical settings, dependency hashes, scientific stage gates, legacy summary handling, absolute artifact paths.
- `pyar/reaction_analysis.py`: trace and candidate paths; `pyar/workflows/reaction_characterization.py`: product traces and route-specific persistence.
- `pyar/workflows/scan_bond.py`: `request.json` signature, orientation result caches, and `ts_candidates/` paths.
- `pyar/benchmarks/conformer.py`, user-facing CLI scripts, docs and tests: consumers of concrete paths.
- Current top-level names include `aggregates/`, `reaction/`, `conformers/`, `grow/`, `microsolvation/`, `scan_bond/`; legacy `solvate` uses `solvation/` with historical sibling `aggregate_NNN` work directories.

**PR 0 deliverables:** A path writer/reader/cache/test/doc matrix; exact repository commit SHA; v1 fixture inventories; end-to-end current-state tests and known failures; incompatibilities and proposed fixes; an explicit list of which paths must retain their original spellings. PR 0 must be reviewable without any workflow implementation code changes.

## 4. Filesystem layout contract

For new v2 runs, keep the *existing public default root spelling*. Within each root, use roles consistently and only when applicable:

```
<existing-default-run-root>/
  request.json                 # immutable resolved scientific identity where supported
  state.json                   # only for workflows that have restart state
  summary.json                 # derived report; not a restart authority
  summary.csv                  # optional table
  pyar-run.toml                # existing CLI audit/provenance record
  inputs/                      # immutable snapshots of original inputs
  checkpoints/                 # if restart-only geometry snapshots are needed
  <domain-stages>/             # steps, pathways, cycles, orientations, etc.
  results/                     # final user-facing structures, and manifest.json
```

- `candidates/`: generated, not accepted; `selected/`: intermediate survivors; `jobs/`: calculator executions; `backend/`: third-party input/output files; `trace/`: trajectories and scientific trace metadata; `results/`: final user-facing outputs only.
- Avoid copying raw trajectories and backend scratch files into `results/`. Final XYZ files may be **copied** there so the manifest provides stable, convenient entry points; avoid hard links that might allow accidental mutation of committed results.
- Preserve established lower-level optimizer job naming and caching (`job_<name>/result_<name>.xyz`) unless independently proven safe to change.
- Relative artifact references use one canonical run-root base, reject path traversal and unexpected symlink escapes, and include a SHA-256 hash where integrity is required.
- Artifact IDs must not depend on energy ranking, file basename, or directory parsing. A manifest is an index, not authoritative state; it should be rebuildable from committed scientific records.

### Workflows

**grow/**: `steps/step_000/selected/` holds initial seed; `steps/step_NNN/{candidates,jobs,selection.json,selected}/` retains staged additions; `results/{structure_NNN.xyz,manifest.json}` exposes final survivors. Keep snapshot paths and restart validation working.

**microsolvation/**: `surface/initial_surface.json` documents the original solute surface; `steps/step_000/` and `steps/step_NNN/` document each newly added solvent, coverage, and selection. Preserve site restrictions, original-solute targeting, post-optimization shell filtering, and geometry-only behavior. `results/` contains final selected clusters. Do not migrate deprecated `pyar-cli solvate` as part of this project.

**aggregates/**: `pathways/pathway_NNN/pathway.json` records exact addition sequence and composition; `steps/` and `selected/` represent intermediate attempts and survivors. `results/aggregate_NNN.xyz` and `manifest.json` expose final clusters. Do not depend on `ag_*` folder names for identity, but retain backend job files and compatibility lookups as necessary. Do not invent redundant stoichiometry nesting for single-target runs.

**scan_bond/**: `orientations/orientation_NNN/` may be introduced only if compatible with established CLI consumers; `scan/` retains profile, trajectory, final frame, and backend outputs; `scan_candidates/` or a new carefully named user-facing directory distinguishes scan maxima from proven TS geometries. Preserve `ts_candidates/` internal legacy lookup where required. `reaction_path/` retains dependency-validated stage summaries; do not mutate a validated stage in place during layout conversion. Keep scan request identity distinct from changeable continuation settings.

**reaction/**: separate `discovery/cycles/` from `characterization/products/product_NNN/pathway_NNN/`. Discovery contains orientations, traces, unbiased relaxation evidence and *path candidates* (not certified TS). Characterization contains NEB, TS, frequency, IRC and endpoint validation stages, with independently verifiable provenance and input hashes. A source product need not have a validated pathway. A new v2 immutable discovery commit record should be used for future scientific dependency hashes, instead of hashing mutable discovery progress state. Do not rewrite v1 `reaction/state.json` or its downstream hash chain.

**conformers/**: `generated/` (embedded/torsion sampling), `refinement/` (calculator jobs), and `results/` (ranked selected structures). Do **not** add undocumented generic resume support. Preserve matching of energy, structure, source conformer identity, and ranking.

**solvation/**: legacy implementation stays on v1. Only independent correctness/safety patches with their own tests; do not silently map v1 restart state to `microsolvation/`.

Detailed folder depth and individual filenames are a **design input** to PR 0, not a license to rename backend-native artifacts. PR 0 must freeze exact v2 trees and old-to-new mappings before code changes.

## 5. Metadata and provenance ownership

| Artifact | Source of truth for | Mutable? |
|---|---|---|
| `request.json` | resolved scientific run identity (where a workflow has this contract) | No after creation |
| `state.json` | current/recoverable workflow progress | Yes, atomic |
| Per-stage records | stage parameters, physical method, artifact hashes, dependency hashes, completion/gate state | Committed results should be immutable; retries have separate attempts |
| `summary.json`, `summary.csv` | derived human/machine-readable report | Rebuildable |
| `results/manifest.json` | indexed stable IDs, output roles, run-relative paths, checksum, source stage | Rebuildable |
| `pyar-run.toml` | existing CLI argument/config/profile/execution provenance | Append history using existing semantics |

Version independently: layout format, workflow state schema, stage-result schema, and manifest schema. Do not overload a single `version` field. JSON serialization must be deterministic when its bytes participate in a checksum. Include value + **units and energy convention** for scientific energies; never equate biased AFIR objective with physical potential energy or barrier.

## 6. Minimal shared API (not a workflow engine)

Candidate names only; PR 0 may choose existing package location if cleaner:

```
pyar/persistence/
  layout.py       # immutable RunLayout/path resolution; no I/O in construction
  artifacts.py    # run-root-relative, validated ArtifactRef and hash helper
  json_store.py   # small atomic JSON read/write primitives
  versions.py     # v1/v2 detection/dispatch; no implicit migrations
```

Do **not** implement `BaseRunState` hierarchies, a generic scientific state machine, database, plugin system, or workflow scheduler. Workflows continue to control stage transitions. Avoid expanding existing APIs until the pilot proves the abstraction needed.

## 7. Crash consistency and locking

A complete stage should be committed only after necessary outputs have been successfully written and verified. New attempts should not overwrite artifacts referenced by a committed state. Persist committed stage state atomically **after** publishing/validating artifacts; on restart ignore/recover incomplete attempt files and never treat their presence alone as completion.

Use a single-writer guard per run directory for v2 workflows; stale-lock recovery must be explicit and safe, not a blind deletion. A lock protects run consistency and does **not** imply arbitrary thread-safe `os.chdir()` execution. No concurrent mutation of one run directory should be supported in the first release.

## 8. Compatibility protocol

- v1 run: retain existing code/path handling; never automatically convert. A completed v1 run is inspectable; where v1 already supported resume/continuation, that contract continues. Do not claim features v1 never supported.
- v2 run: identified unambiguously by a versioned metadata marker. Never infer v2 solely from a `steps/` or `results/` folder.
- Corrupt/ambiguous run: fail closed with actionable error. Do not create a new run in a nonempty unrecognized directory.
- Preserve old CLI default directories through initial migration; use `--output` when supported, and migrate public defaults only through a separate approved deprecation process.
- Keep legacy compatibility `pyar-cli solvate` distinct from modern `pyar solvate`/`pyar microsolvate` semantics.
- Do not silently alter behavior of `pyar inspect`, `pyar-run.toml`, returned `WorkflowResult.run_directory`, `state_path`, or `selected_paths`.

## 9. Proposed PR schedule, one PR at a time

0. **Read-only audit and contract freeze** — matrix of path producers/consumers/hashes; exact tree and schema; baseline fixtures; no science code edits. Stop for architecture review.
1. **Minimal primitives** — layout paths, safe ArtifactRef, version detection, atomic JSON utilities, concurrency guard only as justified by pilot; no change to existing v1 outputs.
2. **grow pilot** — v2 on fresh runs, existing v1 stays readable/resumable; preserve scientific deterministic fixtures, backend naming and selected pool checksum behavior. Stop for review.
3. **microsolvate** — preserve original-solute-centred sampling and shell coverage; no legacy solvate migration.
4. **aggregate** — restructure pathways, selected survivors and final manifest; retain old reuse semantics and legacy pathway import.
5. **scan-bond** — preserve scan request identity and stage extension; ensure `--through scan` to `--through all` does not rescan.
6. **react** — coordinated discovery and characterization with explicit scientific dependency contract; verify entire chain.
7. **conformer** — isolate generated/refinement/results; no new resume claims.
8. **release** — docs, public API/CLI references, packaging, full optional-backend smoke tests, compatibility and deprecation documentation.

Each PR must cite this specification, list exact changed paths, validate that no scientific decisions were touched, and provide test commands, outcomes, known limitations, and expected follow-ups. A PR is not complete until CI passes; do not automatically start the next PR while a review gate remains open.

## 10. Verification requirements

- Freeze representative v1 fixture runs *before changing paths*: aggregate alternative pathways, interrupted grow, interrupted microsolvation, scan-bond with stages, reaction discovery+characterization, conformer completion, legacy solvate state.
- Compare fixed-seed/mocked-backend scientific outputs: coordinates within appropriate tolerances, energies, rank ordering, candidates, survivor identities, coverage, connectivity and pathway acceptance. For stochastic runs capture RNG seeds and input geometry hashes.
- Test kill/crash between artifact publish and state commit, duplicate invocation to same root, stale/partial results, corrupt state, altered scientific request, unknown layout version, path traversal/symlinks, missing artifact and checksum mismatch.
- Verify v1 continuation and v2 stage extension. A stage hash reused for NEB/TS/IRC must not be invalidated merely by reports being regenerated.
- Verify read-only preflight and current CLI run-record semantics, `pyar inspect`, returned API result paths, postprocessing and benchmark file consumers.
- Enforce filesystem contracts on fresh v2 runs but do not assume all backend scratch files are already confined unless executor code has been checked; failures should be explicit and documented.
- Test with optional backends absent. Add an opt-in smoke suite with real xTB and ORCA where available; no requirement for external programs in universal unit tests.
- Use repository-standard checks from `AGENTS.md`: pytest, Sphinx with warnings as errors, and build; report exact versions and test commands. Do not invent passing results.

## 11. Codex execution rules

- Read root `AGENTS.md`, this document, and affected workflow docs/tests before making changes.
- For PR 0: strictly read-only on Python workflow code; write audit/spec findings only. Flag any inconsistency, update draft layout decisions, then stop for human review.
- For each implementation task: select exactly one numbered PR; create an isolated branch/worktree; keep code changes scoped; write minimal tests first where feasible; run targeted tests and all feasible global checks; include before/after output trees and scientific equivalence evidence in PR description.
- If observed code conflicts with this plan, **do not invent a workaround or widen scope silently**. Record the conflict and propose a reasoned amendment for approval.
- Do not bulk rename files, change chemical algorithms, alter optional dependencies, introduce new frameworks, change legacy storage, or begin the next milestone without approval.
- Save approved spec amendments and PR links in a small progress section (or linked GitHub issue), rather than keeping crucial decisions only in chat.

## 12. PR 0 approval checklist

Before accepting PR 0, confirm: the audit includes every filesystem path consumer and producer; the exact immutable/mutable file contracts are frozen; complete vs intermediate artifact semantics are explicit; v1/modern alias distinctions are correct; stage dependency and hash rules preserve known continuation; crash/locking design is narrowly scoped; CLI/API outputs and optional install paths are protected; deterministic fixtures and test commands are identified; no speculative scientific claim has been added.

**Status:** Not implemented. No test outcomes, compatibility claims, or GitHub PRs are implied by this document.
