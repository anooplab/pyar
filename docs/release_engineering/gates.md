# Mandatory PyAR release gates

Policy established 2026-10-18; updated for R3 on 2026-10-19. These gates govern a future public release, not
ordinary feature completion. They do not reopen the frozen ORE scientific scope
or turn its separate benchmark project into unfinished implementation work.
The new release-engineering R0–R5 stages are distinct from the ORE R0–R5 series.

This is the canonical gate definition. The [preparation skill](../../.agents/skills/prepare-release/SKILL.md)
supplies the procedure; the [plan](plan.md) identifies remaining automation.
R3 adds a deterministic scientific evidence comparator and Tier 1/Tier 2
execution surfaces. It does not freeze or approve the release corpus.

## Evidence and decision rules

- Each required gate has status `PASS`, `FAIL`, `BLOCKED` or `NOT_RUN`.
  Only executed evidence can yield `PASS`. Missing tools, credentials, licensed
  backends, baselines or validators mean `BLOCKED` or `NOT_RUN`, never success.
- Record the candidate Git SHA, proposed version, exact command, environment,
  dependency/backend/model versions, input/reference hashes, seeds/thread counts,
  start/end timestamps, exit code, structured results and raw-log locations.
  Record wheel/sdist filenames and SHA-256 hashes. Retain failed attempts too.
- A zero process exit is necessary but insufficient when structured output
  contains scientific failure, partial completion or unsupported cases. Require
  both successful execution and the gate's explicit numerical/schema checks.
- Freeze the required job list, backend matrix, benchmark cases and numeric
  tolerances before evaluation. Never move thresholds after seeing results.
  Do not invent universal chemical tolerances; absent approved values block the
  scientific gate. Benchmarks use tolerances and controlled seeds, not assumed
  bitwise reproducibility across hardware.
- Any candidate or artifact change invalidates its approval and applicable gate
  results. Reused reference evidence needs explicit provenance; old CI results,
  local `dist/` files and historical audit prose do not qualify a new candidate.
- Preparation reports gate evidence, remaining failures and **NO-GO** until all
  mandatory checks pass. It may report **READY FOR OWNER APPROVAL** afterward;
  only an explicit approval of the exact candidate permits publishing.

## Gate matrix

| Gate | Mandatory acceptance and evidence | Existing execution surface / remaining gap |
| --- | --- | --- |
| G0 Candidate | Clean isolated candidate; intended diff reviewed; version/import/distribution metadata agree; changelog and supported-scope statement match the candidate. No unrelated files or unapproved scientific changes. | Git commands, `pyproject.toml`, `CHANGELOG.md`, metadata tests. Unified candidate manifest/version/changelog validator is missing. |
| G1 Code and regressions | Full supported Python CI matrix passes on the exact SHA, including current coverage floor. Relevant negative, interruption/replay, v1/v2 compatibility and optional-dependency tests execute. Review skipped/xfail cases against an approved issue inventory; new unexplained omissions block. False-success and objective-corruption defects affecting release claims block until corrected or the affected capability is explicitly excluded from those claims. | `.github/workflows/ci.yml`, pytest and coverage configuration. Automated skip/xfail inventory and aggregate required-check validator are missing. |
| G2 Scientific and backend validation | Run the preapproved representative corpus on the candidate: optimization, scan/reaction, conformer/selection, microsolvation and TS/IRC routes as selected in the release support matrix. Compare against approved references/baseline within declared energy/gradient/geometry/identity tolerances; verify units, electronic state, objective, constraints, budgets and failure classification. TS claims require saddle/frequency, both IRC branches and endpoint connection checks. Include invalid/stale/interrupted negative controls and paired cost evidence for recovery/performance claims. | Tier 1 is mandatory in PR CI. `.github/workflows/scientific-tier2.yml` runs real xTB GFN2 scan smoke weekly and manually. `tools/scientific_validation.py` checks manifest integrity, hashes, methods, electronic states, versions, raw-log references and numeric tolerances. The catalogue has no owner-approved release corpus; required Tier 2/Tier 3 cases and thresholds remain a release blocker. A small smoke run or mocked solver cannot satisfy the whole gate. |
| G3 Documentation | Strict Sphinx passes; user-facing commands, installation/extras/backend instructions, changelog and version boundaries match actual behavior. Modern and legacy examples have appropriate automated help/preflight/smoke coverage. Generated documentation is reviewed without silently discarding unrelated edits. | Existing `Documentation` CI job and CLI tests. A systematic executable example/CLI-to-doc contract and consolidated release-note check are missing. |
| G4 Installation and optional compatibility | Fresh non-editable core, declared workflow/extra compositions, wheel and sdist installs run outside the checkout; `pip check`, import-origin checks and real CLI/behavior smoke checks pass. Missing extras, models and executables fail actionably before work. Advertised backend versions have explicit positive/negative evidence or a clearly disclosed, owner-approved unvalidated/excluded status in the support matrix; optional absence cannot be relabelled a positive calculation. | `tools/installation_smoke.py`, installation tests and CI matrices exist. Full final-artifact/Python/extra coverage, dependency constraint snapshots and licensed-backend evidence are incomplete. |
| G5 Distribution | Build fresh wheel and sdist from the frozen candidate with recorded build tools; archive inspection and `twine check` pass; required runtime files/licenses exist; private data, generated caches, unlicensed datasets and external model weights are excluded. Both artifact hashes are recorded. Install and test these exact artifacts; promotion must upload the approved bytes. | `python -m build`, `tools/check_distributions.py` and clean-install CI exist. The publish workflow rebuilds independently and lacks the approved-artifact promotion contract and full metadata/license checks. |
| G6 Approval and publish safety | G0–G5 evidence is complete; project owner explicitly approves version, SHA, artifact hashes and target registry. Repository/environment restrictions and required checks are verified. Publication cannot be triggered as a side effect of preparation, generic merge approval or a green CI result. No silent overrides of failed mandatory gates. | `pypi` environment and OIDC Trusted Publishing exist. The observed environment has no protection rules; publishing lacks aggregate release-gate enforcement. These are implementation blockers, not permissions granted by this document. |
| G7 Post-release verification | After separately authorized publishing, version-specific registry metadata contains the expected wheel/sdist hashes; clean registry installs and supported CLI/doctor/smoke checks pass, release notes/tag/commit agree, and logs and artifacts are retained. Publication failure uses a documented incident decision; do not force-move a released tag or replace published bytes. | Manual version-specific registry verification has historical precedent; a deterministic post-release verifier and incident procedure are missing. This gate runs after publication and cannot be claimed in advance. |

G0–G5 are prerequisites to requesting final release approval; G6 is required
before tag push, publish dispatch or upload. G7 closes the release after the
authorized publication. A preparation task stops at the approval boundary.

## Backend and scientific scope

The release support matrix must list each advertised workflow/backend/optimizer
combination, supported Python and backend/model versions, required assets and
licenses, which environment executes it, expected negative behavior, and whether
positive scientific validation is required. Exclusions must be explicit before
evaluation and reflected in user documentation and the owner's approval.
Do not quietly drop a failed mandatory case or describe a route as validated
because its executable was discovered or its Python module imported.

Scientific methods, thresholds and target semantics remain owned by the approved
scientific plan. Release engineering measures that contract; it does not choose
new methods, weaken chemistry gates or switch calculators to make tests pass.

## Release evidence record to implement

R1–R4 define and validate a machine-readable record containing:

- candidate SHA/version and scope; artifact hashes and build environment;
- the required gate/job/case matrix, thresholds and approved exclusions;
- per-gate status, exact commands, exits, run IDs, raw and structured artifacts;
- test/coverage/skip/xfail summaries and benchmark/reference provenance;
- owner approval identity/time, target registry and approved hashes;
- post-publication registry/install results and any incident record.

The R3 scientific comparator is not the aggregate release-record validator.
Missing mandatory cases remain blockers; agents must report the gaps rather than
substitute a subjective "looks ready" verdict.
