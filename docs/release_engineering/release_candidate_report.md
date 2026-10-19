# Release candidate engineering report

**Candidate:** master `4c94f87217cf9a85af29077884fe2efefcf0fc20` (19 October
2026). **Package metadata:** 1.5.1. **Decision:** NOT READY. No final release
artifacts were frozen and no publication or tag was made.

## Engineering checks

PR #78 added optimization cache backend provenance and merged after exact-head
review. Cache identity includes backend/version provenance when safely available
(xTB executable/version, AIMNet2 model hash and Torch version, hybrid provenance,
and geomeTRIC version when selected). If the version is unknown or its probe
fails, cache reuse/write is disabled while calculation can proceed. Focused
local tests on the PR implementation: 43 passed; Ruff and `git diff --check`
passed. The complete relevant local suite on the same implementation reported
1,743 passed, 3 skipped, 717 warnings and 194 unittest subtests. Exact-head CI
passed the supported Python 3.10–3.13 matrix, installation combinations,
documentation and R1 required-check aggregate:
[PR run 39192224563](https://github.com/anooplab/pyar/actions/runs/39192224563),
[push run 39192219869](https://github.com/anooplab/pyar/actions/runs/39192219869).

The post-merge local full-suite attempt under system Python 3.14.7 showed many
import/fixture failures and was interrupted; it is not supported-matrix evidence.
The first isolated local build and Sphinx attempts hit the temporary filesystem quota. Retried with build and output temporary paths on the home filesystem: strict Sphinx completed with zero warnings; `python -m build` created wheel and sdist; `twine check` passed both artifacts. These exploratory artifacts identify as version 1.5.1 and were built from an uncommitted docs worktree, so they are not final release artifacts. A fresh local install attempt stalled on dependency downloads and was stopped; CI clean installation checks passed.

CI built and clean-installed wheel and sdist and tested core plus optional
combinations. These are CI checks, not frozen release artifacts. No approved version 2.0
wheel/sdist pair or uploaded artifact exists. The exploratory 1.5.1 artifact hashes
were wheel `c020c3a1905f8b1e678342b3394139bbdadd232e4bfbd511a6707a7970c09c58`
and sdist `33376d6990437fa50786fe33b3b2165e5d1e5b633bbaac599f0f873746b3a906`;
they are local at `/home/anoop/pyar-release-artifacts/` and are not eligible for
publication.

## Scientific result

[Scientific Tier 2 run 39192735168](https://github.com/anooplab/pyar/actions/runs/39192735168)
passed at the exact candidate SHA. It records xTB 6.7.1 GFN2, Open Babel 3.2.1,
Python 3.11.17 and PyAR 1.5.1. The artifact includes candidate SHA, input and
request data, scan profile/trajectory, structured result and reactant identity.
Three H₂ scan energies were finite and the final constrained distance was
0.7400000005 Å. The scan maximum is explicitly not a confirmed TS and no product
was registered. The workflow does not preserve separate backend stdout/stderr.
The complete provenance and scientific-owner review table are in
[scientific validation final](scientific_validation_final.md).

## Publication and recovery procedure

Before any future publication: freeze the exact candidate, release notes,
migration/backend compatibility notes, approved scientific corpus and thresholds;
build wheel/sdist once, inspect contents, run `twine check`, clean-install those
exact bytes and record SHA-256; verify mandatory gates and obtain owner approval
of version, SHA, hashes and registry; publish only those bytes; then verify
registry hashes and run clean registry install/CLI checks.

If publication partially fails, stop promotion, preserve logs/artifacts, inspect
registry state and obtain an owner-approved recovery decision. Never move a
published tag or replace uploaded bytes. No recovery action applies to this
candidate because nothing was published.
