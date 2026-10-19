# PyAR 1.5.2 candidate and publication report

## Release identity and provenance

- Version: `1.5.2` (`pyar-chem`)
- Exact tagged candidate: `f621275c1ef19950bc562264f1319112eb3bc6ad`
- Versioned owner-approval record: `docs/release_engineering/releases/1.5.2.json` (merged before promotion)
- Approval-record merge commit: `a19e43e38f460299a3321fc85d27277dc4ee72a8`
- Post-release verifier merge commit: `7383f9140d82fef3b4b6d27667e3410498e7efca`
- PyPI: <https://pypi.org/project/pyar-chem/1.5.2/>
- GitHub release and attached files: <https://github.com/anooplab/pyar/releases/tag/v1.5.2>
- Exact candidate CI run `39211560723`: <https://github.com/anooplab/pyar/actions/runs/39211560723>
- Exact candidate Tier 2 run `39211573722`: <https://github.com/anooplab/pyar/actions/runs/39211573722>
- Promotion run `39212474899`: <https://github.com/anooplab/pyar/actions/runs/39212474899> (verification and publish jobs passed; its first post-job was canceled after the version-specific PyPI JSON shape issue was identified).
- Corrected post-publication verification run `39213822582`: <https://github.com/anooplab/pyar/actions/runs/39213822582> (PASS).

## Distribution artifacts

| File | SHA-256 | Published registry status |
|---|---|---|
| `pyar_chem-1.5.2-py3-none-any.whl` | `c63927e0aaa74d91b8d17239b71b57028e2fac6e901974383f90fac75b4a244b` | Present on PyPI; exact digest verified. |
| `pyar_chem-1.5.2.tar.gz` | `e0fedb6878f438ba0882f9d91049b69e8e9e1d70657e3398f0d6a823fd0d127a` | Present on PyPI; exact digest verified. |

The same files are attached to the GitHub release and their asset digests match. Candidate CI artifact ID is `12168250384` (archive SHA-256 `b56933a264ad30cd084fa2cdf0f638c48d7838f9508b0da52e9801c816f17042`). Final scientific Tier 2 artifact ID is `12167911527` (archive SHA-256 recorded in the approval record).

## Automated checks

The exact candidate CI run completed successfully with 28 successful jobs, including the required R1 aggregate, Python 3.10–3.13 test matrix, documentation, lint/format and fast regression checks, wheel/sdist builds and installations, Sella/geodesic interface checks, and 18 optional-extra installation jobs across Python 3.10 and 3.13.

Test evidence artifacts report:

- Full pytest suite: 1,948 tests on each of Python 3.10, 3.11, 3.12 and 3.13; zero failures and errors; 13 optional-dependency skips per version (`irmsd` and `geodesic-interpolate` unavailable).
- Fast regressions: 312 tests, zero failures/errors and one skip because that step does not have the built wheel fixture.
- Deterministic scientific Tier 1: 79 tests, zero failures/errors/skips.
- Distribution metadata: `twine check` passed locally for both final artifacts.
- Post-release clean install: GitHub Actions Python 3.12 and a local Python 3.14 environment both installed `pyar-chem==1.5.2`, passed `pip check`, confirmed `pyar` imports from the fresh venv, and passed `pyar-cli --help` and `pyar --help`.

## Release notes and scientific scope

The release addresses backend optimization result validation, stale markers and cache provenance, per-job failure isolation, geometry restoration on failure, ORCA HF input, per-input reporting and atomic configuration/profile writing. User-facing CLI and installation documentation passed the strict documentation job.

Scientific claims are limited to the owner-approved maintenance scope recorded in `releases/1.5.2.json` and described in `scientific_validation_final.md`. The broader scientific benchmark corpus and major-release validation remain incomplete. In particular, no broad claim is made for reaction exploration, NEB, transition-state/frequency/IRC validation, aggregation/conformer ranking, microsolvation or non-xTB real backends.

## Publication outcome

The guarded promotion workflow verified the exact candidate SHA, successful exact-SHA CI and Tier 2 runs, job inventories, owner approval and artifact hashes before publication. The PyPI environment is restricted to the `master` branch and administrator bypass is disabled. Tag `v1.5.2` points to the exact candidate commit. The version-specific PyPI JSON API now returns the two matching files; the corrected verification-only workflow and clean installation checks passed.
