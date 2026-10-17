# Repository Guidelines

## Project Structure & Module Organization

PyAR is a Python package distributed as `pyar-chem` while preserving the `pyar`
import namespace. Core source lives in `pyar/`: CLI wiring in `pyar/cli.py`,
console helpers in `pyar/scripts/`, backend adapters in `pyar/backends/`, sampling
logic in `pyar/sampling/`, bundled model/runtime assets under `pyar/AIMNet2/`,
the external MLatom compatibility adapter in `pyar/mlatom/`, and workflow/state
code in the corresponding subpackages.
Tests live in `tests/` as `test_*.py`. Sphinx documentation lives in `docs/`.
Avoid editing generated artifacts such as `build/`, `pyar_chem.egg-info/`,
`__pycache__/`, and local log files.

## Build, Test, and Development Commands

- `python -m pip install -e ".[test,docs]"`: install the package in editable
  mode with test and documentation tools.
- `python -m build`: build source and wheel distributions from `pyproject.toml`.
- `pytest`: run the pytest suite.
- `python -m unittest discover -s tests -v`: alternate full test discovery used
  in prior repository checks.
- `python -m sphinx -E -b html -W --keep-going docs docs/_build/html`: build docs
  and fail on warnings.
- `pyar-cli --help` and entry points such as `pyar-react --help`: smoke-test CLI
  registration after packaging changes.

## Coding Style & Naming Conventions

Use Python 3.10+ syntax and 4-space indentation. Prefer small, explicit
functions over broad side effects, especially in CLI and backend code. Keep
module names lowercase with underscores. Classes use `CamelCase`; functions,
variables, and tests use `snake_case`. `pyproject.toml` is the authoritative
packaging surface; keep `setup.py` as a compatibility shim only.

## Testing Guidelines

Use pytest-style files named `tests/test_<feature>.py`. Add focused tests near
the changed behavior, for example CLI tests in `tests/test_cli.py` or packaging
metadata checks in `tests/test_packaging_metadata.py`. For workflow tests that
would call external chemistry programs, patch the boundary and leave direct
negative-path tests for fail-fast validation.

For scientific refactors, compare behavior under identical inputs, settings,
pathway plans and random state. Equal RNG seeds alone do not control iteration
over sets of identity-hashed Molecule objects. Control that ordering in the
comparison harness rather than changing production sampling to make a test pass.
Check selected structures, ordering, energies and composition as applicable;
geometry-only checks do not establish live external-backend correctness.

## Workflow Persistence Changes

Use `docs/run_layout_v2/source_plan.md` and `contract.md` for the output-layout
series, and check live GitHub state before deciding the next milestone.
Preserve scientific defaults, pathway sampling, selection, backend arguments
and native cache discovery. Keep v1 continuation with its existing owner; do
not silently migrate existing runs. Fresh v2 runs use explicit layout markers,
relative hashed artifact references and immutable selection records. Final
result copies are not restart authority.

Exercise interrupted publication and recovery, input/request mismatch, artifact
tampering and durable v1 fixtures. Shared helpers must preserve the caller's
workflow-specific exception contract, including resume paths. For frozen legacy
functions, an AST comparison is useful alongside runtime compatibility checks.

## Authorized PR Delivery

For an authorized implementation series, complete one PR at a time: implement,
review, fix, review again, validate, then merge only when authorized and all
applicable CI checks pass for the current head. Reassess the next milestone
after each merge. Do not infer merge permission from a review-only request.
Use the installed `pyar-pr-delivery` skill for the detailed delivery procedure.

Use a repository-provisioned environment and retain raw validation logs outside
the tracked tree. Report final exit status and test counts, not partial progress
output. Run focused regressions and appropriate full-suite, strict Sphinx and
package checks. Inspect Sphinx's generated-file side effects and exclude them
from the commit unless requested; preserve any pre-existing user edits.

Respect an explicit stop-on-failure instruction: preserve the branch/worktree,
inspect the actual traceback, stop delivery and later milestones, and send the
authorized notification. A subsequent instruction to continue fixing that PR
authorizes repair iterations within that PR; it does not permit merging a
failing head or automatically restarting later milestones. Never send mail
without authorization in the task/session. Verify send success before claiming
notification, and report inability to send if the connector fails.

Verify remote merge state and merge commit, then fast-forward a clean local
main branch and confirm it matches the remote. A passing local suite or a
successful push alone does not establish delivery.

## Commit & Pull Request Guidelines

Recent commits use short imperative summaries such as `Align packaging metadata`
and `Expand xTB backend guide with installation checks`. Keep commits scoped and
avoid mixing scientific workflow changes with packaging or documentation cleanup.
Pull requests should describe the behavior changed, list commands run, mention
external backend assumptions, and link issues when applicable.

## Security & Configuration Tips

Do not commit machine-specific paths, credentials, generated logs, or large new
runtime assets without provenance and license notes. Backend executables such as
`orca`, `g16`, `psi4`, `mopac`, `xtb`, `define`, and OpenBabel tools are external
system dependencies; document assumptions instead of hard-coding local paths.
