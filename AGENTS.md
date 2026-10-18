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

- `python -m pip install -e ".[selection,xtb]" --group dev`: install the CI
  test capabilities and developer tools using pip with dependency-group support.
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
packaging surface; do not reintroduce a second metadata source in `setup.py`.

## Testing Guidelines

Use pytest-style files named `tests/test_<feature>.py`. Add focused tests near
the changed behavior, for example CLI tests in `tests/test_cli.py` or packaging
metadata checks in `tests/test_packaging_metadata.py`. For workflow tests that
would call external chemistry programs, patch the boundary and leave direct
negative-path tests for fail-fast validation.

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

## Permanent Engineering Policies

- Prevent regressions with focused tests and the applicable full suite. Preserve
  legacy APIs and continuation contracts; report skips and expected failures.
- Preserve the approved scientific scope, units, objectives and acceptance
  criteria. Check raw backend and workflow evidence; process exit alone cannot
  establish scientific success. Scientific changes need explicit owner approval.
- Synchronize user-facing behavior, CLI examples, compatibility/version notes
  and changelog entries with implementation. Validate documentation commands
  and strict Sphinx output for the affected scope.
- Keep optional dependencies optional. Test supported backend/extra combinations
  and actionable missing-dependency failures; distinguish discovery/import tests
  from calculations using real executables and model assets.
- For packaging changes and release qualification, verify fresh wheel and sdist
  artifacts in clean, non-editable installations outside the checkout. Check
  dependency consistency, imports, entry points, required files and
  excluded/private/generated assets.
- Release preparation follows [.agents/skills/prepare-release/SKILL.md](.agents/skills/prepare-release/SKILL.md)
  and the [release gates](docs/release_engineering/gates.md). Record executed
  evidence against the exact candidate and artifact hashes. Missing mandatory
  evidence blocks readiness. Passing checks never grants permission to tag,
  dispatch publishing or upload packages; obtain explicit owner approval of the
  version, commit and artifacts before those actions.
