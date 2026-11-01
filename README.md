# PyAR

PyAR searches molecular structures and reaction candidates, runs electronic-structure calculations through selected backends, and helps inspect the resulting evidence. Use one command: `pyar`.

Python 3.10–3.13 is the classified, CI-tested range. The distribution is `pyar-chem`; the Python namespace is `pyar`.

## Start here

Read the [installation and first-calculation guides](https://pyar.readthedocs.io/). They take you from a fresh environment to an energy calculation, geometry optimization, minimum qualification and result interpretation. The user guide covers conformers, aggregates, growth, explicit microsolvation, reaction exploration, scans, NEB/TS/IRC, and analysis. Developer/API details are separate.

These sources describe **1.6 development**, with an intentionally new CLI and data compatibility boundary. Preserve the original environment/version for older work. No automatic migration of older runs is provided.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "pyar-chem[selection,xtb,conformer] @ git+https://github.com/anooplab/pyar.git@master"
pyar --help
pyar doctor
```

For reproducibility pin a source commit instead of `master`. Python extras do not install external executables such as xTB or ORCA; install the selected backend separately and validate the actual route with `--check`.

```bash
pyar energy water.xyz --backend xtb --check
pyar optimize water.xyz --backend xtb
pyar frequency job_water/result_water.xyz --backend xtb
```

Download the water input and follow the complete worked example in the documentation. A successful optimizer does not by itself establish a minimum, saddle or reaction connection. Retain raw backend outputs, settings, units and scientific qualification with reported results.

## Development

```bash
python -m pip install -e ".[selection,xtb]" --group dev
python -m pytest -q tests/test_optimizer_input_state.py tests/test_backend_trial_input.py
python tools/check_documentation.py
python -m sphinx -E -b html -W --keep-going docs docs/_build/html
```

See the developer guide and `AGENTS.md` for contribution and release procedures. PyAR is distributed under GPLv3; retain component licenses and consult the publication references in the documentation when citing its chemistry workflows.
