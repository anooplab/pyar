# PyAR

PyAR is a chemistry-focused structure-search package for aggregation, reaction discovery, solvation growth, and bond scans.

## Install

```bash
python -m pip install pyar-chem
```

## Quick Start

```bash
pyar-cli --help
pyar-cli -a C H -as 1 4 -N 8
pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16
```

## Supported Workflows

- `aggregate` for clusters, aggregates, and noncovalent complexes
- `react` for AFIR-style reaction searches between two reactants
- `solvate` for microsolvation, ligand addition, and growth around a core
- `scan-bond` for a simple bond-distance probe
- `pyar-reaction-trace` for reaction-trace analysis

## External Program Requirements

Some workflows rely on external executables such as xTB, ORCA, Gaussian, Psi4, MOPAC, Turbomole, OpenBabel, MLatom, and DFT-D4. See [docs/external_programs.rst](docs/external_programs.rst) for the official project websites and installation notes.

## Documentation

Full documentation: [docs/](docs/) and https://pyar.readthedocs.io/

## Citation

If you use PyAR, cite the paper that matches your chemistry problem. See [docs/publications.rst](docs/publications.rst) for the current publication map. For general cluster-building use, start with:

- Nandi et al., *Computational and Theoretical Chemistry* 1111, 69-81 (2017)
- Khatun et al., *Frontiers in Chemistry* 7:644 (2019)
