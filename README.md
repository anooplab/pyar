# PyAR
PyAR stands for "Python program for aggregation and reaction"

Documentation:

[![Documentation Status](https://readthedocs.org/projects/pyar/badge/?version=latest)](https://pyar.readthedocs.io/en/latest/?badge=latest)

* Local source: [`docs/`](docs/)
* Hosted docs: https://pyar.readthedocs.io/

# Installation

The PyPI distribution name is `pyar-chem`. The import namespace remains
`pyar`, and the main command-line entry point is `pyar-cli`.

From PyPI:

```bash
python -m pip install pyar-chem
```

From a local checkout, install the package for development with:

```bash
python -m pip install -e .
```

This installs the `pyar` package, the `pyar-cli` command line tool, and the
bundled runtime assets used by the MLatom and AIMNet2 interfaces.

Optional Python extras install feature-specific dependencies:

```bash
python -m pip install "pyar-chem[selection]"
python -m pip install "pyar-chem[ml]"
python -m pip install "pyar-chem[aimnet2]"
python -m pip install "pyar-chem[openbabel]"
python -m pip install "pyar-chem[test]"
python -m pip install "pyar-chem[docs]"
python -m pip install "pyar-chem[all]"
```

For a local checkout, use the same extras with `-e`, for example
`python -m pip install -e ".[test,selection]"`.

## Python API

PyAR 2.0 has a narrowed public Python API. Use `import pyar`,
`pyar.workflows`, `pyar.sampling`, `pyar.selection`, and `pyar.backends` for
supported imports. Legacy paths such as `pyar.trial_generation`,
`pyar.orientation_sampling`, `pyar.aggregator`, `pyar.reactor`,
`pyar.interface`, and `pyar.data_analysis.clustering` are not supported in
2.0; see `docs/migration_2_0.rst` for the import mapping.

# Features:
* Automated prediction of unknown reactions between two reactants (A+B)
* Automated prediction of the geometries of aggregates, atomic clusters etc.
* Automated search for reaction for bond forming between two atoms in two different molecules.



## Setting Up Environment

To set up your environment for the tasks, you can create and edit your `.bashrc` or `.bash_profile` file, depending on your system configuration. After that source ~/.bashrc to run those changes.

```bash
# Create an alias for "mndo2020"
alias mndobin="mndo2020"

# Add Gaussian 16 to your PATH
export PATH=$PATH:/home/apps/g16

# Create an alias for GAUSS_EXEDIR
alias GAUSS_EXEDIR="g16"

# Create an alias for MLatom.py
alias mlatom="MLatom.py"

```

```bash
# Install dftd4 executable in this way
conda config --add channels conda-forge
conda install dftd4
conda install -c conda-forge openbabel
alias dftd4bin="dftd4"
```

## Requirements 
Core Python dependencies are declared in `pyproject.toml` and installed
automatically with `python -m pip install pyar-chem`. The default install is
kept small enough for `import pyar`, `pyar-cli --help`, basic molecule helpers,
and sampling/trial-generation utilities.

Optional extras:

* `[selection]` for DScribe/MBTR descriptors, HDBSCAN, MDAnalysis,
  matplotlib, pandas, scikit-learn, and ASE-backed clustering/descriptor
  helpers
* `[ml]` for MLatom, TorchANI, PyTorch, H5MD, and related ML interfaces
* `[aimnet2]` for AIMNet2 runtime support with PyTorch and ASE
* `[xtb]` for Python-side geomeTRIC/ASE support used by xTB-style reaction
  optimization paths
* `[openbabel]` for the OpenBabel Python binding (`openbabel>=3.2.0`)
* `[test]` for pytest, build, and twine
* `[docs]` for Sphinx documentation builds
* `[all]` for the combined optional Python dependency set

### External programs required by backend

External executables are not Python dependencies and must be installed
separately for the corresponding backend:

* ORCA: `orca`
* Gaussian: `g16`
* Psi4: `psi4`
* MOPAC: `mopac`
* xTB: `xtb`
* Turbomole: `define` and related Turbomole tools
* OpenBabel tools: `obabel`, `babel`, `obminimize`, `obenergy`
* Dispersion correction utilities: `dftd4` when using workflows that call it

# Interfaced with electronic structure theory programmes
- mlatom_aiqm1
- aimnet2
- Mopac
- Turbomole
- Psi4
- Xtb
- Orca

# Molecule generations 

```bash
pyar-cli -a C H -N 8 -as 6 6 --software aiqm1_mlatom -m 1 2
```

You can also generate from a formula in aggregate mode:

```bash
pyar-cli --aggregate --formula C5H4 -N 8
```

# Molecular clusters

## XTB
```bash
pyar-cli -s water.xyz water.xyz --software xtb -ss 10 -N 16 -c 0 0 -m 1 1
```
## AIMNet2
```bash
pyar-cli -s water.xyz water.xyz --software aimnet_2 -ss 10 -N 16 -c 0 0 -m 1 1
```

PyAR bundles AIMNet2 model assets for the AIMNet2 interfaces. AIMNet2 is a
third-party project from the Isayev Lab and is MIT licensed upstream. See
`THIRD_PARTY_LICENSES/AIMNet2-LICENSE` and
`THIRD_PARTY_LICENSES/AIMNet2-PROVENANCE.md` for license and provenance details.

This will generate molecules up to a maximum of 6 carbon and 6 hydrogens with **aiqm1_mlatom** potential using 8 trial orientations.
Here `C` and `H` are element-symbol inputs. XYZ files are still accepted when you want to provide explicit starting coordinates:
```bash
1
carbon
C  0.0  0.0   0.0
```
```bash
1
hydrogen
H  0.0  0.0   0.0
```
# Reaction

To study the reaction between two reactants A and B using ORCA software interface, with force from 100 to 1000 using N=8 trial orientation, the commandline line argument is,

```bash
pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software orca
```

A.xyz and B.xyz are the cartesian coordinate files of the reactants.

After a successful geomeTRIC-backed AFIR reaction run, analyse one orientation
job directory with:

```bash
pyar-reaction-trace reaction/gamma_0100/orientation_xxxxx --plot
```

or, from inside the orientation directory:

```bash
pyar-reaction-trace . --plot
```

This writes `path_summary.csv`, candidate geometries in `candidate_ts/`, and
optional plots in `trace_plots/`.

For `pyar-cli`:

- `--react` requires exactly two XYZ input files.
- `--scan-bond` requires exactly two XYZ input files.
- `--solvate` requires at least two XYZ input files.
- `--formula` is only supported together with `--aggregate` in this CLI.

## pyar-cli
The main program can be used as below:

```
pyar-cli options files
```

There are other scripts for a few automation tasks.

* `pyar-reaction-trace` analyzes a reaction trace in a job directory or
  `reaction_trace/` directory and can write PNG plots with `--plot`.
  Use `--plot-only` to skip rewriting the summary files.
* `pyar-cli trace` runs the same reaction-trace analysis from the main CLI.

## pyar-trial-generation
`pyar-trial-generation` can be used for
* for making different orientations of two molecules.
* Making a composite molecule containing a _seed_ molecule and __N__ number of monomer molecules.
* Orient two molecules such that _i_'th atom of one molecule and _j_ 'th atom of second molecule have shortest distance
between them

## pyar-clustering
* for a clustering analysis of n input molecules to find unique molecules.

## pyar-optimiser
* for the bulk optimisation of several molecules



# References

1. "A Global Optimizer for Nanoclusters ", Maya Khatun, Rajat Shubhro Majumdar, Anakuthil Anoop <a href="https://www.frontiersin.org/articles/10.3389/fchem.2019.00644/full">Frontiers in Chemistry 2019, 644</a>
1. "A tabu-search based strategy for modeling molecular aggregates and binary reactions" S Nandi, SR McAnanama-Brereton, MP Waller, A Anoop, <a href="https://www.sciencedirect.com/science/article/pii/S2210271X17301627">Computational and Theoretical Chemistry 2017, 1111, 69-81</a>  
1. "AIMNet2: A Neural Network Potential to Meet your Neutral, Charged, Organic, and Elemental-Organic Needs" D. M. Anstine, R. Zubatyuk, O. Isayev, <a href="https://doi.org/10.1039/D4SC08572H">Chemical Science 2025, 16, 10228-10244</a>
