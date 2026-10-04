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

The task-oriented CLI pilot supports clustering:

```bash
pyar clustering *.xyz
pyar clustering --help
```

`pyar-clustering *.xyz` continues to work with the same options, defaults,
and implementation. The shell expands `*.xyz` into input paths. Other
workflows continue to use their existing commands.

Optimize every supplied XYZ independently with the same computational settings:

```bash
pyar optimize *.xyz --backend xtb
pyar optimize *.xyz --backend xtb --check
pyar optimize *.xyz --backend orca --method r2scan-3c --nprocs 8
```

`--check` validates all inputs, electronic states, options, and requirements
without starting calculations or creating job directories. A backend is always
explicit. Charge defaults to zero; omitted multiplicity is inferred as singlet
for even electron counts and doublet for odd counts, separately for each input.
xTB uses GFN2-xTB without a basis; ORCA/Gaussian/Turbomole use BP86/def2-SVP
unless overridden (ORCA built-in composite methods supply their own basis).
Processors default to one and the geometry optimizer to native (geomeTRIC for
Gaussian, whose existing native adapter emits a single-point route); geomeTRIC is
available for backends with an existing energy/gradient route.

`pyar-optimiser` remains available with its original options and defaults.
The modern command performs minimum optimizations only and omits reaction
bias controls (`--gamma`, `--site`). Options the selected adapter cannot consume
are rejected rather than silently ignored. Distinct input basenames are required
because they identify independent optimization jobs.

The modern bond scan uses two XYZ fragments and 0-based, fragment-local atom
indices, with eight starting orientations by default:

```bash
pyar scan-bond A.xyz B.xyz --atoms 3 4 --backend xtb
pyar scan-bond A.xyz B.xyz --atoms 3 4 --backend orca --method BP86 --basis def2-SVP
pyar scan-bond A.xyz B.xyz --atoms 3 4 --backend xtb --through all
pyar scan-bond A.xyz B.xyz --atoms 3 4 --backend xtb --through all --check
```

The backend is explicit; `--orientations` (or `-N`) overrides the orientation
count. Standalone xTB defaults to GFN2-xTB; `--xtb-model` selects another
supported model. ORCA retains its existing method/basis and external g-xTB
wrapper validation. A single charge or multiplicity applies to both fragments;
two values apply separately. Omitted multiplicities follow electron parity.
The default `--through scan` stops after the scan. Later stages keep the existing
continuation controls and output hierarchy (`--output scan_bond`).

All inputs, atom indices, settings and stage-dependent requirements are checked
before orientation generation. `--check` performs this validation without
launching programs, creating output directories, or running any stage. It checks
dependency availability and importable APIs, rather than probing executables.
`pyar-cli scan-bond` remains supported with its original `--software` option,
ORCA default, and required `-N` argument.

- `aggregate` for clusters, aggregates, and noncovalent complexes
- `react` for AFIR-style reaction searches between two reactants
- `solvate` for microsolvation, ligand addition, and growth around a core
- `scan-bond` for a simple bond-distance probe
- `pyar-reaction-trace` for reaction-trace analysis

## External Program Requirements

Some workflows rely on external executables such as xTB, ORCA, Gaussian, Psi4, MOPAC, Turbomole, OpenBabel, MLatom, and DFT-D4. The optional extras install Python dependencies only; they do not bundle large model or vendor files into the main `pyar-chem` wheel. AIMNet2 `.jpt` models, AIQM1 `.pt` models, and vendored MLatom binaries must come from the upstream project or another separate model/package source. See [docs/external_programs.rst](docs/external_programs.rst) for the official project websites and installation notes.

## Documentation

Full documentation: [docs/](docs/) and https://pyar.readthedocs.io/
Changelog: [CHANGELOG.md](CHANGELOG.md)

## Citation

If you use PyAR, cite the paper that matches your chemistry problem. See [docs/publications.rst](docs/publications.rst) for the current publication map. For general cluster-building use, start with:

- Nandi et al., *Computational and Theoretical Chemistry* 1111, 69-81 (2017)
- Khatun et al., *Frontiers in Chemistry* 7:644 (2019)

### Inspect stored energies and compare structures

```bash
pyar energies *.xyz
pyar compare reactant.xyz product.xyz
pyar compare A.xyz B.xyz --atom-mode all
pyar compare A.xyz B.xyz --json
pyar energies *.xyz --json
```

`energies` ranks stored XYZ comment-line energies from lowest to highest, using
Eh for absolute energies and kcal/mol relative to the global minimum. Every
input must contain a readable energy. `pyar-energy-table` remains supported.

`compare` reports composition, stored energies (when available), graph-aware
aligned RMSD, component counts, and inferred edge changes. Its energy sign is
ΔE = E(B) − E(A). Energy differences between different compositions are not
relative isomer/conformer energies. The default RMSD uses heavy atoms, with
all atoms used for hydrogen-only systems. No default geometry-equivalence
threshold is applied; `--rmsd-threshold` requests an explicit threshold.

Connectivity uses an element-labelled coordinate graph with covalent-radius
adjacency (default `--bond-scale 1.15`), not formal bond orders. Indexed edge
changes use 0-based indices and assume row correspondence when element order
matches; XYZ cannot verify correspondence of repeated elements. When element
order differs, indexed changes are omitted while graph-isomorphism and
permutation-aware RMSD remain available. `--maximum-mappings` controls the
existing enumeration limit, and incomplete comparisons are reported explicitly.
Both commands read files without running calculations. The legacy
`pyar-similarity` pool deduplication command remains available unchanged.
