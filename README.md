# PyAR

PyAR is a chemistry-focused structure-search package for aggregation, reaction discovery, solvation workflows, and bond scans.

## Install

```bash
python -m pip install pyar-chem
```

## Quick Start

```bash
pyar-cli --help
pyar-cli -a C H -as 1 4 -N 8
pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
pyar microsolvate solute.xyz solvent.xyz --count 10 --backend xtb
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

- `aggregate` searches a requested final composition with alternative build pathways
- `grow` repeatedly adds one species to a specific seed and retains every selected stage
- `react` for AFIR-style reaction searches between two reactants
- `microsolvate` for original-solute-centred first-shell sampling
- legacy `solvate` for backward-compatible solvation restarts
- `scan-bond` for a simple bond-distance probe
- `pyar-reaction-trace` for reaction-trace analysis

## Modern conformer, aggregate and grow commands

```bash
# RDKit conformer search, optionally followed by backend refinement
pyar conformer "CCO"
pyar conformer "CCO" --backend xtb

# Requested final composition; inputs are building-block types, not a fixed seed
pyar aggregate A.xyz B.xyz
pyar aggregate water.xyz --size 6 --backend xtb
pyar aggregate water.xyz ammonia.xyz --size 4 2 --backend xtb

# Fixed seed, repeated addend; every selected intermediate stage is retained
pyar grow metal.xyz ligand.xyz --count 4 --backend xtb
pyar grow cluster.xyz H.xyz --count 12

# First-shell microsolvation keeps targeting the original solute surface
pyar microsolvate methane.xyz water.xyz --count 8
pyar microsolvate methane.xyz water.xyz --count 8 --backend xtb
pyar microsolvate solute.xyz water.xyz --count 3 --site 7 --backend xtb

# Validate, save and inspect reusable task setups
pyar optimize molecule.xyz --backend xtb --dry-run
pyar react A.xyz B.xyz --backend xtb --bias-max 100 --write-config reaction.toml --check
pyar run reaction.toml --check
pyar info
pyar backends
pyar doctor

# Read-only input, settings, requirements and restart validation
pyar conformer "CCO" --check
pyar aggregate A.xyz B.xyz --check
pyar grow seed.xyz ligand.xyz --count 4 --check
```

All three commands allow omission of `--backend`: conformer uses its established
RDKit search, while aggregate and grow use bounded structural diversity selection
without fabricated energies. Aggregate explores requested final compositions and
alternative build orders, with no permanently privileged seed. Grow repeatedly
adds one species to a fixed seed and retains its growth series; it is not limited
to solvation. Legacy `solvate` remains separate.

Aggregate and grow default to eight orientations (`--orientations`/`-N`), with
survivor budgets of eight and twelve respectively. Aggregate `--size` needs one
positive count per input; omission means one of each for multiple inputs. A
single building block requires an explicit size. Existing element/formula inputs
and `--formula` expansion retain their coordinate/atomic-cluster semantics.

Backend settings follow the modern minimum-optimization policy: standalone xTB
uses GFN2 without DFT method/basis defaults. Charge defaults to neutral for XYZ;
omitted multiplicity follows electron parity. Conformer uses a reliable input
formal charge when available. Explicit charge/spin overrides are validated.
RDKit remains optional: install `pip install "pyar-chem[conformer]"` for conformer
search. XYZ conformer input also needs the external OpenBabel `obabel` executable.
`--check` validates XYZ coordinates and converter availability without invoking
OpenBabel conversion, embedding, optimization or writing state.

Existing output layouts and restart policies remain authoritative: `conformers/`
(requires a fresh run), `aggregates/` and `grow/` (or grow's `--output DIR`).
`pyar-conformer`, `pyar-cli conformer`, `pyar-cli aggregate`, `pyar-cli grow` and
legacy boolean workflow flags remain supported with their existing contracts.

`grow` places repeated addends against the complete current structure. The new
`microsolvate` workflow instead samples the original solute's accessible surface
through first-shell construction; prior solvents block occupied regions and
remain in the collision/energy environment but do not become new placement
targets. Geometry-only microsolvation is bounded and does not invent energies.
Legacy `pyar-cli solvate` remains on its existing restart-compatible path and
prints a deprecation warning. See [microsolvation documentation](docs/microsolvate.rst).

## Aggregate and Grow

```bash
# Final-composition search, without a permanent privileged seed
pyar-cli aggregate A.xyz B.xyz --aggregate-size 2 3 -N 16 --software xtb

# Sequential additions to a specified seed, with bounded survivors at each step
pyar-cli grow metal.xyz ligand.xyz --count 4 -N 16 --software xtb --xtb-model gfn2
```

Omit `--software` for bounded geometry-only growth. Intermediate stages, resolved
request, restart state, and final summary are retained under `grow/`.
See [growth documentation](docs/grow.rst) for selection and restart policy.

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

`compare` also attempts optional chemical bond-order perception with RDKit's
maintained `rdDetermineBonds` implementation of the xyz2mol approach, followed
by canonical isomeric SMILES. This is separate from coordinate adjacency and
RMSD. SMILES depend on the geometry and the charge used; XYZ alone does not
encode charge, multiplicity, bond orders, or complete chemical identity.

```bash
pyar compare anion_a.xyz anion_b.xyz --charge -1
pyar compare neutral.xyz anion.xyz --charge 0 -1
```

One `--charge` value applies to both structures; two values apply to A and B.
Without this option, the XYZ interface assumes neutral charge and reports it
as **assumed**, rather than known. The reusable molecular API can also accept
reliably known charges from other input/state. PyAR never searches other
charges to force successful perception. Perception failures retain their
reason while coordinate and energy analysis continue. Canonical SMILES are
model-dependent evidence; matching SMILES do not establish complete chemical
or stereochemical identity.

RDKit remains optional. To enable SMILES perception, install the existing
extra with `pip install "pyar-chem[identity]"`. Without RDKit, the report shows
chemical perception as unavailable and continues with all coordinate-based
analysis. JSON includes each perception result, charge and source, failure
reason, and canonical-SMILES match status.

### Composable XYZ utilities

```bash
# Generate eight deterministic encounter geometries (no optimization)
pyar orient A.xyz B.xyz
# Inspect a unique geometry set without writing files
pyar deduplicate *.xyz
# Energy window first, then graph-first duplicate removal
pyar select *.xyz --within 5 --unique
# Inspect atoms, Hill formula, energy, components, and optional SMILES
pyar identify *.xyz
# Extract disconnected coordinate components
pyar split complex.xyz
# Use the existing reaction-trace analysis and plotting
pyar trace reaction_job --plot
```

`orient` writes `orientations/` with `orientation_000.xyz`, etc., and the exact
`trial_vectors.dat`. It uses the established Fibonacci directions, Halton
quaternion rotations (disabled for atomic incoming fragments), and contact
placement with scale 1.2 for atoms or 1.5 for molecules. Use `--orientations`
(or `-N`), `--sequence-offset`, `--distance-scale`, and `--output` to override.

`deduplicate` uses PyAR's existing adaptive RMSD threshold and conservative
graph-first policy. It retains lower-energy representatives when every energy
is available; otherwise input order is preserved. Failed or uncertain
comparisons retain structures. `--rmsd-threshold` provides an explicit override.

`select` requires `--within` (inclusive kcal/mol window) or `--top`. With both,
it applies the window, then top N, then optional `--unique` pruning. Every
input must contain a readable energy. Both `select` and `deduplicate` inspect
only by default; `--output DIR` copies retained originals. Non-empty destinations
and colliding input basenames are rejected.

`identify` accepts one `--charge` for all inputs or one per input and supports
`--json`. Missing charge is labelled assumed neutral. It shares RDKit
DetermineBonds chemical perception with `compare`, without charge searching.
RDKit is optional; `pip install "pyar-chem[identity]"` enables SMILES for both
commands. Missing RDKit or failed perception leaves coordinate information
available.

`split` uses covalent-radius adjacency (`--bond-scale 1.15`) and preserves
original atom order and coordinates. It writes `<stem>_fragments/` or `--output
DIR`; one component produces no redundant file. Components are geometric
adjacency groups, not formal chemical fragments. No energies, charges or spin
states are assigned. Geometry-only outputs explicitly mark energies unavailable.

`orient` and `split` refuse non-empty destinations. Inputs are never modified.
`trace` retains legacy artifacts and options including `--plot`, `--plot-only`,
`--plot-directory`, `--max-force`, and `--exclude-energy-outliers`; an omitted
path defaults to the current directory. `pyar-reaction-trace` remains supported.

### Modern reaction search

```bash
pyar react A.xyz B.xyz --backend xtb --bias-max 100
pyar react A.xyz B.xyz --backend xtb --bias-max 100 --check
pyar react A.xyz B.xyz --backend xtb --bias-max 100 --orientations 32
pyar react A.xyz B.xyz --backend xtb --bias-max 150 --bias-alpha-margin 0.001
pyar react A.xyz B.xyz --backend xtb --bias-max 100 --bias-potential softmin --softmin-beta 1.0
pyar react A.xyz B.xyz --backend xtb --bias-controller fixed --bias-min 100 --bias-max 1000
```

Exactly two XYZ reactants, an explicit backend, and `--bias-max` are required.
The bias ceiling is a chemical-problem-specific choice in kJ/mol; the example
values above are syntax examples, **not universal recommended bias strengths**.
Modern defaults are AFIR, accepted-step adaptive control, geomeTRIC/TRIC biased
optimization, eight deterministic orientations, and the established proximity
factor 2.3. Fixed and scheduled controllers retain the canonical bias-range
semantics and require `--bias-min`. Use `--bias-scheduled-alpha` with an explicit
scheduled controller. `--site I J` uses 0-based indices local to A and B.

Charge defaults to zero per fragment. Multiplicity is inferred from electron
parity when omitted. One `--charge`/`--multiplicity` value applies to both inputs;
two values apply separately. Both fragment and merged states are validated.
Standalone xTB uses the existing GFN2 Cartesian provider, without DFT method or
basis settings or a Turbomole requirement. ORCA requires an explicit method and,
for DFT, basis; its built-in xTB methods must omit basis. External ORCA g-xTB and
composite methods are rejected for this reaction gradient route. Gaussian is
currently rejected because its legacy native adapter does not request geometry
optimization, so a valid native release cannot be guaranteed. Legacy behavior
is unchanged. AIMNet2 uses the existing portable package-relative model/runtime
asset validation; missing assets fail preflight, without consulting historical
machine-specific global defaults.

Every run validates inputs, controller/settings, both biased and native-release
capabilities/dependencies, OpenBabel product identity, and restart compatibility
before execution. `--check` runs that same read-only validation and prints the
resolved request without creating directories, generating orientations, running
programs, or changing state. Pip-installable geomeTRIC/ASE support is supplied by
`pyar-chem[xtb]`; xTB, ORCA and the OpenBabel executable are external requirements.

Biased candidates still undergo the existing native unbiased relaxation and
canonical chemical-identity gate before product acceptance. `reaction/state.json`
retains the existing restart protocol; compatible running states resume and
incompatible/completed states are rejected. Existing legacy checkpoint policies
are preserved; checkpoints lacking physical-energy provenance cannot resume on
this route. No automatic NEB/TS/frequency/IRC pipeline is added.

Inspect an orientation's existing reaction trace with
`pyar trace reaction/gamma_adaptive/orientation_000_geom/job_adaptive_000_geom`.
The completion report includes workflow status, product count/paths, run directory
and state path. A completed zero-product search succeeds. Legacy `pyar-react`,
`pyar-cli react`, and `pyar-cli -r` retain their existing parser/default semantics.
