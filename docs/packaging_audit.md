# Packaging audit and implementation record

Audited `master` before edits on 2026-10-13; implementation branch:
`packaging-workflow-extras`. Scientific defaults, algorithms, console entry
points, Python >=3.10 and setuptools remain the compatibility contract.

## Current state

Core declares NumPy, SciPy, NetworkX, Autograd and conditional tomli. Workflow
extras exist only for conformer (RDKit); selection is an internal broad extra.
`all` includes developer tools and fragile OpenBabel source requirements.
The modern CLI already dispatches lazily, implements doctor and check-only
execution, and delegates backend preflight to shared legacy helpers.
`setup.py`, `build_backend.py`, and root `requirements.txt` have no active
consumers; ReadTheDocs still consumes `docs/requirements.txt`.

## Dependency matrix and code evidence

| Dependency | Actual use | Proposed ownership |
| --- | --- | --- |
| NumPy | Molecule, geometry, I/O and all workflows | core |
| SciPy | property distances; RMSD assignment; CLI selection distances | core |
| NetworkX | coordinate/chemical graphs; CLI selection policy | core |
| tomli | package/config TOML fallback on Python 3.10 | conditional core |
| Autograd | no imports remain; bias gradients use NumPy analytical expressions | remove unused requirement |
| ASE | representations; energy-gradient providers; NEB; AIMNet2 | aggregate, reaction, microsolvate, neb, geometric, aimnet2 |
| DScribe | lazy SOAP/MBTR in representations; default growth selection | aggregate, reaction, microsolvate |
| HDBSCAN | first auto clustering attempt in selection/clusterers.py | aggregate, reaction, microsolvate |
| scikit-learn | DBSCAN/OPTICS selection and HDBSCAN dependency | aggregate, reaction, microsolvate |
| RDKit | conformer generation and chemical identity | conformer, identity |
| geomeTRIC | reaction biased relaxation; optional geometry/TS route | reaction, neb, geometric |
| Torch | bundled AIMNet2 calculator; TorchANI and ML optimization | aimnet2, ml |
| aimnet | official current distribution requires Python >=3.11 and a different calculator/model interface | defer migration; preserve supported 3.10 and legacy model semantics |
| Sella | selected NEB TS optimizer, imported lazily | sella |
| geodesic-interpolate | selected NEB path interpolation, imported lazily | geodesic |
| irmsd | optional structure comparison adapter | structure-comparison |
| pandas | legacy clustering CSV and descriptor utility | selection compatibility/analysis |
| MDAnalysis | descriptor utility only | selection compatibility/analysis |
| matplotlib | analysis plots including Turbomole plot helper | analysis |
| MLatom | external compatibility adapter, AIQM1 and mlopt | ml |
| TorchANI | ANI backend; existing version constraint requires compatibility review | ml |
| h5py / pyh5md | no direct imports after MLatom vendor removal | delegate to external ML package metadata after validation |
| OpenBabel | executable identity/conversion, MOPAC, force-field backend; Python bindings used by native AIMNet2 helper | executable discovery; remove native helper's avoidable binding requirement for explicit XYZ/charge |
| xTB | subprocess backend, model capability probe | external executable; no API migration |
| ORCA / Gaussian | native and Cartesian gradient routes | external executables |
| Turbomole | define, ridft, rdgrad, statpt, jobex routes | external executables, route-specific |
| MOPAC | executable plus obabel input conversion | external executables |
| Psi4 | executable adapter, no imported Python API | external executable |
| DFT-D4 | no direct import or standalone adapter found | backend method concern, no invented extra |
| pytest / coverage / build / twine / Sphinx | tests, distribution checks, docs | dependency groups |

Workflow evidence: `workflows/aggregate.py` and `growth/service.py` share
selection; `selection/policy.py` selects SOAP for aggregates and MBTR for
molecules, with HDBSCAN first in `selection/clusterers.py`. Reaction uses that
selection plus Cartesian providers and geomeTRIC. Conformer uses RDKit and
graph comparisons, with separately selected refinement. Microsolvation shares
growth/selection. NEB imports ASE and optimizer stacks at execution time.

## Proposed extras and compatibility

Public workflows: aggregate, grow (alias), conformer, reaction, microsolvate,
solvate (compatibility alias), neb. Python capabilities: geometric, aimnet2,
sella, geodesic, structure-comparison, identity, analysis, ml.

Keep useful selection and xtb aliases with migration notes (xtb installs Python
optimizer support only). Remove misleading/unreliable openbabel pip extra;
document conda-forge/system installation. Move test/docs/release/dev tooling
to standard dependency groups, update CI and ReadTheDocs first. Retain `all`
as self-referencing end-user extras only.

The runtime registry will own module-to-extra hints and executable installation
sources/probes, shared by preflight and doctor. Metadata consistency tests will
check declared extras against this registry. Optional components do not make
doctor fail. Missing requirements must fail before calculation or output files.

## Authoritative dependency references

- AIMNet package: https://pypi.org/project/aimnet/ (current 0.2.0 requires 3.11+).
- AIMNet source/API: https://github.com/isayevlab/aimnetcentral/
- Dependency groups: https://packaging.python.org/en/latest/specifications/dependency-groups/
- pip group support (25.1+): https://pip.pypa.io/en/stable/user_guide/#dependency-groups
- ORCA installation: https://www.faccts.de/docs/orca/6.1/tutorials/first_steps/install.html

## Validation record

Implementation and clean-install validation are pending. This audit records
proposals, not completion claims. Remaining requirements include all selected
Python capability availability checks, route-specific external discovery,
backwards compatibility tests, clean workflow installations, wheel and sdist
inspection/installations, pip check, CLI smoke tests, and documentation build.
