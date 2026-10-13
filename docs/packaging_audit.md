# Packaging audit and implementation record

Audited `master` before edits on 2026-10-13; implementation branch:
`packaging-workflow-extras`. Scientific defaults, algorithms, console entry
points, Python >=3.10 and setuptools remain the compatibility contract.

## Baseline state (before implementation)

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
| Autograd | no production imports remain; tests retain an independent derivative reference | test dependency group |
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

## Implemented extras and compatibility

Public workflows: aggregate, grow (alias), conformer, reaction, microsolvate,
solvate (compatibility alias), neb. Python capabilities: geometric, aimnet2,
sella, geodesic, structure-comparison, identity, analysis, ml.

Keep useful selection and xtb aliases with migration notes (xtb installs Python
optimizer support only). Remove misleading/unreliable openbabel pip extra;
document conda-forge/system installation. Move test/docs/release/dev tooling
to standard dependency groups, update CI and ReadTheDocs first. Retain `all`
as self-referencing end-user extras only.

The runtime registry owns module-to-extra hints and executable installation
sources/probes, shared by preflight and doctor. Metadata consistency tests
check declared extras against this registry. Optional components do not make
doctor fail. Missing requirements must fail before calculation or output files.

## Authoritative dependency references

- AIMNet package: https://pypi.org/project/aimnet/ (current 0.2.0 requires 3.11+).
- AIMNet source/API: https://github.com/isayevlab/aimnetcentral/
- Dependency groups: https://packaging.python.org/en/latest/specifications/dependency-groups/
- pip group support (25.1+): https://pip.pypa.io/en/stable/user_guide/#dependency-groups
- ORCA installation: https://www.faccts.de/docs/orca/6.1/tutorials/first_steps/install.html

## Final installation contract

Core is NumPy, SciPy, NetworkX, and tomli on Python 3.10. Autograd was removed
because no runtime imports remain. Python >=3.10, setuptools, version 1.5.1,
console entry points, scientific defaults, and the subprocess xTB architecture
are preserved.

| Workflow extra | Installed Python stack, excluding core |
| --- | --- |
| aggregate | ASE, DScribe, HDBSCAN, scikit-learn |
| grow | aggregate alias |
| conformer | RDKit |
| reaction | aggregate + geometric |
| microsolvate | aggregate alias |
| solvate | microsolvate compatibility alias |
| neb | geometric |

| Capability extra | Installed Python stack |
| --- | --- |
| geometric | ASE, geomeTRIC >=1.1 |
| aimnet2 | ASE, Torch; retained legacy calculator with external model |
| sella | Sella >=2.6.0 |
| geodesic | geodesic-interpolate >=1.0.0 |
| structure-comparison | irmsd >=0.1.0 |
| identity | RDKit |
| analysis | ASE, MDAnalysis, matplotlib, pandas |
| ml | MLatom, Torch, TorchANI 2.0, ASE, setuptools >=69,<81, requests, joblib |

Clean ML imports exposed two upstream metadata omissions: MLatom imports
requests and joblib eagerly without declaring them. TorchANI 2.0 imports
pkg_resources, removed in newer setuptools. The ML extra supplies these
requirements explicitly; h5py and pyh5md remain MLatom transitive requirements.
No dependency classification remains unresolved for the audited current routes.

### Compatibility mapping

| Previous installation | Current installation |
| --- | --- |
| core | core, with unused Autograd removed |
| conformer | conformer (preserved) |
| selection | aggregate + analysis (compatibility alias retained) |
| xtb | geometric (alias retained); install the xTB executable separately |
| openbabel | removed unreliable source-build extra; install system/conda tools |
| test, docs | standard test/docs dependency groups |
| all | composable end-user capabilities; developer tools excluded |
| setup.py / requirements.txt | pyproject.toml; obsolete unconsumed duplicates removed |
| docs/requirements.txt | docs dependency group; CI and ReadTheDocs consumers updated |

`pyar.dependencies` owns workflow Python checks, executable sources, supported
probes, and installation hints. The backend registry owns selected route needs;
metadata consistency tests tie both registries to declared extras. Normal
execution and existing check-only commands fail early on missing dependencies.
Doctor reports core, optional Python modules, external programs, and versions;
missing optional components do not make a core installation unhealthy.

External discovery covers xTB, ORCA, Gaussian g16, Psi4, MOPAC, OpenBabel
(obabel/obminimize/obenergy), Turbomole (define/ridft/rdgrad/statpt/jobex), and
geometric-optimize. Installation sources and safe verification guidance are in
[external_programs.rst](external_programs.rst). Vendor programs without a safe
version flag use PATH discovery and a vendor sample calculation; PyAR does not
invoke those programs without inputs as a version probe.

### Canonical examples

```bash
python -m pip install "pyar-chem[aggregate]"
# Install xTB separately, then:
pyar aggregate molecule.xyz --size 2 --backend xtb

python -m pip install "pyar-chem[conformer,aimnet2]"
# Configure PYAR_AIMNET2_MODEL as documented in installation.rst, then:
pyar conformer CCO --input-format smiles --backend aimnet_2

python -m pip install "pyar-chem[reaction]"
# Install ORCA separately, then check the selected route:
pyar react first.xyz second.xyz --backend orca --method BP86 --basis def2-SVP --check
```

The model remains external to distributions; the existing checkout model and
adapter are retained. Native XYZ optimization with an explicit charge uses ASE
and does not require OpenBabel Python bindings. Legacy format/charge inference
still reports its OpenBabel requirement. Official aimnet migration was deferred
because its current Python >=3.11 requirement and different calculator/model
contract do not preserve PyAR's Python 3.10 and legacy model behavior.

## Validation record

- Full suite: 1268 passed (including the existing subtest coverage); 685 warnings. Coverage run:
  71.53%, above the existing 50% gate.
- Clean non-editable core, aggregate, reaction, conformer+aimnet2, and
  neb+sella+geodesic+structure-comparison installs passed import/CLI checks.
  Grow, microsolvate, and solvate passed against the clean aggregation stack,
  which is exactly their declared dependency closure.
- Python 3.10: core, aggregation, conformer+AIMNet2 and pip check passed.
  Python 3.13: core, reaction, conformer+AIMNet2 and pip check passed.
  Python 3.11 supplied the other representative local checks. CI covers
  eight clean installation modes on both 3.10 and 3.13; the original run passed all clean-install jobs; full-suite collection exposed the
  missing test-only Autograd dependency, corrected in the review fixes.
- Actual xTB aggregation completed, selected a finite four-atom geometry,
  and its raw backend log reported normal termination and a finite energy.
- Actual AIMNet2 methane conformer refinement completed with a finite energy
  and selected five-atom geometry using the existing model, with no OpenBabel
  Python bindings installed.
- Missing-xTB and missing-ORCA checks exit with actionable official guidance,
  no traceback, and no calculation outputs. SOAP/HDBSCAN are exercised directly
  in the clean aggregate install, without Torch or geomeTRIC.
- Final wheel and sdist: archive inspection rejects caches, model weights,
  obsolete binaries/machinery, and wheel development files; 212 wheel entries
  and 363 sdist entries passed. Twine metadata checks passed. Core installs
  from both archives passed help, doctor and pip check. Representative workflow
  wheel checks exercise aggregation, reaction, grow, solvate and microsolvation.
- Strict Sphinx builds passed in both the existing development environment and
  a clean docs-group environment, with warnings treated as errors. The standard
  pip group command and pip check passed. CI and ReadTheDocs YAML parsed cleanly.
- The complete all extra resolved on Python 3.13 with CPU Torch. Pip check,
  external MLatom/TorchANI imports, analysis entry requirements, every workflow
  requirement, and Sella/geodesic/iRMSD/AIMNet2 calculator imports passed.
  Upstream syntax warnings remain; these are not missing dependencies.

## Files and delivery

Changes are scoped to pyproject metadata, capability/preflight and optional
import handling, the AIMNet2 adapter's explicit-charge installation path,
installation/migration/backend docs, ReadTheDocs/CI, distribution/smoke tooling,
and focused packaging/CLI tests. Removed setup.py, build_backend.py and duplicate
requirements files only after consumer inspection. Generated API files and
build artifacts are excluded from the change. The branch is
`packaging-workflow-extras`; no release, version bump, or remote push is included.

Changed paths (relative to the repository):

- `.github/workflows/ci.yml`
- `.readthedocs.yaml`
- `README.md`
- `build_backend.py`
- `docs/backends.rst`
- `docs/external_programs.rst`
- `docs/installation.rst`
- `docs/migration.rst`
- `docs/packaging_audit.md`
- `docs/requirements.txt`
- `pyar/AIMNet2/calculators/aimnet2_ase_opt.py`
- `pyar/backend_capabilities.py`
- `pyar/backends/__init__.py`
- `pyar/backends/aimnet2_assets.py`
- `pyar/backends/aimnet_2.py`
- `pyar/cli.py`
- `pyar/dependencies.py`
- `pyar/energy_gradient_providers.py`
- `pyar/neb.py`
- `pyar/optimization_request.py`
- `pyar/optional_dependencies.py`
- `pyar/reaction_request.py`
- `pyar/scripts/cli_meta.py`
- `pyar/scripts/conformer.py`
- `pyar/scripts/descriptor.py`
- `pyar/scripts/grow.py`
- `pyar/scripts/modern_aggregate.py`
- `pyar/scripts/modern_conformer.py`
- `pyar/scripts/modern_grow.py`
- `pyar/scripts/modern_microsolvate.py`
- `pyar/scripts/modern_path_stage.py`
- `pyar/scripts/modern_react.py`
- `pyar/scripts/react.py`
- `pyproject.toml`
- `requirements.txt`
- `setup.py`
- `tests/test_cli.py`
- `tests/test_installation_contract.py`
- `tests/test_modern_react.py`
- `tests/test_packaging_metadata.py`
- `tools/check_distributions.py`
- `tools/installation_smoke.py`

## Follow-up issues

1. Evaluate official aimnet migration with an explicit Python-support decision,
   model provenance/license review, and calculator/energy equivalence tests.
2. Evaluate the Python xTB API separately; preserve subprocess model/capability
   checks until behavior, distribution, and failure-state equivalence are proven.
3. Upgrade TorchANI and remove the pkg_resources/setuptools compatibility pin
   only after validating ANI model/API behavior. Report MLatom's missing requests
   and joblib metadata upstream.
4. Run vendor-licensed ORCA/Gaussian/Turbomole integrations on an authorized host;
   this work verifies their discovery/error paths, not vendor calculations.
5. Run the configured remote CI and platform-specific installation matrix before
   publishing a release. Local validation is Linux; CPU Torch was used to avoid
   installing unrelated GPU runtime packages.

## Requirement completion audit

Sections 1–5: baseline import audit, minimal core, complete workflow closures,
composable capabilities, and actual clean installations. Section 6: official
AIMNet compatibility evaluated; current adapter retained with tested Python stack.
Sections 7–11: external executable registry, unchanged xTB architecture, shared
metadata, automatic/check-only preflight, and healthy optional-aware doctor.
Sections 12–15: standard developer groups, runtime-only all, consumer-safe cleanup,
and retained setuptools/layout. Sections 16–20: real-install CI and archives,
unchanged Python support, task docs and migration aliases, focused acceptance
tests. Sections 21–22: canonical stories verified through positive xTB/AIMNet and
negative ORCA paths; scientific scope preserved. Sections 23–24: isolated branch,
coherent commits, concrete verification and this detailed implementation report.
All local implementation and validation gates above passed. Vendor-licensed calculations remain follow-ups. The review fixes add test-only
Autograd, exercise Python 3.10 conditional core requirements, normalize relative
AIMNet2 model paths before directory changes, and preserve every explicit-charge
XYZ frame. Updated full-suite and distribution checks passed; the corrected
remote full-suite run is checked separately.
