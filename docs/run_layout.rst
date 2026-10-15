Run outputs and compatibility
=============================

Fresh grow, microsolvation, aggregate, scan-bond, reaction and conformer runs use
layout v2. Default run directory names and scientific defaults are unchanged.
Existing recognized v1 runs retain their original owners and paths; running a
new PyAR version does not move or migrate them. Unrecognized nonempty directories
and mixed or corrupt metadata fail closed. Preserve an old run and choose a
fresh directory when a restart is unsupported.

Finding outputs
---------------

Use the returned workflow result and the run's ``results/manifest.json`` rather
than a recursive filename search. Manifests contain run-relative artifact paths,
opaque UUIDs, roles, SHA-256 hashes and source provenance. Numeric filenames are
display ranks, not scientific identity. A role describes the evidence: a
scan-maximum candidate is not a validated TS, and an accepted unbiased reaction
product does not itself establish a validated reaction path.

.. list-table:: Fresh v2 runs
   :header-rows: 1
   :widths: 18 40 42

   * - Run root
     - Scientific records and native work
     - Public results and continuation
   * - ``grow/``
     - ``inputs/``; ``steps/step_NNN/selection.json`` and execution attempts
     - ``results/structure_NNN.xyz``; verified selected pools drive continuation
   * - ``microsolvation/``
     - Original solute/surface, selected-pool and coverage commits per solvent
     - Flat structure copies; original-solute targeting and coverage are retained
   * - ``aggregates/``
     - Explicit pathways, addition selections and native source snapshots
     - Flat aggregate copies, retaining legacy cross-path composition groups;
       a completed aggregation is not resumable
   * - ``scan_bond/``
     - ``orientation_NNN/scan_commit.json``; native scan/candidate/relax attempts
     - Candidate refs in the manifest; ``--through`` can extend a verified scan
       without rescanning
   * - ``reaction/``
     - ``discovery/commit.json`` and accepted products; separate characterization
     - Product copies and source IDs; completed discovery can feed later path
       characterization without rerunning discovery
   * - ``conformers/``
     - ``generated/records.json``; ``refinement/execution``; ``completion.json``
     - Ranked conformer copies; inspection metadata provides no resume API

Scan and reaction path stages share a cumulative native execution directory
within each attempt. Per-stage commits index those files and their existing
scientific gates; the executor's sibling paths and hashes are preserved.
Standalone NEB, TS, frequency and IRC commands retain their previous layouts.
Deprecated legacy ``solvate`` also retains its existing layout and behavior.

For example, verify a manifest artifact without reconstructing its filename:

.. code-block:: python

   from pathlib import Path
   from pyar.persistence import ArtifactRef, read_json

   root = Path("grow").resolve()
   manifest = read_json(root / "results/manifest.json")
   entry = manifest["artifacts"][0]
   reference = ArtifactRef.from_dict({
       key: entry[key] for key in ("artifact_id", "role", "path", "sha256")
   })
   result_path = reference.verify(root)

Use the owning workflow's loader and scientific validation for restart decisions;
a checksum alone does not establish convergence or stage validity. Native files
inside execution trees can retain their historical names and absolute internal
references. The new authoritative artifact refs are run-relative; relocating a
complete native calculation is not a general portability guarantee.

Authority and reports
---------------------

Committed selections, inputs, discovery and completion records are scientific
evidence. Do not edit them or overwrite their artifacts to repair a report.
Progress state and summaries have workflow-specific roles: reaction
characterization depends on the discovery commit, while an in-progress grow
still needs its progress state and verified selections. A final result copy is
not automatically a restart seed. Do not assume that deleting arbitrary state
files is a supported recovery procedure.

For conformers, RDKit force-field energies remain in kcal/mol and physical
backend energies in Hartree, in separate records. Benchmark consumers dispatch
on the explicit marker and verify completion/manifest provenance before reading
selected coordinates. Legacy benchmarks still read the old selected directory.

Reanalysis tools may write new candidate files. To analyze a completed reaction
trace with different settings, work on a copy rather than modifying evidence
referenced by a discovery commit. Trace and candidate locations are recorded in
the product and stage records; use those references instead of assuming a gamma
or orientation directory name.

Writer recovery and compatibility
---------------------------------

V2 writers use a local POSIX advisory lock next to the run root. An active lock
cannot be overridden. Grow, microsolvation, aggregate, scan-bond and modern react
provide ``--recover-stale-lock`` for explicit stale-owner recovery; Python scan
and reaction entry points also accept ``recover_stale=True``. Recovery retains
owner evidence and never turns a failed scientific gate into a successful stage.
``--check`` remains read-only, including when that flag is present. Conformer
runs remain nonresumable: preserve an interrupted run and use a new directory.
Do not delete the lock file to bypass a writer.

There is no automatic v1-to-v2 migration or announced removal of v1 continuation
in this series. CLI aliases, parser defaults, physical backend parameters,
sampling, selection and validation gates are preserved. Fresh-v2 selected paths
change with the layout, so downstream code should use returned paths/manifests;
old runs continue returning their established paths.

Validation and backend scope
----------------------------

The regression suite covers v1/v2 scientific comparisons, interrupted commit
publication, checksum rejection, scan extension, reaction discovery through path
gates and real seeded RDKit generation. Numerical boundary tests are distinct
from live external-calculator runs. CI checks Python 3.10--3.13, clean archives,
core/optional installations, Sella and geodesic paths. Installed aggregate/grow
geometry-only and conformer generation smokes run outside the checkout.

Install Python extras separately from external programs. ``pyar install status
--json`` inventories the current machine; workflow ``--check`` is the relevant
preflight. ORCA/Gaussian/Turbomole and other external backends still require their
own installation and, where applicable, licenses. AIMNet2 model files remain
external and are not newly bundled by this layout series. See :doc:`installation`
and :doc:`persistence` for installation and persistence details.
