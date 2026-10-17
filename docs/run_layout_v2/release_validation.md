# Run-layout v2 release validation

PR 8 closes the run-layout v2 implementation series; it does not change
PyAR's version, publish a package, retire v1 continuation or distribute new model
weights. Public behavior and compatibility are described in `../run_layout.rst`.
Live CI and final merge status are authoritative on the linked GitHub PR.

## Local evidence (2026-10-15)

- Full suite: **1,345 passed, 1 skipped, 194 subtests passed**. The one skip is
  the in-tree wheel-content test when distributions are built outside `dist/`;
  both final archives were separately inspected by `tools/check_distributions.py`.
- Focused CLI/reaction/path review: 211 passed. Scan reuse/backend review:
  128 passed. The full suite retains the integrated discovery-to-endpoint gate
  tests, negative frequency gate, seeded real RDKit comparisons, legacy
  continuation fixtures and interrupted publication recovery.
- Strict Sphinx with warnings as errors passes. Source and wheel builds pass.
  Distribution inspection checks all v2 owner modules and excludes caches,
  development files and unprovenanced large model assets.
- Clean wheel and sdist core installations were exercised outside the checkout.
  Core-only imports include the new layout/state modules without forcing ASE,
  RDKit, Torch, geomeTRIC, Sella or descriptor extras. CLI help, doctor,
  installation inventory and missing-extra errors pass; dependency checks pass.
- Installed-wheel aggregate/grow geometry-only workflows, microsolvation, and
  real conformer generation pass outside the checkout. Manifest artifact checks
  pass. RDKit and aggregate dependencies were added only to the dedicated smoke
  environment, separate from the clean core environments.
- All **787 original v1 artifact hashes across 44 snapshots** still match the
  PR 0 baseline inventory. Original snapshot bytes were read, not rewritten.
- Public `inspect` reads a real v2 scan run's existing run metadata successfully.

The local development interpreter was Python 3.14. CI additionally covers Python
3.10--3.13, clean core/optional installations, Sella and geodesic stages. PR 8
adds an explicit lightweight conformer-installation matrix and installed v2
geometry/conformer smoke execution; CI results must pass before merge.

## Live backend scans

Water plus helium, one orientation, two scan points ending at 2.0 Angstrom, with
single-thread execution:

```bash
python -m pyar.modern_cli scan-bond water.xyz helium.xyz --atoms 0 0 \
  --backend xtb --xtb-model gfn2 --orientations 1 \
  --scan-end 2.0 --scan-points 2 --output xtb_scan
python -m pyar.modern_cli scan-bond water.xyz helium.xyz --atoms 0 0 \
  --backend orca --method BP86 --basis def2-SVP --orientations 1 \
  --scan-end 2.0 --scan-points 2 --output orca_scan
```

Both complete using real executables (GFN2-xTB and ORCA 6.1.1). Repeating the
same requests preserves request/scan-commit and artifact bytes **and mtimes**:
13 checked xTB files and 27 ORCA files. This uncovered and corrected an ORCA
loader side effect: committed v2 scans now use read-only normalization reuse;
v1 retains its original publication behavior. The regression test now checks
both hashes and mtimes. No backend method, grid or scientific gate was changed.

## Limits and observed backend issue

- The live smoke establishes constrained scan execution and reuse, not a TS or
  validated reaction path. The full path-chain regression executes orchestration,
  trace analysis, product decisions and scientific gates with numerical and
  identity-tool boundaries replaced; it is not a live external reaction search.
- An initial HF/STO-3G ORCA probe failed because the existing keyword builder
  adds `RI` without specifying HF-exchange handling. ORCA reported that RI was
  enabled but exchange needed RIJDX/RIJCOSX/RIJK treatment. The method builder is
  unchanged from the pre-series baseline; this layout PR does not change that
  scientific input policy. The documented BP86/def2-SVP route passed. HF input
  handling requires a separately scoped backend fix.
- Gaussian, Psi4, MOPAC and Turbomole executables were unavailable locally.
  Torch/AIMNet2 model inference and MLatom/TorchANI were not locally exercised;
  CI checks their optional installation/runtime imports where supported. Model
  weights remain external. Their absence is not evidence of backend correctness.
- This series validates native relative artifact refs, not arbitrary relocation
  of completed calculator trees: native cache files can retain absolute paths.
- Existing scientific behavior remains the compatibility target. In particular,
  legacy aggregate cross-path output can contain intermediate composition groups;
  the manifest reports actual composition rather than silently filtering them.
