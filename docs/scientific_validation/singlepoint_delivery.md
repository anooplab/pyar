# Standalone calculation delivery validation

Engineering validation performed on 2026-10-25 against base `529ec6b` plus the
batch energy/engrad/frequency implementation. This is implementation evidence,
not a release qualification or a scientific accuracy benchmark.

## Executed checks

- Full suite: **1845 passed, 10 skipped, 194 subtests passed**. The 10 skips
  are explicitly opt-in real single-point tests; the optional-module inventory
  and `check_test_outcomes.py` accept and report them. No expected failures.
- Full branch-coverage run: **76.03%**, exceeding the existing 50% suite gate.
- Focused deterministic suite: **90 passed, 10 skipped**.
- Enabled real tests: **10 passed** across the minimum/gradient/frequency and
  solvent/open-shell test selections. xTB 6.7.1 and ORCA 6.1.1 were executed;
  CPU PyTorch, TorchANI 2.0 and the installed AIMNet2 model were executed.
- Strict Sphinx, documentation contract (227 CLI examples), Ruff and
  `git diff --check` passed.
- Wheel and sdist builds passed. Both were installed non-editably in separate
  temporary environments outside the checkout, with core dependencies only.
  Real xTB energy, engrad and numerical frequency ran in both installations.
  The deliberately unoptimized H2 frequency completed but correctly returned
  a nonstationary scientific gate failure.
- Distribution inspection found all seven single-point package modules and
  no coverage or complexity-tool state in either artifact.

## Complexity and coverage

Radon 6.0.1 measured all **61 new functions**, including the evaluation closure, in `pyar/singlepoint/` and
`pyar/scripts/singlepoint.py`. Maximum cyclomatic complexity is **6**.
Function CRAP uses `CC^2 * (1 - executed_statement_fraction)^3 + CC`;
maximum CRAP is **6.0**. All functions meet the local target of 6.
This measure uses executed statements, not a claim of full branch coverage.

The focused coverage run measured **534/537 statements (99.44%)** and
**168/172 branches (97.67%)** in the new modules. Coverage runs use an external
data file to keep operational data out of the source tree. The existing
coverage tool's module-specific source probing caused NumPy import errors in
Python 3.14; using `--source=pyar` avoids that probe and measures successfully.

## Real scientific contract checks

- H2 minima for xTB, ORCA HF/STO-3G and AIMNet2: independent numerical location,
  gradient versus central energy differences, projected Hessian and positive
  vibrational mode. ANI uses optimized water because H2 was not a suitable
  minimum smoke case for ANI1x in the selected interval.
- Water xTB ALPB and GBSA; ORCA CPCM and SMD: finite energies/gradients,
  solvent/gas differences, and consistency between energy-only and gradient
  evaluations at identical geometry and settings.
- Charged doublet H2+ gradients with xTB and ORCA.
- Deterministic analytic diatomic Hessian/frequency checks, linear/nonlinear
  rigid-mode counts, nonlinear rotation/translation invariance, imaginary-mode
  classification and unit/sign conversions.

## Limits

Gaussian SCF HF/DFT and ASE's conventional Turbomole dscf/grad route have
manual-based input and mocked adapter coverage, not local executable tests.
Their presence does not assert installation-specific validation. Other native
optimizer adapters do not automatically acquire standalone capabilities.

Frequencies use central differences of gradients, not analytical Hessians.
No thermochemistry or automatic composite/solvent corrections are assembled.
Minimum qualification is limited to the recorded gradient and imaginary-mode
criteria; reaction TS qualification and broad method/element benchmarks remain
outside this delivery. All small negative modes remain visible for review.
