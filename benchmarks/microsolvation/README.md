# Microsolvation MS0 benchmark

MS0 freezes the existing microsolvation implementation and evaluates its
geometry-only and xTB-optimized candidate pools against external structure
generators. It does not change the search algorithm, scoring, defaults, or
restart format. Primary comparisons use four waters: CH4(H2O)4 (17 atoms),
NH3(H2O)4 (16 atoms), formamide(H2O)4 (18 atoms), and acetate(H2O)4 (19 atoms).
Smaller clusters are scaling controls; larger clusters are extensions when
the tools complete within the stated resource budget.

The fixtures are project-generated singlet geometries. Acetate is an isolated
anion with charge -1 and no counterion; keep that limitation explicit in any
chemical interpretation. The fixed solute atom order is the prefix and each added water is O,H,H. All
cross-method energies are vacuum GFN2-xTB single points on the final
coordinates. Compare energies only at the same solute and water count. This
common rescore does not make the search potentials equivalent: PyAR ranks by
coverage then energy/diversity, CREST QCG uses xtb-IFF growth, and ORCA
SOLVATOR uses its docking wall (with `VACUUMSEARCH TRUE` here). Preserve those
constraints in any interpretation. The harness runs CREST QCG grow, a
post-QCG conformer search, and a QCG ensemble attempt with one 0.2 ps MD
segment (`-mdlen`, the installed CREST 2.12 spelling) and a 100 fs dump
interval. Record each mode's status separately; grow and post-QCG conformers
remain usable if the QCG ensemble mode fails.

Method references: [CREST QCG options](https://crest-lab.github.io/crest-docs/page/documentation/qcg.html),
[ORCA 6.1 SOLVATOR](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/solvator.html),
and the [official xtb-IFF release](https://github.com/grimme-lab/xtbiff/releases/tag/v1.1).

## Reproduce

Use the repository's test/development Python environment and keep large raw
outputs outside the checkout. Set `XTBEXE` for ORCA if its xTB executable is
not next to the ORCA installation. The exact commands, paths, versions,
resource limits, failures and result table belong in `REPORT.md` after a run.
`protocol.json` is the machine-readable budget and case matrix.

For the complete n=4 matrix, choose a fresh scratch directory and run:

```bash
python benchmarks/microsolvation/run.py \
  --workdir /path/to/scratch/ms0 --counts 4
python benchmarks/microsolvation/analyze.py /path/to/scratch/ms0
```

The runner records failed methods and returns non-zero if any method fails;
it continues through the remaining cases and preserves each log. The analyzer
rescales saved XYZ candidate pools with each case's charge and multiplicity.

The primary benchmark targets all four solutes at n=4. Do not silently drop a
case if a method fails; record it as `not_run` or `failed` and retain the raw
log. A comparison should use candidate pools, report counts, and separate
geometry-only from common-method-rescored results. Extend to CH4(H2O)8 and
NH3(H2O)6 only if the n=4 runs are robust and the fixed budget allows it. No
bulk thermodynamics or general ranking conclusion follows from this small
benchmark.

## Current status

See `REPORT.md` for the observed pilot evidence and remaining limitations.
