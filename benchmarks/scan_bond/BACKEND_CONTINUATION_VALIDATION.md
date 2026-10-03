# Backend scans and reaction-path continuation validation

Validated on 2026-10-03 using geomeTRIC 1.1.1, ASE 3.29.0,
xTB 6.7.1 and ORCA 6.1.1. These are implementation checks, not optimizer
performance benchmarks or evidence of a chemical reaction mechanism.

## Interface decisions

`--software` selects the energy/gradient backend. ORCA remains the default
and uses native relaxed scans; the other registered providers use constrained
geomeTRIC internal-coordinate minimization at each scan point. Constrained
coordinates use geomeTRIC's orthogonalized constraint treatment (`conmethod=1`).
A three-atom analytic pair potential verifies that both the prescribed bond
distance and the freely optimized other distances are correct; nonzero forces
along the constrained distance are allowed.

`--through` selects the last cumulative stage. Its values are `scan` (the CLI
default), `neb`, `ts`, `frequency`, `irc`, `endpoints`, `endpoint-frequency`,
and `all`. `endpoints` stops after both IRC endpoint optimizations;
`endpoint-frequency` and `all` additionally verify minima and reaction
connection. The Python API's omitted `through` retains its previous
scan-plus-relaxation behavior.

The scan waypoint is only an NEB initializer. The TS stage consumes the
standardized NEB candidate artifact, preferring the final climbing image.
All frequency, IRC and endpoint confirmation gates remain independent.

## Automated coverage

- Registered xTB, ORCA, Gaussian and AIMNet2 routing with unchanged backend,
  merged charge/multiplicity, process count and requested optimizer.
- Every cumulative stopping point and every scientific gate failure.
- Complete scan-to-validation workflow with expensive numerical boundaries
  mocked, plus extension from a completed NEB without rescanning or rerunning NEB.
- Real constrained geomeTRIC optimization against an analytic potential.
- Incomplete/modified scan artifacts rejected during recovery.
- Nonconverged constrained points do not produce a completed scan.
- Endpoint frequencies consume optimized endpoint artifacts without repeating
  optimization or mutating the upstream artifacts.
- Failed endpoint optimization blocks subsequent frequency validation.
- Stage reuse respects optimizer parameters, physical settings and artifact
  and upstream summary hashes.
- ORCA built-in xTB gradient keywords omit DFT-only settings; unrestricted
  singlet ORCA DFT remains unrestricted in continuation.
- CLI scientific failures have a nonzero exit status.

Run from the repository environment:

```bash
python -m pytest -q tests/test_bond_scan_backends.py tests/test_scan_bond.py \
  tests/test_neb.py tests/test_energy_gradient_providers.py
python -m pytest -q
python -m sphinx -E -b html -W --keep-going docs /tmp/pyar-scan-docs
python -m build
pyar-neb --help
pyar-cli scan-bond --help
git diff --check
```

## Real executable smoke checks

A three-point H2 compression scan from 0.96 to 0.74 angstrom completed using
standalone GFN2-xTB and native ORCA GFN2-xTB. All point energies and geometries
were finite and the retained trajectories passed the scan-grid checks.

| Scan backend | Points | Final energy / Hartree | Result |
| --- | ---: | ---: | --- |
| Standalone xTB, `--gfn 2` | 3 | -0.981983694723 | Complete |
| ORCA, `GFN2-xTB` | 3 | -0.98198369 | Complete |

A separate ORCA GFN2-xTB provider call returned energy -0.98198369872 Hartree
and gradient norm 0.029902748285794926 Hartree/Bohr, both finite. ORCA required
its documented `XTBEXE` environment configuration.

The installed xTB executable does not advertise `--gxtb`. Its g-xTB scan
was explicitly rejected. It did not silently run GFN2-xTB. Real g-xTB smoke
validation therefore remains dependent on installing an executable that
supports that Hamiltonian. Gaussian/AIMNet2 continuation is tested through
provider/workflow boundaries, without real Gaussian/AIMNet2 smoke calculations.
A complete real scan-to-IRC reaction calculation was not performed here.

Repeat the real smoke checks with the saved script:

```bash
# First configure XTBEXE for your local ORCA/xTB installation.
python benchmarks/scan_bond/validate_backend_scans.py \
  --output /path/to/new/scan-smoke-results \
  --backend xtb-gfn2 --backend orca-gfn2
```

The script retains inputs, per-backend scan artifacts and `smoke_results.json`
in the supplied directory. Add `--backend xtb-gxtb` only with a capable xTB
executable; an unsupported model is reported as a failed check.

## Review and fixes

Reviewed scientific gates, model dispatch, physical settings, artifact
consistency, restart behavior, CLI failure reporting and backend isolation.
Fixed constraint projection after the analytic test exposed incorrect free
bond relaxation with the upstream default constraint method. Split endpoint
optimization from frequency checks while preserving existing combined
`pyar-neb --stage endpoints` behavior. Preserved upstream optimized geometries
and summary hashes when adding frequencies. Added an explicit gate against
nonconverged endpoint optimizations. Made reuse respect explicitly supplied new starting structures.
Carried SCF type into continuation and
preserved unrestricted singlet ORCA DFT. Prevented one corrupt orientation
from aborting other orientations. Added nonzero CLI failure reporting after
live artifact inspection exposed a misleading zero exit from the old CLI.
