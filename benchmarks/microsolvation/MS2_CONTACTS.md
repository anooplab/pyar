# MS2 contact-reporting pilot

**Purpose:** check that the opt-in post-relax atom-contact report is attached to
the existing selection artifacts and does not alter legacy `--site` sampling,
energy ranking, or selected structures.

**System and settings:** the repository acetate anion fixture plus one water,
with fragment charges `-1, 0` and singlet multiplicities `1, 1`; xTB 6.7.1;
one solvent, four orientations, two survivor slots, and the same initial
`--site 5 6` in both arms. The reporting run additionally specified
`--site-group carboxylate=5,6 --contact-distance 3.4`. Both runs completed.

| Rank | xTB energy (Eh) | Legacy structure SHA-256 prefix | Site-report run SHA-256 prefix |
|---:|---:|---|---|
| 1 | -19.260996660132 | `357b168bff6f` | `357b168bff6f` |
| 2 | -19.263795987094 | `cdd7625668f4` | `cdd7625668f4` |

Both selected structures had one water in geometric contact with the named
carboxylate group at the 3.4 Å cutoff. The closest reported site contacts
were 1.624 and 1.759 Å. The covalent-radii proximity check found the original
solute connectivity unchanged in both selected structures. These distances
are geometric evidence only; this small pilot does not establish a chemical
bond, a general contact threshold, or improved sampling. Selected file hashes
and energies match exactly between runs, so additive reporting left this
legacy targeting case unchanged. One water does not assess second-shell
population or multi-site allocation.

Reproduce from the repository root (outputs remain outside the tree):

```sh
pyar microsolvate benchmarks/microsolvation/fixtures/acetate.xyz \
  benchmarks/microsolvation/fixtures/water.xyz --count 1 --site 5 6 \
  --orientations 4 --maximum-number-of-seeds 2 --backend xtb \
  --charge -1 0 --multiplicity 1 1 --output /home/anoop/tmp/ms2-legacy

pyar microsolvate benchmarks/microsolvation/fixtures/acetate.xyz \
  benchmarks/microsolvation/fixtures/water.xyz --count 1 --site 5 6 \
  --orientations 4 --maximum-number-of-seeds 2 --backend xtb \
  --charge -1 0 --multiplicity 1 1 \
  --site-group carboxylate=5,6 --contact-distance 3.4 \
  --output /home/anoop/tmp/ms2-contacts
```
