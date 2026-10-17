# MS0 pilot report

Date: 2026-10-17. This is a finite-cluster comparison pilot, not a bulk
solvation study and not a change to PyAR's production ranking. All primary
systems contain four waters: methane (17 atoms), ammonia (16), formamide (18),
and acetate (19). Energies are common vacuum GFN2-xTB single-point values in
hartree. Compare values only within one solute and water count.

## Protocol and software

The four project-generated fixtures are in `fixtures/`; they are neutral
singlets except acetate, which is an isolated -1 singlet without a counterion.
Solute atom order is the prefix and each water is ordered O,H,H. PyAR used 8
orientations per growth step, up to 4 seeds, and 32 Fibonacci surface points
per solute atom. PyAR was run in both geometry-only and xTB-optimized modes.
ORCA SOLVATOR used GFN2-xTB with `NSOLV 4`, `FIXSOLUTE FALSE`, and
`VACUUMSEARCH TRUE`. CREST QCG grow used GFN2-xTB, `-nofix`, and `-samerand`;
the separate conformer run used the completed methane QCG structure. Every
external process used one thread. The resource budget was one repeat per
method/case and a 3600-second per-job timeout.

Common rescoring used GFN2-xTB in vacuum with the case charge and `uhf=0`.
The analyzer validated atom counts, solute-prefix identity, and the O,H,H
fragment order. It rescored all saved PyAR step selections (17 structures per
case/mode) and the available final structures. PyAR coverage and shell-filter
diagnostics are reported separately; hydrogen-bond counts and aligned RMSD
diversity were not evaluated in this pilot.

Installed versions: PyAR 1.5.1, Python 3.14.7, CREST 2.12 (2022 build), xTB
6.7.1, xtb-IFF 1.1, and ORCA 6.1.1. ORCA used `XTBEXE` set to the installed
xTB executable. CREST's banner lists xTB 6.4.0 as its compatible version; the
successful QCG grow runs on xTB 6.7.1 are host-specific evidence, not a general
compatibility claim. Raw structures and logs remain outside the repository.

## n=4 results

Each PyAR value is the lowest energy among its four final selected structures;
the spread is the maximum-minus-minimum within that pool. ORCA and CREST values
are single-structure observations and have no pool spread.

| Solute + 4 waters | PyAR geometry-only: candidates; minimum; spread (kcal/mol) | PyAR xTB: candidates; minimum; spread (kcal/mol) | ORCA SOLVATOR: candidates; minimum | CREST QCG grow: candidates; minimum |
|---|---|---|---|---|
| CH4(H2O)4, 17 atoms | 4; -24.465593115450; 4.355 | 4; -24.501062487567; 5.919 | 1; -24.493933023481 | 1; -24.503624429488 |
| NH3(H2O)4, 16 atoms | 4; -24.714386968916; 6.274 | 4; -24.759703560400; 9.533 | 1; -24.758598750859 | 1; -24.758325078245 |
| Formamide(H2O)4, 18 atoms | 4; -30.936198625890; 21.192 | 4; -30.982996213817; 7.680 | 1; -30.982676687535 | failed during CREST initialization |
| Acetate-(H2O)4, 19 atoms, charge -1 | 4; -34.413433349272; 28.556 | 4; -34.541062555129; 5.373 | 1; -34.540304060982 | 1; -34.533598482986 |

The methane post-QCG CREST conformer run completed in 107 seconds. CREST
pre-optimized 686 structures and retained one unique conformer; its common
rescore was -24.504074853282 Eh. The external pools therefore remain
singletons in these tests. No ranking or uncertainty claim follows from the
small number of independent candidates.

At the last PyAR growth step, each case generated 32 placements and selected
four shell-valid structures. Shell filtering rejected one methane placement,
two ammonia placements, two formamide placements, and one acetate placement.
Mean selected solute-surface coverage was 0.3195, 0.3075, 0.2247, and 0.1935,
respectively. PyAR geometry-only structures are unrelaxed, which is reflected
in their wider energy spreads, especially for formamide and acetate.

## External-tool failures and recovery

The initial CREST QCG attempt could not find xtb-IFF. The official xtb-IFF 1.1
release was installed at user scope and its version check passed; the archive
SHA-256 was `6d5dee10ba39ecb1e26374b7bf82e2364d8ebed5b6f35cdce4ced73245c296f5`.
Methane, ammonia, and acetate QCG grow runs then completed and were rescored.

CREST 2.12's QCG ensemble path failed after growth with a Fortran
sequential-read-after-EOF runtime error. The first attempt used the newer
`-mdtime` spelling; retrying with the installed version's `-mdlen 0.2` spelling
reproduced the same parser failure. Formamide QCG grow separately aborted with
`free(): invalid size` during initialization. Both errors and full logs are
preserved under the scratch benchmark run directories; the unaffected ORCA
and PyAR cases completed. The methane post-QCG CREST conformer run provided a
successful fallback ensemble workflow, but it ended with one unique candidate.

## Interpretation and next gate

The runs establish that solvated final structures of comparable size can be
generated for all four planned cases and rescored under one Hamiltonian. The
acetate case explicitly represents an isolated finite anion, without a
counterion. Search potentials differ across PyAR, CREST, and ORCA; common
rescoring does not erase those constraints. The external candidate pools are
too small for method ranking, RMSD-diversity claims, or statistical uncertainty
estimates, and no bulk thermodynamics are inferred. MS0 changes no production
defaults or ranking. MS1 requires separate scientific approval.
