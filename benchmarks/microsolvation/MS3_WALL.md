# MS3 outer-sphere wall pilot

This small pilot validates the optional wall implementation path; it is not a
benchmark of search quality or evidence for a preferred wall setting.

## Protocol

- Solute: methane, five atoms; solvent: one water molecule.
- One added solvent; two deterministic placement orientations; one retained
  seed; shell objective; same xTB GFN2 and geomeTRIC backend for each run.
- xTB 6.7.1, Python 3.14.7. The wall run used a 3.0 Angstrom radius and
  0.5 hartree/Angstrom squared force constant. The control left the wall off.
- Raw runs and logs were kept outside the repository under
  `/home/anoop/tmp/ms3-pilot/` and are not part of this commit.

## Observed result

Both runs completed and retained a shell-valid methane-water structure after
the final unbiased relaxation. The wall run reported search-objective energy
and wall energy separately from the final unbiased physical energy. The wall
contribution at the recorded search endpoint was approximately `3.7e-5` and
`4.5e-5 Eh` for the two candidates. This controlled example does **not** show a better minimum or
improved solvent-contact retention with the wall; it establishes that the
two-stage path executes and that physical ranking excludes the wall term.

The pilot is deliberately too small to support a scientific recommendation.
The contrived methane-water escape probe started the water center 12 Å from
methane. The no-wall xTB optimization converged with the water 8.38 Å from the
sampled methane surface, outside the 3.5 Å shell tolerance. With a 4.5 Å
fragment-center wall (`k = 1.0 Eh/Å²`), the wall-guided stage moved the solvent
into the shell. The subsequent unbiased relaxation retained it 0.54 Å from the
surface; the final physical energy was `-9.246958 Eh`, compared with
`-9.245768 Eh` for the escaped control. The wall contribution at the biased
endpoint was zero because the solvent had reached the flat-bottom region. This
is one deliberately contrived rescue, not evidence of general improvement or
a recommended parameter set. The wall remains experimental and off by default.
