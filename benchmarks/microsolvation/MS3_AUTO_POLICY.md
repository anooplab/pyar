# MS3 case-dependent automatic wall policy probe

This probe verifies that ``--wall auto`` derives a geometry-specific envelope,
leaves a candidate alone when its unbiased optimization remains shell-valid,
and invokes a wall-guided rescue after an unbiased candidate escapes. It is a
single deliberately difficult case, not a benchmark of search quality or a
scientific validation of universal wall parameters.

## Policy

The trigger is the existing post-optimization shell-validity check. The wall
activates only after an unbiased optimization returns usable coordinates and
fails that check. Its radius is

``max accessible sampled solute surface radius from solute COM``
``+ max solvent atom radius from solvent COM + 0.5 Å``.

Its harmonic constant makes a 0.5 Å overrun cost 3 kBT at a 298.15 K
reference, or approximately 0.02266044 Eh/Å². The reference sets a soft search
stiffness scale; it is not an implicit simulation temperature or a solvent
free-energy model. A wall search is followed by an unbiased relaxation, and
only the latter physical energy is used for candidate ranking. Cluster
objective runs disable this policy because a finite aggregate is not a shell
confinement target.

## Live probe

- Solute: methane (5 atoms); solvent: water (3 atoms); one solvent fragment.
- xTB 6.7.1 with geomeTRIC, 150 optimization cycles, shell tolerance 3.5 Å.
- The water was deliberately initialized 12 Å from methane to force an
  unbiased escape; this is a rescue-path smoke test, not a representative
  placement distribution.
- Repository-provisioned Python 3.14.7. Raw driver and output are outside the
  checkout at `/home/anoop/tmp/ms3-pilot/auto_policy_escape.py` and
  `/home/anoop/tmp/ms3-pilot/auto-policy.log`.

The policy derived a 4.984375 Å radius and 0.02266044 Eh/Å² force constant.
Unbiased optimization failed shell validity, so the wall was applied; the
following unbiased relaxation produced a shell-valid structure with the
maximum solvent-center-to-surface distance 0.7522 Å. The wall-guided endpoint
reported 0.5558 Eh maximum wall energy during the trajectory and zero at the
endpoint. The final physical energy was -9.24728878 Eh; the biased objective
was recorded separately. This verifies the state transition and energy
separation for one xTB case only. A separate xTB run with the deterministic
near-surface placement remained shell-valid after its unbiased optimization;
it completed with one optimization, made no wall-search call, and returned
physical energy -9.24722880 Eh. Automated tests also confirm that the radius
changes with solvent fragment extent.

Logs and generated structures are not tracked. No conclusion about improved
minimum recovery, solvent thermodynamics, or a universally optimal radius or
force constant follows from this probe.
