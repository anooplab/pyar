Using PyAR 1.6 after an older version
=====================================

PyAR 1.6 requires Python 3.11 or newer; Python 3.10 is no longer supported.
Use a fresh environment with a supported interpreter.

PyAR 1.6 is an intentionally breaking cleanup. The single executable is ``pyar``.
Standalone commands such as ``pyar-react``, ``pyar-cli`` and ``pyar-neb`` and
compatibility-only Python import paths are removed. Use the corresponding
``pyar react``, ``pyar neb`` and other commands in :doc:`cli`.

Old run/config/checkpoint generations are not migrated. Use the environment and
version recorded by the original work to reproduce it. Keep that environment
separate from 1.6. A historical XYZ can be used as a new input, but its comment
line does not replace explicit method, charge/spin, mapping and provenance.
Do not relabel old JSON/TOML with the new schema to bypass validation.

Useful specialist utilities now live under subcommands; see :doc:`tools`.
Current data integrity, restart and scientific acceptance checks remain required.
PyAR 1.6 release/package publication is separate from the repository development
version described by these pages.

Geometry-analysis corrections
-----------------------------

``pyar.property.calculate_angle(A, B, C)`` now measures the bond angle at the
middle atom B in degrees. Earlier code incorrectly measured at C. This corrects
internal-coordinate bond angles and donor-H-acceptor angles in the existing
hydrogen-bond analysis. Its distance and angle cutoffs are unchanged. Hydrogen-
bond distance matrices now retain fractional distances rather than truncating
them to integers.

Internal-coordinate torsions now pass Cartesian positions to the dihedral
calculator and return signed degrees for valid four-atom paths. Previously,
these paths failed because atom indices were passed instead of positions.

If earlier analyses used these routines, recompute the angle and hydrogen-bond
reports from the original geometries. Stored reports are not repaired
automatically; the corrections do not rewrite input geometries.

Undefined geometry is an error
------------------------------

The angle and dihedral helpers raise ``ValueError`` when a bond-angle ray has
zero length, the central torsion bond has zero length, or either torsion plane
is undefined because its bond vectors are collinear. Internal-coordinate and
hydrogen-bond analysis propagate these errors. Inspect coincident atoms or
linear four-atom paths before using those analyses; an undefined torsion is
not a zero-degree torsion.

Defined zero- and 180-degree bond angles and signed torsions retain their
existing conventions. No near-degeneracy tolerance or scientific cutoff is
introduced, and no stored calculation or geometry is rewritten.
