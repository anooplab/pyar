Molecule
========

PyAR uses :class:`pyar.Molecule.Molecule` as its lightweight geometry container.
It stores the atomic symbols, Cartesian coordinates (Angstrom), and a small set
of workflow metadata (charge, multiplicity, fragments, energy).

XYZ IO
------

* :meth:`pyar.Molecule.Molecule.from_xyz` reads an ``.xyz`` file into a
  :class:`pyar.Molecule.Molecule`.
* :meth:`pyar.Molecule.Molecule.mol_to_xyz` writes the current structure out.

Internally, PyAR provides an exception-based parser :func:`pyar.Molecule.parse_xyz`
and preserves the historical fatal-exit behavior via :func:`pyar.Molecule.read_xyz`.

Geometry Transforms
-------------------

PyAR provides both in-place and non-mutating transform styles.

In-place (mutates the object)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* :meth:`pyar.Molecule.Molecule.translate`
* :meth:`pyar.Molecule.Molecule.rotate_3d`
* :meth:`pyar.Molecule.Molecule.align`

Non-mutating (returns a new object)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* :meth:`pyar.Molecule.Molecule.translated`
* :meth:`pyar.Molecule.Molecule.rotated_3d`
* :meth:`pyar.Molecule.Molecule.aligned`

ASE Interop (Optional)
----------------------

If ASE is installed, you can convert between PyAR and ASE objects:

* :meth:`pyar.Molecule.Molecule.to_ase_atoms`
* :meth:`pyar.Molecule.Molecule.from_ase_atoms`
