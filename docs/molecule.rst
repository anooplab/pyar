Molecule
========

PyAR uses :class:`pyar.Molecule.Molecule` as its lightweight geometry container.
It stores the atomic symbols, Cartesian coordinates (Angstrom), and a small set
of workflow metadata (charge, multiplicity, fragments, energy).

Object Contract
---------------

The number of atom symbols must match the number of coordinate rows. Fragment
definitions, when provided, must refer to valid atom indices. Invalid
geometries fail at construction time with a clear ``ValueError``.

Two molecules can be combined explicitly with
:meth:`pyar.Molecule.Molecule.merged_with`. The ``+`` operator remains an
equivalent shorthand for existing PyAR workflows.

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

``rotate_3d`` rotates around the molecule centroid. ``align`` centers the
molecule at its center of mass before aligning it with its principal axes.
Reading :attr:`pyar.Molecule.Molecule.moments_of_inertia_tensor` does not
modify coordinates.

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
