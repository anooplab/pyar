Sampling
========

PyAR trial geometry generation first samples directions on the unit sphere
and, when the monomer has more than one atom, samples monomer orientations
with quaternion-based rotations. This spatial sampling layer is separate from
backend optimization and is intentionally not a proposal-memory or adaptive
sampling system.

Direction Methods
-----------------

The supported direction methods are:

* ``fibonacci``: deterministic equal-area spherical Fibonacci sampling.
* ``lhs``: Latin-hypercube sampling on the sphere, seeded when requested.
* ``random``: random unit vectors on the sphere, seeded when requested.
* ``lhs_maximin``: oversampled Latin-hypercube candidates filtered by greedy
  maximin selection.
* ``fibonacci_maximin``: oversampled Fibonacci candidates filtered by greedy
  maximin selection.

The default direction method is ``fibonacci``.

Rotation Methods
----------------

The supported rotation/orientation methods are:

* ``halton``: deterministic low-discrepancy unit quaternions.
* ``random``: random unit quaternions, seeded when requested.

The default rotation method is ``halton``.

Benchmarking
------------

Run the coverage benchmark with:

.. code-block:: bash

   pyar-benchmark-orientations -N 8 12 20 28

The table reports:

* ``min_sep_deg``: smallest angular separation between sampled directions.
  Higher values reduce duplicate approaches.
* ``cover_radius_deg``: approximate worst unsampled angular gap, measured on
  a dense deterministic reference grid. Lower values indicate better
  full-sphere coverage.
* ``mean_cover_deg``: average angular distance to the closest sampled
  direction. Lower values are better.
* ``centroid_norm``: directional imbalance. Values near zero are preferred.

Determinism and Sequence Offsets
---------------------------------

The sampling layer is reproducible for fixed inputs:

* repeated calls to deterministic methods return identical results;
* seeded random methods return identical results for the same seed;
* ``sequence_offset`` advances the deterministic sampling sequence so
  restartable multi-population runs can generate distinct but reproducible
  point sets.

Trial-vector generation uses Fibonacci directions by default and Halton
quaternion rotations by default. The generated trial-vector array contains
three unit direction components followed by three Z-X-Z Euler angles in
radians. When the monomer is monatomic, the angle columns are zeroed because
the legacy geometry-merging path does not need rotational degrees of freedom.

``trial_vectors.dat`` is written as a stable whitespace-separated table with
one trial vector per row and six columns. The file is overwritten on each
write so repeated calls with the same input produce identical output.

The implementation is exposed by ``pyar.sampling`` together with the
``pyar-trial-generation`` command.

Literature Basis
----------------

Spherical Fibonacci point sets are established deterministic point sets for
nearly uniform sampling on the unit sphere and can be generated for arbitrary
numbers of points:

* R. Marques, C. Bouville, M. Ribardiere, L. P. Santos, and K. Bouatouch,
  *Spherical Fibonacci Point Sets for Illumination Integrals*, Computer
  Graphics Forum **32** (2013) 134-143, DOI: 10.1111/cgf.12190.
* B. Keinert, M. Innmann, M. Saenger, and M. Stamminger, *Spherical
  Fibonacci Mapping*, ACM Transactions on Graphics **34** (2015), article
  193, DOI: 10.1145/2816795.2818131.

The present benchmark focuses on separation and approximate covering radius
because PyAR needs directionally distinct approaches with no large unsampled
region of the sphere.
