Trial Direction Sampling
========================

Trial geometry generation begins by choosing directions on the unit sphere.
This spatial sampling problem is independent of the element types and should
be assessed separately from molecular rotation and local optimization.

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

Methods
-------

``fibonacci`` is deterministic equal-area spherical Fibonacci sampling.
``lhs`` and ``random`` are randomized baselines. ``lhs_maximin`` and
``fibonacci_maximin`` select a farthest-point subset from an oversized
candidate pool.

The aggregate workflow default uses spherical Fibonacci placement directions.
For multi-atom monomers it pairs those directions with deterministic
low-discrepancy quaternion rotations, converted to the existing Euler-angle
rotation API.

Standalone commands that request multiple generated populations use a
deterministic ``sequence_offset`` for each population. This retains
reproducibility while avoiding repeated copies of the same point design.

The implementation is exposed by ``pyar.trial_generation`` and the
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
