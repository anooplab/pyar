geomeTRIC reaction-path validation
==================================

``pyar-neb`` runs an unbiased reaction-path workflow using a registered PyAR
energy/gradient backend. Endpoint and NEB optimizations use geomeTRIC; the TS
stage also uses geomeTRIC by default and can optionally use Sella. Atom count,
element order, and atom mapping must agree between reactant, product, and TS
waypoint files. Supply one finite XYZ frame per file.

The complete sequence is:

#. ``relax``: independently optimize both input endpoints without reaction bias.
   Both must converge and retain their covalent-radius connectivity. A temporary
   topology change in an adaptive trajectory is insufficient to launch this workflow.
#. ``neb``: construct an initial band through the supplied TS waypoint and
   optimize it. The default ``linear`` initialization uses the current
   piecewise Cartesian interpolation. Optional ``idpp`` initialization
   refines each side separately while keeping the supplied TS waypoint fixed.
   Optional ``geodesic`` initialization applies IDPP first, then smooths each
   side separately with redundant-internal-coordinate geodesic interpolation.
   In all cases, geomeTRIC optimizes the band; it must converge and its
   highest-energy image must be an interior image.
#. ``ts``: optimize that image as a first-order saddle candidate. geomeTRIC is
   the default optimizer; an optional Sella saddle optimizer can be selected.
#. ``frequency``: calculate a fresh finite-difference Cartesian Hessian and
   verify stationarity and exactly one significant imaginary frequency.
#. ``irc``: trace forward and backward directions separately, recording convergence
   for each branch. A converged backward branch cannot conceal a failed forward branch.
#. ``endpoints``: optimize both IRC endpoints again, then calculate their Hessians
   and frequencies. Confirm the reaction only if both IRC branches converged,
   both endpoints are minima, and they match the two relaxed references one-to-one.

Final endpoint optimization is unconstrained: it may change connectivity. Matching
uses the resulting minima, covalent-radius graphs, and mapped RMSD after rigid
alignment. Two branches ending at the same minimum cannot confirm a reaction.
The two optimized IRC endpoints must themselves differ in connectivity or by
more than the mapped RMSD tolerance, even when both fit the reference tolerances.
These checks provide numerical evidence; inspect chemical identities and the
appropriateness of the chosen electronic-structure method as well.

Complete and separate runs
--------------------------

For the bundled HCN/HNC isomerization example::

   pyar-neb --start tests/data/neb/hcn.xyz --end tests/data/neb/hnc.xyz \
     --ts-guess tests/data/neb/guess.xyz --software xtb --max-cycles 200 \
     --output hcn_hnc

``--stage all`` is the default. PyAR's ``xtb`` gradient provider defaults to
g-xTB. Select its Hamiltonian explicitly with ``--xtb-model``:

* ``--xtb-model gxtb`` runs ``xtb input.xyz --gxtb --grad`` (the default, for
  backward compatibility).
* ``--xtb-model gfn2`` runs ``xtb input.xyz --gfn 2 --grad``.

The executable must advertise ``--gxtb`` for g-xTB. If it does not, PyAR stops
with an error rather than silently calculating with the executable's default
Hamiltonian.

The selected model is part of physical restart provenance, so changing it
invalidates stages calculated with the other Hamiltonian. The explicit GFN2-xTB
selection allows controlled calculations on the same Hamiltonian family used
by GFN2-xTB benchmark datasets such as RGD1-TSopt-GFN2. This selects the model
through the ``xtb`` executable; it does not claim numerical identity with other
xTB implementations such as tblite. ``--method`` and ``--basis`` do not select
the model for the xTB energy-gradient provider. For other backends, use their
supported method/basis settings consistently; ``--xtb-model`` has no effect.

Each stage is independently runnable. Use the same output directory and physical
settings for consecutive stages::

   pyar-neb --stage relax --start hcn.xyz --end hnc.xyz --software xtb --output hcn_hnc
   pyar-neb --stage neb --ts-guess guess.xyz --software xtb --output hcn_hnc
   pyar-neb --stage ts --software xtb --output hcn_hnc
   pyar-neb --stage frequency --software xtb --output hcn_hnc
   pyar-neb --stage irc --software xtb --output hcn_hnc
   pyar-neb --stage endpoints --software xtb --output hcn_hnc

``ts --ts-geometry guess.xyz`` can optimize an independent TS guess without
endpoint or NEB artifacts. ``frequency --ts-geometry optimized.xyz`` can verify
an independently optimized TS. IRC consumes the verified frequency geometry
and Hessian; final endpoint comparison additionally needs the ``relax`` artifacts.

Artifact checks and convergence
-------------------------------

Every stage writes ``<stage>_summary.json`` with method settings, convergence
evidence, stage parameters, input/output hashes, and hashes of dependency
summaries. ``workflow_summary.json`` records
the full run. Stages reject modified geometries/Hessians, failed prerequisites,
and mismatched physical settings. After replacing an upstream structure, rerun
its dependent stages. Changing an upstream validation threshold or result also
invalidates dependent stages, even if the geometry and Hessian are unchanged.
Dependency acceptance conditions are rechecked on loading. Schema version 1
summaries can be reused without repeating their calculations by adding
``--reuse-legacy-summaries`` to each staged command that consumes them. PyAR
checks their recorded method, input and output hashes, dependency chain, and
scientific gates, then upgrades the summaries in place. Since version 1 did not
record stage-specific options, migrated summaries are explicitly marked with
unverified legacy parameters; new stage summaries record their options. Without
the flag, PyAR refuses legacy reuse and directs you to rerun that stage. The
flag is for staged commands that consume prior summaries; ``--stage all`` runs
the full workflow. Current summaries use version 2.
Process count may change between stages. Separate stages leave the historical
``workflow_summary.json`` unchanged; consult the current stage summary.

``--images`` must be odd and at least three. ``--max-cycles``, ``--ts-max-cycles``,
``--irc-max-cycles`` (per branch), and ``--endpoint-max-cycles`` control the
respective iteration limits. The default imaginary-frequency threshold is
20 cm⁻¹; smaller negative frequencies remain in the output but do not count
toward the Hessian index. Frequency validation also requires maximum and RMS
atomic gradient norms below 4.5e-4 and 3.0e-4 Hartree/Bohr, respectively.

NEB initialization
------------------

``--interpolation linear`` is the default for backward compatibility. It
constructs the same piecewise Cartesian path through the supplied TS waypoint
as previous PyAR versions. ``--interpolation idpp`` asks ASE's IDPP
interpolator to refine the reactant-to-waypoint and waypoint-to-product path
halves independently. The relaxed endpoints and supplied waypoint remain fixed;
IDPP only prepares the initial images and does not replace geomeTRIC NEB or
establish a minimum-energy path. ``--idpp-fmax`` and ``--idpp-steps`` control
ASE's IDPP initialization (defaults 0.1 and 100). These settings are recorded
in the NEB stage summary and must match when reusing that stage.

``--interpolation geodesic`` additionally smooths the deterministic IDPP seed
on each side of the supplied waypoint using the optional
``geodesic-interpolate`` package. Install it with
``pip install 'pyar-chem[geodesic]'``. ``--geodesic-tol`` and
``--geodesic-max-iter`` control smoothing (defaults 0.002 and 15). PyAR keeps
the requested image count and restores each smoothed image to the rigid-body
frame of its seed; the relaxed endpoints and supplied waypoint remain exact.
The package's random ``redistribute()`` initializer is not used. Geodesic
initialization only prepares coordinates: using a supplied TS waypoint does
not verify it, and smoothing does not establish a minimum-energy path,
transition state, converged reaction path, or barrier. The installed package
version is recorded for provenance. See the
`geodesic-interpolate project <https://github.com/virtualzx-nad/geodesic-interpolate>`_.

TS optimization
---------------

``--ts-optimizer geometric`` is the default and preserves the existing
geomeTRIC TS optimization. Optionally, install ``pip install 'pyar-chem[sella]'``
and select ``--ts-optimizer sella`` to run Sella's first-order saddle optimizer
with the same unbiased PyAR energy/gradient calculator. ``--sella-fmax`` sets
Sella's force convergence threshold in eV/angstrom (default 0.05), and
``--ts-max-cycles`` sets its maximum optimizer steps. These controls are recorded
only when Sella is selected; changing inactive Sella options does not invalidate
a geomeTRIC TS stage.

Sella convergence records optimizer convergence only. It does not confirm a
transition state or replace the independent frequency validation, which remains
the source of ``first_order_saddle_confirmed``. Completed TS artifacts can be
used by frequency, IRC, and endpoint stages without Sella installed. The Sella
method is described by `Hermes et al., J. Chem. Theory Comput. 2022
<https://doi.org/10.1021/acs.jctc.2c00395>`_.

Sub-1e-5 angstrom deviations from a best-fit line are removed before frequency
evaluation, with the correction recorded. Both the gradient and Hessian are
evaluated at that corrected geometry. This avoids losing a bending mode at a
numerically bent linear minimum. Final minima require zero significant imaginary
frequencies and successful optimization. The default endpoint RMSD tolerance is
0.5 angstrom (``--irc-endpoint-rmsd-tolerance``).

Outputs include ``neb_path.xyz``, ``ts_optimized.xyz``, ``frequency_geometry.xyz``,
``ts_hessian.txt``, individual IRC paths, ``irc_path.xyz``, and
``irc_forward_relaxed.xyz`` / ``irc_backward_relaxed.xyz`` with their Hessians
and frequency files. Iteration-limit trajectories are retained, but do not count
as converged. The CLI exits nonzero when its requested scientific validation fails.
