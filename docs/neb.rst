geomeTRIC reaction-path validation
==================================

``pyar-neb`` runs an unbiased reaction-path workflow using a registered PyAR
energy/gradient backend. Every geometry optimization uses geomeTRIC. Atom count,
element order, and atom mapping must agree between reactant, product, and TS
waypoint files. Supply one finite XYZ frame per file.

The complete sequence is:

#. ``relax``: independently optimize both input endpoints without reaction bias.
   Both must converge and retain their covalent-radius connectivity. A temporary
   topology change in an adaptive trajectory is insufficient to launch this workflow.
#. ``neb``: interpolate through the TS waypoint and optimize the band. The band
   must converge and its highest-energy image must be an interior image.
#. ``ts``: optimize that image as a first-order saddle candidate.
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

``--stage all`` is the default. PyAR's ``xtb`` gradient provider currently invokes
``xtb --gxtb``; ``--method`` and ``--basis`` do not change that provider's model.
For other backends, use their supported method/basis settings consistently.

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
