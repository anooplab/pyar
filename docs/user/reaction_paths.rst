Validate a reaction pathway with NEB, TS and IRC
================================================

Prepare matching ordered atoms in the two endpoint geometries and a reasonable
waypoint/path seed. Check charge, multiplicity, model and mapping. Download the
small illustrative :download:`HCN endpoint </examples/hcn.xyz>`,
:download:`HNC endpoint </examples/hnc.xyz>` and
:download:`waypoint </examples/guess.xyz>` to explore the interface.

.. code-block:: bash

   pyar neb hcn.xyz hnc.xyz --ts-guess guess.xyz --backend xtb --check
   pyar neb hcn.xyz hnc.xyz --ts-guess guess.xyz --backend xtb

The commands illustrate the available workflow; this HCN path was not newly
qualified during the documentation rewrite. Inspect each stage summary and
native output. A converged band with an endpoint maximum or no suitable climbing
image is insufficient TS evidence. The TS stage consumes the accepted NEB
candidate, then frequency/Hessian validation establishes the saddle gate.

``pyar ts`` operates on an explicit geometry or supported prior path evidence.
``pyar irc`` follows supported saddle evidence in both directions. Read each
command's help for stage inputs and output locations. A complete reaction
connection requires converged branches, accepted endpoint minima, matching atom
mapping, distinct references and the required identity/connectivity/geometry
checks. Optimizer success is not a substitute for these gates.

Do not edit stage hashes, reuse a stale Hessian after changing geometry, or treat
one successful IRC direction as a connection. Keep the request, stage parameters,
artifact dependencies and endpoint assignment with any reported barrier.
Electronic barrier reports are not activation free energies unless the necessary
thermochemical terms have been calculated consistently.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
