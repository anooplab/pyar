Build and rank molecular aggregates
===================================

Aggregation searches a requested final composition. Inputs are fragment species;
``--size`` gives their counts. Unlike ``grow``, no input is a permanently privileged
seed. Begin with a small orientation and pathway budget to check the chemistry.

Download :download:`water.xyz </examples/water.xyz>` and
:download:`sodium.xyz </examples/sodium.xyz>`.

.. code-block:: bash

   pyar aggregate water.xyz --size 2 --backend xtb --orientations 2 --maximum-number-of-seeds 2 --check
   pyar aggregate water.xyz --size 2 --backend xtb --orientations 2 --maximum-number-of-seeds 2

For a charged example:

.. code-block:: bash

   pyar aggregate sodium.xyz water.xyz --size 1 1 --charge 1 0 --multiplicity 1 1 --backend xtb --orientations 1 --maximum-number-of-seeds 1

The sodium–water example was exercised through native xTB for this change;
its backend command received charge +1. Each molecule carries the electronic
state passed centrally to optimization. Inspect actual backend output rather
than assuming an XYZ comment conveys the charge.

Results are published below ``aggregates/results``; pathway/step records show
which candidates survived and why. Composition, identity, connectivity policy
and energy comparability all matter. Use ``pyar energies`` on the selected XYZs,
then view representative structures. A tiny orientation count is an installation
or teaching example, not a converged global-minimum search. Evaluate sampling
and independent methodology before research conclusions.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
