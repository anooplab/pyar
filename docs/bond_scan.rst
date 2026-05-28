Bond Scan
=========

Use a bond scan when you want a simple structural probe along one distance
coordinate between two fragments. This is not a full reaction-path search, but
it is useful for quickly checking whether a chosen bond-forming or
bond-breaking coordinate behaves sensibly before running a more expensive
reaction search.

Typical chemistry questions
---------------------------

Bond scans are useful when you want to ask questions such as:

* What happens as two selected atoms are brought closer together?
* Does a close-contact structure look chemically reasonable before a full AFIR
  run?
* Which distance range should be explored more carefully?
* Is a simple coordinate scan enough, or is a full reaction search needed?

Basic command
-------------

Run a scan between atom ``1`` in the first structure and atom ``2`` in the
second structure:

.. code-block:: bash

   pyar-cli scan-bond 1 2 A.xyz B.xyz -N 8
   pyar-cli --scan-bond 1 2 A.xyz B.xyz -N 8

The two numbers identify the atoms involved in the distance scan. The two XYZ
files provide the starting fragments.

How to use the results
----------------------

A bond scan should be treated as a preliminary probe. If the scan suggests a
chemically interesting close-contact region, use the resulting structures as
starting points for a reaction search, higher-level optimisation, or manual
mechanistic analysis.

A scan is not a substitute for transition-state confirmation. It does not, by
itself, prove a minimum-energy path or a transition state.
