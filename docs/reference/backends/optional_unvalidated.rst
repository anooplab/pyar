Optional backends without native validation
===========================================

MLatom/AIQM1, Gaussian and Turbomole remain optional, user-configured
integrations. Their native scientific calculations are **NOT_TESTED** for this
development milestone. Users with the required installations may try the
existing adapters; availability in PyAR is not a claim of validated numerical
results or suitability for a particular system.

These programs are not requirements of the core PyAR installation. Install and
license external software separately. Selecting another backend does not
require installing them. Existing unit, parser, synthetic-output and
failure-handling tests remain in place; they do not establish native scientific
validation. Historical evidence, where recorded, applies only to its recorded
versions, settings and calculations.

MLatom and AIQM1
----------------

Supply MLatom and its compatible Python dependencies, a licensed MNDO
installation for the AIQM1 semiempirical contribution, and the DFT-D4
executable. Configure MLatom's ``mndobin`` and ``dftd4bin`` environment variables
for those executables. Obtain MNDO access from its distributors; installing
MLatom does not grant an MNDO licence or install that program.

The current ``pyar mlopt`` driver explicitly uses Gaussian for optimization,
so it also requires a working, licensed Gaussian installation. The
``aiqm1_mlatom`` adapter invokes this driver; its ``mlatom_aiqm1`` alias does
not select an ASE optimizer. No optimizer migration is implied by installing
ASE. This milestone does not upgrade the ML dependency pins, change optimizer
selection or qualify a new combination of package versions.

Consult the `MLatom installation documentation
<https://aitomistic.com/mlatom/installation.html>`_ and the
`MNDO licence application <https://mndo.kofo.mpg.de/license.php>`_.

Gaussian
---------

Supply a licensed Gaussian installation, expose ``g16`` as required by the
selected PyAR route, and configure the vendor environment. See
:doc:`gaussian_backend` for the existing interface. Native execution remains
NOT_TESTED here; synthetic log and process tests protect software failure
contracts only.

Turbomole
----------

Supply a licensed Turbomole installation and configure its environment and
executables for the selected route. Native optimization uses ``define``,
``ridft``, ``rdgrad``, ``statpt`` and ``jobex``; other routes may require
``dscf`` and ``grad``. The ``xtb_turbo`` compatibility route additionally needs
xTB. Consult the `Turbomole documentation
<https://www.turbomole.org/turbomole/turbomole-documentation/>`_.

Current-run status markers, convergence checks, accepted geometry rules and
retained diagnostics still apply. User-managed installation does not bypass
those safeguards. Native Turbomole calculations remain NOT_TESTED here.

Trying an installation
----------------------

Use ``pyar install status`` and the selected command's ``--check`` option where
available to inspect prerequisites. These checks do not run or qualify the
scientific calculation and may not detect every external configuration issue.
Inspect raw native output, convergence, molecular state, units and resulting
geometry against an independently established calculation before relying on
results. Preserve software versions, calculation settings and artifacts when
reporting a failure.
