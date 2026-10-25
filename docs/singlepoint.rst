Single-point energy, gradients and frequencies
==============================================

Use ``pyar energy``, ``pyar engrad`` and ``pyar frequency`` on one or more XYZ
files. Each input is evaluated independently at its supplied coordinates;
these commands do not optimize or add reaction bias. ``gradient`` and ``engra``
are aliases of ``engrad``. ``pyar energies`` instead reads energies already
stored in XYZ files.

Short commands and defaults
---------------------------

.. code-block:: bash

   pyar energy *.xyz
   pyar engrad *.xyz
   pyar frequency *.xyz

The default backend is external xTB, using GFN2-xTB, one processor, accuracy
0.1 and at most 250 electronic iterations. Charge defaults to zero;
multiplicity follows the same electron-parity singlet/doublet inference as
``pyar optimize``. Specify the intended spin explicitly for other states.

All inputs are preflighted before calculations begin. Duplicate basenames are
rejected, including files from different directories. A runtime failure is
recorded and subsequent independent inputs continue. The default output root
is ``energy/``, ``engrad/`` or ``frequency/``. Use a new ``--output`` root for
each level or solvent; existing roots are never overwritten or used as a cache.

.. code-block:: text

   energy/
     request.json
     summary.json
     summary.csv
     pyar-run.toml
     job_molecule1/
       input.xyz
       request.json
       result.json
       backend/reference/
         backend.out
         backend.err
         backend_command.json
         backend_exit.json
         evaluation.json
     job_molecule2/
       ...

Gradient jobs also write ``gradient.txt``. Frequency jobs write ``hessian.txt``
and ``normal_modes.json`` and retain every displaced evaluation under
``backend/displacement_NNNN_plus/`` and ``backend/displacement_NNNN_minus/``.
Python calculators retain their worker request, result and output streams;
Turbomole retains ASE's native files. Results contain geometry hashes, source
file hashes, electronic state, resolved theory and units. Native program
versions are extracted where available, with raw output retained; ML metadata
includes package/model identity and AIMNet2 model SHA256. An arbitrary XYZ's
optimization level is unknown and is never inferred from its comment or energy.

Choosing a backend and theory
-----------------------------

.. list-table:: Standalone support (separate from optimizer support)
   :header-rows: 1

   * - Backend
     - Defaults and controls
   * - ``xtb``
     - ``--xtb-model gfn2``; also gfn0, gfn1, gfnff, gxtb when the executable supports it.
   * - ``orca``
     - PBE/def2-SVP, TightSCF. Choose ``--method`` and ``--basis``.
   * - ``gaussian`` (or ``g16``)
     - PBE1PBE/def2-SVP. Supported SCF HF/DFT methods only; unsupported correlated methods are rejected.
   * - ``turbomole``
     - PBE/def2-SVP through ASE's conventional dscf/grad adapter, one processor. Choose an ASE-supported functional or HF.
   * - ``aimnet_2``
     - Installed legacy AIMNet2 wB97M-D3 model; singlets only in this adapter. Requires the aimnet2 extra and model assets.
   * - ``ani`` (or ``torchani``)
     - ANI1x; choose ANI1ccx or ANI2x only if the installed TorchANI supplies it. Neutral singlets; requires the ml extra.

Other optimizer adapters are not automatically available for standalone
calculations. ``pyar backends --json`` reports standalone capabilities separately.
Gaussian and Turbomole inputs follow public documentation and have mocked
adapter tests; a locally unexecuted adapter is not a claim of real-backend
qualification. Unsupported options fail rather than being silently ignored.

ORCA composite and xTB methods define their own basis; omit ``--basis`` for
these. Additional dispersion defaults to ``none``; corrections intrinsic to a
chosen method remain present. ``--dispersion D3BJ``, ``D3ZERO`` or ``D4`` is
available for ORCA; Gaussian supports D3BJ and D3ZERO. No optimizer-specific RI,
ECP or extra dispersion choices are injected. Choose theory/basis settings
appropriate to the elements and electronic state.

Advanced controls include ``--nprocs``, ``--memory-mb`` (ORCA per process,
Gaussian total), ``--scf-cycles`` and xTB ``--accuracy``. Molecules execute
sequentially, with processors assigned within each calculation. See command
``--help`` for all options. Cancellation propagates; completed per-input results
and incremental summaries remain available.

Higher-level single points and solvent calculations
---------------------------------------------------

Optimize first using the intended geometry level, then supply the resulting
XYZ geometries to a new calculation root:

.. code-block:: bash

   pyar energy optimized/*.xyz --backend orca --method PBE0 --basis def2-TZVP --output higher_level
   pyar energy optimized/*.xyz --backend orca --method PBE0 --basis def2-TZVP --output gas
   pyar energy optimized/*.xyz --backend orca --method PBE0 --basis def2-TZVP --solvent-model CPCM --solvent water --output water

Solvent models are xTB ALPB/GBSA, ORCA CPCM/SMD, and Gaussian PCM/SMD. Backend
solvent names must be valid for the installed program. Supplying ``--solvent``
alone selects ALPB, CPCM or PCM respectively. The same settings apply to
gradients and numerical frequencies, including all displaced calculations.

A paired difference ``E(solvent) - E(gas)`` requires matching geometry,
composition, electronic state and compatible theory. These commands retain
both evaluations; they do not automatically pair results or claim a universal
solvation free energy. Solvent methods can include different non-electronic
terms. Likewise, combining higher-level electronic energies with lower-level
thermal corrections is a separate, explicit analysis; no composite free
energies are assembled automatically.

Frequency analysis and minimum qualification
--------------------------------------------

.. code-block:: bash

   pyar frequency optimized/*.xyz --backend orca --method PBE --basis def2-SVP
   pyar frequency optimized/*.xyz --step 0.005 --imaginary-threshold 20 --output tighter_check

Frequency analysis first evaluates energy and gradient, then builds a Cartesian
Hessian by central finite differences of gradients. The cost is ``6N + 1``
energy/gradient evaluations per N-atom geometry. No analytical Hessian is
claimed. ``--step`` is in bohr. The same theory and state apply to every
calculation. No input geometry is linearized or optimized.

Mass-weighted translations and rotations are projected out using SVD: three
rigid modes for an atom, five for a linear molecule and six for a nonlinear
molecule. Atomic masses and the observed rigid-mode count are recorded. All
vibrational frequencies, including small negative values, are retained in
cm^-1; negative values denote imaginary modes. The analysis is harmonic, with
no frequency scaling, isotope override or thermochemistry calculation.

A minimum requires both stationarity and no significant imaginary frequencies.
Defaults are a maximum per-atom gradient norm below ``4.5e-4`` Hartree/bohr,
RMS of per-atom norms below ``3e-4`` Hartree/bohr, and no frequency below
``-20`` cm^-1. These can be changed using ``--max-gradient``, ``--rms-gradient``
and ``--imaginary-threshold``; the selected criteria are recorded.

Results distinguish ``minimum_confirmed``, ``nonstationary``,
``first_order_saddle_candidate`` and ``higher_order_saddle_candidate``. A first-order
saddle candidate is not a validated reaction transition state. Near-zero and
small negative modes require review, particularly in floppy systems. An
``engrad`` job reports stationarity but always leaves the minimum unassessed.

Unconfirmed minima produce a warning against use as validated minima or for
minimum-based thermochemistry, and a nonzero frequency batch exit status.
Calculation completion is recorded separately from scientific qualification.
A geometry optimized at another level may not be stationary at the frequency
level; results are retained for diagnosis rather than silently reoptimized.

Preflight, configuration and tests
----------------------------------

.. code-block:: bash

   pyar energy *.xyz --check
   pyar frequency *.xyz --dry-run --json
   pyar engrad *.xyz --check --write-config engrad.toml
   pyar run engrad.toml

``--check`` and ``--dry-run`` do not create calculation directories. Normal
runs support profiles, saved configurations, machine-readable JSON and the
shared PyAR run record. A batch returns nonzero if any calculation failed;
frequency also returns nonzero when any completed geometry is not qualified
as a minimum. Results and raw backend evidence remain available.

Run deterministic adapter/science tests with:

.. code-block:: bash

   python -m pytest tests/test_singlepoint.py tests/test_singlepoint_backends.py
   PYAR_TEST_REAL_SINGLEPOINT=1 python -m pytest tests/test_singlepoint_backends.py -k real -v

The opt-in real tests require xTB, ORCA, AIMNet2 and ANI with their dependencies
and model assets. They locate H2 minima (water for ANI), check gradients against energy
differences, independently test frequency qualification, and exercise solvent
and charged doublet calculations for xTB/ORCA. They are small
contract checks, not scientific accuracy benchmarks.

Backend references
------------------

* `xTB single-point manual <https://xtb-docs.readthedocs.io/en/latest/sp.html>`_
* `ORCA 6.1 manual <https://www.faccts.de/docs/orca/6.1/manual/>`_
* `Gaussian Force <https://gaussian.com/force/>`_ and `SCRF <https://gaussian.com/scrf/>`_
* `Turbomole single-point modules <https://manual.turbomole.org/v8.0/how-to-run-turbomole/a-quick-and-dirty-tutorial-for-the-define-input-generator/>`_
