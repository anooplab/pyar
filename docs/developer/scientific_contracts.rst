Scientific and backend contracts
================================

Ordinary software regression evidence and scientific validation are separate.
A mocked backend establishes control/data-flow behavior; native results require
raw executable/model/settings/output provenance and predefined acceptance.

Preserve atom order/mapping, fragment ownership, charge/spin, energy/gradient
units and signs, Hamiltonian/model, constraints and optimization objectives.
An accepted geometry must be the validated frame, not merely the last file a
program happened to write. Preserve physical/unbiased energy separately from
bias objectives. Check finite values, identity and genuine convergence.

Minimum qualification requires the implemented gradient and vibrational rules.
A transition-state claim requires the first-order-saddle/Hessian protocol; a
reaction connection also requires both IRC directions and qualified endpoint
assignment. Near-linear mode handling, thresholds and endpoint identity are
scientific policy. Structural refactoring does not authorize changes to them.

Caches must bind to current request, ordered geometry, electronic state,
provider/model version and accepted result evidence. Checkpoints must bind to
a compatible supported protocol and verified artifacts. Exit codes, optimizer
booleans, stale XYZs and manually edited summaries do not establish authority.
Recovery keeps its bounded approved policy and owner chemistry gates.

Use approved analytic/reference benchmarks and their recorded units/tolerances.
Do not invent new reference values from a single successful run. Retain raw
failure evidence; label unexecuted required calculations NOT_TESTED or BLOCKED
as appropriate. Numerical changes need explicit scientific review.

Numerical tolerances
----------------------

Every numerical benchmark case must state its metric, unit, expected value and
nonnegative tolerance in a frozen manifest. There is no universal cross-method
tolerance. Keep absolute and relative tolerances separate when near-zero values
make relative error unstable. Record finite-difference step size and whether a
derivative comparison uses individual components or a norm.

Exact invariants (atom mapping, hashes, schema identity), analytic/numerical
derivatives, same-method replay, structural equivalence, reaction-path evidence
and stochastic search outcomes require different comparisons. Numerical replay
must name its backend/model/runtime; stochastic comparisons must specify seeds,
resources and outcome distributions. Do not impose byte identity on genuinely
stochastic chemistry or replace exact provenance checks with broad tolerances.

Reference provenance and catalogue
-----------------------------------

Before a case is release eligible, preserve its source/version, DOI or URL,
license, acquisition record, exact input ordering, charge, multiplicity,
backend/method versions, settings, reference values, hashes and an independently
specified validation procedure. A reference update requires scientific rationale,
old/new hashes and owner review. Current PyAR output cannot be its own independent
reference. Retain external data outside Git when redistribution is not licensed.

The machine-readable inventory is
``benchmarks/scientific_validation/catalogue.json``. Its entries distinguish
constructed structural invariants, externally optimized conformer ensembles,
external RGD1 transition states, microsolvation protocols and bounded scan
interface evidence. Read each entry's limitations and eligibility; none alone
qualifies every PyAR workflow. Archived local reruns are operational evidence,
not replacements for the licensed inputs/manifests and source attribution kept
in the repository.
