# Workflow characterization evidence

Source commit: `9f1e567ce22b5d382d1facc9ff1da624de5d55a9` (PR94 merge).
These full JSON documents are recorded from the original checkout, not inferred
from the refactoring. `provenance.json` binds them to that source. Never update a
golden in a pure structural refactor.

The harness controls backend execution and methane embedding. It records calls
in order, including input-coordinate hashes, settings, electronic state and job
names. It captures every file: complete JSON/JSONL and CSV contents, exact XYZ
bytes/hashes, locks, diagnostics, immutable receipts and manifests. MS4 includes
proposal PCG64 states before/after, physical-evaluation reservations and results.
The fixtures establish software behavior; they do not validate physical models.
Existing independent analytical tests and real-backend comparisons remain required.

Only absolute run roots, execution timestamps/durations and lock PID/host are
normalized. UUID generation is deterministic with distinct valid identifiers;
identity/reference relationships and call ordering remain checked. Diagnostic
UUID filenames get stable role/name aliases. Artifact-reference hashes for JSON
are recalculated from the normalized referred documents, recursively. Hashes
embedded in CSV representations use the same mapping. Path-bearing optimization
request fingerprints are first verified against their original complete request,
then recomputed with only the root normalized. Scientific settings, coordinates,
thresholds, atom order, acceptance, counters, and every selected XYZ byte remain
exact. List order is preserved. No missing/failing cases are filtered out.

Reproduce from a clean original worktree (choose your own directory):

```sh
git worktree add --detach ../pyar-original 9f1e567ce22b5d382d1facc9ff1da624de5d55a9
python tools/record_workflow_goldens.py --baseline-root ../pyar-original --check --repeat 5
python -m pytest tests/test_workflow_golden.py tests/test_function_shape_measurement.py -q
python tools/measure_function_shape.py --revision 9f1e567ce22b5d382d1facc9ff1da624de5d55a9
```

Use the repository development environment, including the selection, RDKit and
identity capabilities used by these tests. The source checkout must be clean at
the stated commit; recording refuses a different revision or production edits.
The original/current trees are also reconciled across the separate explicit-CWD,
constants and typed-backend work; a discovered identity-log ownership regression
was fixed independently before decomposition, without rewriting these baselines.
