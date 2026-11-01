Architecture and ownership
==========================

The sole command entry is ``pyar.cli:main``. It loads task adapters lazily and
observes request preparation/execution/result through ``cli_runtime``. Parsers
translate user intent; request/service owners validate scientific settings;
workflows sequence sampling, execution, selection and persistence; backend
adapters translate physical values and interpret native evidence.

Shared electronic-state/preflight checks live in ``cli_validation`` independently
of the retired parser. ``cli_inputs`` builds task inputs. ``optimiser.optimise``
copies settings and injects the current molecule's charge/multiplicity/SCF type
and xTB unpaired-electron count before fingerprinting/backend construction.
Every backend attempt rewrites its trial XYZ from that request. Result-cache
reuse requires validated input/model/artifact provenance.

Scan orchestration owns route validation, chemistry preparation, immutable
request/recovery, orientation execution, result publication and finalization.
NEB separates ordered request validation, physical settings, cache adoption,
stage-specific execution, publication and failure receipts. Stage handlers
return result/input/artifact/dependency evidence to a single publication owner.
Recursive dependency/hash validation precedes acceptance of upstream evidence.

Current persistence owns supported series/layout identity, lock ownership,
immutable artifact references and workflow checkpoints. Derived reports are not
restart authority. Atomic mutable publication is distinct from an immutable
stage commit; neither creates a multi-file transaction by itself. Process-wide
working directories remain a resource requiring exception-safe restoration;
restoration does not make concurrent threads safe.

Keep scientific kernels and owner acceptance rules explicit. A helper is useful
when it clarifies responsibility and dependencies. Moving branches between
functions alone is not a reduction in aggregate scientific decisions.
