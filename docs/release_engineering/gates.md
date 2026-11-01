# Mandatory release gates

These are current engineering procedures. Release preparation requires explicit evidence for the actual source commit and built artifact hashes.

- G0: intended candidate diff, consistent version/import/distribution metadata, changelog and declared support scope.
- G1: supported Python CI matrix, applicable coverage/negative/interruption/restart/optional-dependency checks; reviewed skips and expected failures.
- G2: strict documentation, validated commands, working downloadable examples and current public APIs.
- G3: fresh wheel and sdist installs, dependency consistency, sole pyar entry point, required licensed assets, no private plans or unrelated runtime output.
- G4: approved scientific matrix, reference provenance, predefined tolerances and real backend evidence for every mandatory advertised calculation route.
- G5: owner review of complete candidate evidence, exact source and artifact hashes, known limitations and support scope.
- G6: explicit owner approval before tag push, publish dispatch or package upload; verify protected publication environment and trusted publishing controls.
- G7: after authorized publication, registry version/hash verification, clean registry installation, release/tag/commit agreement and recorded post-release smoke results.

Each gate is PASS, FAIL, BLOCKED or NOT_RUN. Only executed evidence can pass. Missing tools, approved numerical references, licensed programs, credentials or assets block the relevant claim. Preserve failed attempts and raw logs. A successful exit code does not override a failed scientific gate.

Freeze mandatory cases and tolerances before evaluation. Record commit, environment, executable/model versions, inputs/settings, seed/thread counts, times, outcomes and artifact hashes. Never alter tolerances after seeing results or drop failed mandatory cases silently.

Candidate/artifact changes invalidate approval and applicable evidence. A green PR authorizes neither a release nor package publication. Readiness means ready for owner approval only after all required prepublication gates pass. Keep private release evidence outside public manuals and distributions.
