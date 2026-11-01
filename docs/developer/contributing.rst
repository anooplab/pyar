Contribute and validate changes
===============================

Install developer capabilities in an isolated source environment:

.. code-block:: bash

   python -m pip install -e ".[selection,xtb]" --group dev
   python -m pytest -q tests/test_optimizer_input_state.py tests/test_backend_trial_input.py
   python tools/check_documentation.py
   python -m sphinx -E -b html -W --keep-going docs docs/_build/html

Run focused tests for each change and broader checks when changed dependencies
or unresolved failures justify them. Native backend tests are marked separately
and require the executable/model. Keep missing capability/skip evidence visible.
The final integration CI matrix checks supported Python versions, packaging,
documentation, quality and selected optional integrations.

Keep regression witnesses tied to public behavior and scientific invariants.
Test malformed results, negative validation, interrupted attempts and supported
restarts as well as success. Preserve input state and original exceptions.
Do not weaken scientific gates to reduce complexity. The source-bound Radon
inventory and changed-code decision budget include extracted helpers; explain
ownership and aggregate decisions, not just a coordinator score.

Builds may regenerate the API page. Include generated changes only when
intentional. Never commit private plans, raw operational logs, credentials or
unrelated runtime output. Small licensed fixtures/reference evidence remain
necessary test assets and require provenance.
