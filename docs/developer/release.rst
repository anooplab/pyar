Packaging and release procedures
================================

Validate built artifacts in fresh, noneditable environments outside the checkout.
Inspect wheel/sdist contents and metadata, licenses, optional extras and the sole
``pyar`` console entry. Check dependency consistency, canonical imports, root and
subcommand help, actionable missing dependencies and supported backend routes.

.. code-block:: bash

   python -m build

Run the required CI/engineering checks for the actual candidate and retain
source/artifact hashes and executed evidence. Scientific qualification is
scoped to its actual backend/model/protocol; package tests do not certify every
scientific workflow. Private release evidence and plans should remain outside
public documentation/package assets.

The repository's release skill and engineering gates govern release preparation.
A successful PR or build does not grant permission to tag, dispatch publication
or upload packages. Obtain explicit owner approval of version, commit and
artifacts for those actions. Read the Docs latest follows its configured source;
stable must identify an actually released version rather than an unpublished
candidate. Verify the hosted build and served pages independently of local
Sphinx success.

The :download:`mandatory release gates </release_engineering/gates.md>` define
the evidence required before publication approval.
