# Contributing to PyAR

Use the developer guide in the documentation for architecture, canonical APIs,
scientific contracts and validation procedures. `AGENTS.md` contains repository
conventions. Keep changes scoped and preserve raw evidence for scientific claims.

Install an isolated development environment with
`python -m pip install -e ".[selection,xtb]" --group dev`.
Run focused tests for the affected behavior, then the applicable integration
checks. External backend tests require their actual programs/models; record
skips and avoid calling mocked behavior scientific validation.

Update command help, user examples and changelog with interface changes. Verify
`python tools/check_documentation.py` and a strict Sphinx build. Keep private
plans, operational reports and runtime outputs outside Git. Licensed fixtures,
reference data and provenance needed by reproducible tests must remain available.

Release preparation and publication have separate owner approval gates.
