"""Conservative categories: existing composed tests are integration by default.

Scientific campaigns must opt in explicitly; benchmark-harness tests are not
scientific validation. All categories still execute in the full CI suite.
"""

import pytest

UNIT_FILES = {
    "test_persistence.py",
    "test_formula_generator.py",
    "test_coordinate_graph.py",
    "test_selection_policy.py",
    "test_selection_distances.py",
    "test_solvation_request.py",
    "test_conformer_request.py",
    "test_aggregation_request.py",
    "test_release_evidence.py",
}
REGRESSION_FILES = {
    "test_reaction_review_regressions.py",
    "test_modern_cli_review.py",
    "test_intramolecular_scan.py",
    "test_safe_recovery.py",
    "test_optimization_recovery_audit.py",
    "test_packaging_metadata.py",
    "test_installation_contract.py",
    "test_backend_failures.py",
}


def pytest_collection_modifyitems(items):
    for item in items:
        name = item.path.name
        if not any(
            item.get_closest_marker(m) for m in ("unit", "integration", "scientific")
        ):
            item.add_marker(
                pytest.mark.unit if name in UNIT_FILES else pytest.mark.integration
            )
        if name in REGRESSION_FILES:
            item.add_marker(pytest.mark.regression)
