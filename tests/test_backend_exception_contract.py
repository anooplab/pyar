"""Typed failures are independent of optional backend installations."""

import ast
from pathlib import Path

import pytest

from pyar.backend_errors import (
    backend_operation,
    BackendExecutionError,
    BackendOutputError,
    BackendConvergenceError,
    BackendSCFError,
    OptimizationOutcome,
)


@pytest.mark.parametrize(
    "status", [False, None, "SCFFailed", "GradFailed", "UpdateFailed", "CycleExceeded"]
)
def test_untyped_failure_values_cannot_escape_adapter_boundary(status):
    @backend_operation
    def execute():
        return status

    with pytest.raises(BackendOutputError):
        execute()


def test_unexpected_adapter_error_preserves_cause():
    failure = ValueError("malformed backend output")

    @backend_operation
    def execute():
        raise failure

    with pytest.raises(BackendExecutionError) as caught:
        execute()
    assert caught.value.__cause__ is failure


def test_typed_adapter_error_retains_identity():
    failure = BackendSCFError("SCF rejected")

    @backend_operation
    def execute():
        raise failure

    with pytest.raises(BackendSCFError) as caught:
        execute()
    assert caught.value is failure


def test_budget_exhaustion_stays_provisional_at_optimization_owner(monkeypatch):
    import numpy as np
    from pyar import optimiser
    from pyar.core.molecule import Molecule

    molecule = Molecule(["H"], [[0, 0, 0]], name="incomplete")

    class Geometry:
        optimized_coordinates = np.zeros((1, 3))
        energy = -1.0

        def optimize(self):
            raise BackendConvergenceError("budget exhausted", usable_geometry=True)

    monkeypatch.setattr(optimiser, "build_geometry", lambda *args: Geometry())
    result = optimiser._backend_job(molecule, {})
    assert result.status is OptimizationOutcome.CYCLE_EXCEEDED
    assert not optimiser.is_success(result.status)
    assert optimiser.is_usable(result.status)
    assert result.energy == -1.0


def test_all_adapter_execution_methods_have_typed_boundaries():
    root = Path(__file__).parents[1] / "pyar"
    paths = [
        *sorted((root / "backends").glob("*.py")),
        root / "energy_gradient_providers.py",
    ]
    checked = []
    for path in paths:
        for cls in (
            node
            for node in ast.walk(ast.parse(path.read_text()))
            if isinstance(node, ast.ClassDef)
        ):
            if (
                cls.name == "EnergyGradientProvider"
            ):  # Structural protocol, no executable body.
                continue
            for method in cls.body:
                if isinstance(method, ast.FunctionDef) and method.name in (
                    "optimize",
                    "evaluate",
                    "run_turbomole_jobex",
                ):
                    assert any(
                        isinstance(decorator, ast.Name)
                        and decorator.id == "backend_operation"
                        for decorator in method.decorator_list
                    ), (path.name, cls.name, method.name)
                    checked.append((path.name, cls.name, method.name))
    assert checked  # The adapter inventory must not silently become empty.
