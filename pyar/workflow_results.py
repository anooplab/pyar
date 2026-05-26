"""Structured result objects returned by public workflow entrypoints."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class WorkflowResult:
    """Common metadata for workflow outcomes."""

    workflow: str
    status: str
    run_directory: str
    state_path: str | None = None
    selected_paths: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self):
        """Return a plain dictionary for logging or serialization."""
        return {
            "workflow": self.workflow,
            "status": self.status,
            "run_directory": self.run_directory,
            "state_path": self.state_path,
            "selected_paths": list(self.selected_paths),
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class AggregateResult(WorkflowResult):
    """Structured result for aggregation runs."""


@dataclass(frozen=True)
class SolvationResult(WorkflowResult):
    """Structured result for solvation runs."""


@dataclass(frozen=True)
class ReactionResult(WorkflowResult):
    """Structured result for reaction runs."""
