"""Versioned, inspectable restart state for aggregation workflows."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import numpy as np

from pyar import __version__

STATE_VERSION = 1
STATE_FILENAME = "state.json"


class AggregateStateError(RuntimeError):
    """Raised when an aggregation run cannot be safely resumed."""


def _json_value(value):
    """Return ``value`` in a JSON-serializable representation."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


class AggregateRunState:
    """Persist completed aggregation pathways independently of optimizer state."""

    def __init__(self, root_directory, data):
        self.root_directory = Path(root_directory).resolve()
        self.aggregate_directory = self.root_directory / "aggregates"
        self.state_file = self.aggregate_directory / STATE_FILENAME
        self.data = data

    @classmethod
    def create(cls, root_directory, request, pathway_labels, legacy_import=None):
        """Create and persist state for a new or imported aggregation run."""
        data = {
            "version": STATE_VERSION,
            "workflow": "aggregate",
            "package_version": __version__,
            "status": "running",
            "request": _json_value(request),
            "pathways": [str(label) for label in pathway_labels],
            "completed_pathways": [],
            "final_selected_results": [],
        }
        if legacy_import is not None:
            data["legacy_import"] = str(legacy_import)
        state = cls(root_directory, data)
        state.save()
        return state

    @classmethod
    def load(cls, root_directory, expected_request):
        """Load a running aggregation state file and validate its request."""
        state_file = Path(root_directory).resolve() / "aggregates" / STATE_FILENAME
        if not state_file.is_file():
            return None
        try:
            with state_file.open() as fp:
                data = json.load(fp)
        except (OSError, ValueError) as exc:
            raise AggregateStateError(
                f"Could not read aggregation state file {state_file}: {exc}"
            ) from exc
        if data.get("version") != STATE_VERSION or data.get("workflow") != "aggregate":
            raise AggregateStateError(
                f"Unsupported aggregation state format in {state_file}"
            )
        state = cls(root_directory, data)
        state.validate_progress()
        state.validate_request(expected_request)
        if data.get("status") != "running":
            raise AggregateStateError(
                f"Aggregation state in {state_file} is already {data.get('status')!r}; "
                "start a new calculation in a new directory."
            )
        return state

    def validate_request(self, expected_request):
        """Reject resume attempts that change the aggregation calculation."""
        if self.data.get("request") != _json_value(expected_request):
            raise AggregateStateError(
                "Existing aggregation state does not match this invocation; "
                "resume with the original inputs and settings or use a new directory."
            )

    def validate_progress(self):
        """Reject incomplete or inconsistent persisted pathway progress."""
        pathways = self.data.get("pathways")
        completed = self.data.get("completed_pathways")
        if not isinstance(pathways, list) or not isinstance(completed, list):
            raise AggregateStateError("Aggregation state has invalid pathway progress data.")
        if len(completed) > len(pathways):
            raise AggregateStateError("Aggregation state records more completed pathways than requested.")
        for index, record in enumerate(completed):
            if not isinstance(record, dict):
                raise AggregateStateError("Aggregation state has invalid completed pathway data.")
            if record.get("index") != index or record.get("label") != pathways[index]:
                raise AggregateStateError("Aggregation state completed pathways are out of sequence.")

    def completed_pathway_count(self):
        """Return the number of pathway builds already completed."""
        return len(self.data.get("completed_pathways", []))

    def remaining_pathway_labels(self):
        """Return persisted pathways that still need to be calculated."""
        completed = self.completed_pathway_count()
        return self.data["pathways"][completed:]

    def complete_pathway(self, label, selected_results):
        """Persist completion and selected result paths for the next pathway."""
        index = self.completed_pathway_count()
        if index >= len(self.data["pathways"]) or self.data["pathways"][index] != str(label):
            raise AggregateStateError(
                "Aggregation pathway completion is out of sequence; "
                "the saved state cannot be updated safely."
            )
        self.data["completed_pathways"].append(
            {
                "index": index,
                "label": str(label),
                "status": "completed",
                "selected_results": [str(path) for path in selected_results],
            }
        )
        self.save()

    def finish(self, final_selected_results):
        """Persist terminal workflow state and final selected output paths."""
        if self.remaining_pathway_labels():
            raise AggregateStateError(
                "Cannot complete aggregation state while pathways remain unfinished."
            )
        self.data["status"] = "completed"
        self.data["final_selected_results"] = [
            str(path) for path in final_selected_results
        ]
        self.save()

    def save(self):
        """Atomically write ``aggregates/state.json``."""
        self.aggregate_directory.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=self.aggregate_directory,
            prefix=".state-",
            suffix=".json",
            delete=False,
        ) as fp:
            json.dump(self.data, fp, indent=2, sort_keys=True)
            fp.write("\n")
            temporary_path = Path(fp.name)
        os.replace(temporary_path, self.state_file)
