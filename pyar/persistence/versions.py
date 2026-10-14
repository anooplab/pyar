"""Read-only layout detection; scientific reuse remains workflow-owned."""
from pathlib import Path

from .artifacts import validate_uuid
from .json_store import read_json
from .layout import WORKFLOWS


STATE_SCHEMAS = {"grow": 2, "microsolvation": 2, "aggregate": 2, "reaction": 3, "conformer": 4}


class LayoutError(ValueError):
    """Unknown, corrupt or ambiguous run metadata must not be mutated."""


def detect_version(root, workflow, *, validate_v1=None, validate_v2=None):
    """Return 0 (fresh), 1 (recognized legacy) or 2 (explicit new layout).

    Validators are read-only callbacks owned by workflows. A v1 validator must
    positively recognize and validate its records; filenames alone do not
    suffice. V2 detection checks identity/format, not scientific stage gates.
    """
    if workflow not in WORKFLOWS:
        raise LayoutError(f"Unsupported workflow: {workflow!r}")
    root = Path(root).resolve()
    if root.exists() and not root.is_dir():
        raise LayoutError(f"Run root is not a directory: {root}")
    marker = root / "layout.json"
    if marker.is_symlink():
        raise LayoutError("Layout marker must not be a symlink")
    if not marker.exists():
        if not root.exists() or not any(root.iterdir()):
            return 0
        if validate_v1 is not None and validate_v1(root):
            return 1
        raise LayoutError("Nonempty unrecognized run; use its original workflow or a fresh directory")
    try:
        data = read_json(marker)
        if (not isinstance(data, dict) or data.get("format") != "pyar-run-layout"
                or type(data.get("layout_version")) is not int or data["layout_version"] != 2
                or data.get("workflow") != workflow):
            raise ValueError("Unknown layout version/format or mismatched workflow")
        validate_uuid(data.get("run_id"))
        required = ["state.json"] if workflow == "conformer" else ["request.json"]
        if workflow not in {"scan-bond", "conformer"}:
            required.append("state.json")
        for filename in required:
            path = root / filename
            if path.is_symlink():
                raise ValueError(f"Authority must not be a symlink: {filename}")
            record = read_json(path)
            schema = "workflow_state_schema_version" if filename == "state.json" else "request_schema_version"
            expected_schema = STATE_SCHEMAS[workflow] if filename == "state.json" else 1
            if (not isinstance(record, dict) or "version" in record
                    or record.get("run_id") != data["run_id"]
                    or record.get("workflow") != workflow
                    or type(record.get(schema)) is not int or record[schema] != expected_schema):
                raise ValueError(f"Mixed or mismatched authority: {filename}")
        if validate_v2 is not None:
            validate_v2(root)
    except (OSError, ValueError, TypeError) as exc:
        raise LayoutError(f"Invalid or incomplete v2 run; validate/recover explicitly: {exc}") from exc
    return 2
