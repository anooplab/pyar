"""Atomic, request-validated restart state for fixed-seed growth."""
import json
import hashlib
import os
import tempfile
from pathlib import Path
from pyar.state.aggregate import _json_value


class GrowStateError(RuntimeError):
    """An existing growth run cannot safely be reused."""


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".state-",
                                         suffix=".json", delete=False) as fp:
            temporary = fp.name
            json.dump(_json_value(data), fp, indent=2, sort_keys=True, allow_nan=False)
            fp.write("\n")
            fp.flush()
            os.fsync(fp.fileno())
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)


class GrowRunState:
    def __init__(self, directory, data):
        self.directory = Path(directory).resolve()
        self.data = data
        self.path = self.directory / "state.json"

    @classmethod
    def load(cls, directory, request):
        directory = Path(directory).resolve()
        path = directory / "state.json"
        if not path.exists():
            if directory.exists() and any(directory.iterdir()):
                raise GrowStateError(f"Non-empty growth directory without state: {directory}")
            return None
        try:
            data = json.loads(path.read_text())
        except (ValueError, OSError) as exc:
            raise GrowStateError(f"Cannot read growth state {path}: {exc}") from exc
        if data.get("version") != 1 or data.get("workflow") != "grow":
            raise GrowStateError("Unsupported growth state format")
        if data.get("request") != _json_value(request):
            raise GrowStateError("Growth request differs from saved state; use a fresh output directory")
        completed = data.get("completed_steps", [])
        if ([item.get("step") for item in completed] != list(range(1, len(completed) + 1))
                or len(completed) > request["count"] or data.get("next_step") != len(completed) + 1
                or data.get("status") not in {"running", "completed", "no_candidates"}
                or (data.get("status") == "completed" and len(completed) != request["count"])):
            raise GrowStateError("Invalid growth progress in saved state")
        if completed:
            if [r["path"] for r in data.get("current_seeds", [])] != completed[-1]["selected_paths"]:
                raise GrowStateError("Growth current pool does not match the last completed step")
        for step in completed:
            if step["selected_count"] != len(step["selected_paths"]) or step["selected_count"] > request["maximum_number_of_seeds"]:
                raise GrowStateError("Invalid saved growth survivor budget")
            for saved in step["selected_paths"]:
                path = (directory / saved).resolve()
                if not path.is_relative_to(directory) or not path.is_file():
                    raise GrowStateError(f"Missing completed growth stage: {path}")
        for reference in data.get("current_seeds", []):
            snapshot = (directory / reference["path"]).resolve()
            if not snapshot.is_relative_to(directory) or not snapshot.is_file():
                raise GrowStateError(f"Missing or invalid growth snapshot: {snapshot}")
            if hashlib.sha256(snapshot.read_bytes()).hexdigest() != reference.get("sha256"):
                raise GrowStateError(f"Growth snapshot was modified: {snapshot}")
        return cls(directory, data)

    def save(self):
        write_json(self.path, self.data)
