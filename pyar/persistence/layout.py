"""Immutable path construction without creating calculation output."""
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

WORKFLOWS = frozenset({"grow", "microsolvation", "aggregate", "reaction", "scan-bond", "conformer"})


def relative_path(value):
    """Require a portable, canonical run-relative path."""
    if not isinstance(value, str) or not value or "\\" in value or "\x00" in value:
        raise ValueError("Artifact path must be a nonempty POSIX relative path")
    path = PurePosixPath(value)
    if (path.is_absolute() or any(p in {"", ".", ".."} for p in value.split("/"))
            or ":" in path.parts[0]):
        raise ValueError(f"Unsafe run-relative path: {value!r}")
    return path


@dataclass(frozen=True)
class RunLayout:
    root: Path
    workflow: str

    def __post_init__(self):
        if self.workflow not in WORKFLOWS:
            raise ValueError(f"Unsupported v2 workflow: {self.workflow!r}")
        object.__setattr__(self, "root", Path(self.root).expanduser().resolve())

    def path(self, relative):
        """Resolve against the run root and reject symlink escapes."""
        target = self.root.joinpath(*relative_path(relative).parts)
        if not target.resolve().is_relative_to(self.root):
            raise ValueError(f"Artifact path escapes run root: {relative!r}")
        return target

    @property
    def marker(self):
        return self.path("layout.json")

    @property
    def inputs(self):
        return self.path("inputs")

    @property
    def results(self):
        return self.path("results")

    def step(self, index):
        if isinstance(index, bool) or not isinstance(index, int) or index < 0:
            raise ValueError("Step index must be a nonnegative integer")
        return self.path(f"steps/step_{index:03d}")
