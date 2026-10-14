"""Stable artifact references; hashes describe published bytes, never rank."""
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
from uuid import UUID, uuid4

from .layout import relative_path


class ArtifactError(ValueError):
    """An artifact cannot safely be located or its integrity is invalid."""


def validate_uuid(value):
    if not isinstance(value, str):
        raise ArtifactError("Identity must be a UUID string")
    try:
        parsed = UUID(value)
    except ValueError as exc:
        raise ArtifactError("Identity must be a UUID string") from exc
    if str(parsed) != value:
        raise ArtifactError("Identity must be a canonical UUID string")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class ArtifactRef:
    artifact_id: str
    role: str
    path: str
    sha256: str
    source_record_id: str | None = None

    def __post_init__(self):
        validate_uuid(self.artifact_id)
        if self.source_record_id is not None:
            validate_uuid(self.source_record_id)
        if not isinstance(self.role, str) or not self.role.strip():
            raise ArtifactError("Artifact role must be nonempty")
        try:
            relative_path(self.path)
        except ValueError as exc:
            raise ArtifactError(str(exc)) from exc
        if not isinstance(self.sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", self.sha256):
            raise ArtifactError("Artifact checksum must be a lowercase SHA-256 digest")

    def resolve(self, run_root):
        root = Path(run_root).resolve()
        target = root.joinpath(*relative_path(self.path).parts).resolve()
        if not target.is_relative_to(root):
            raise ArtifactError(f"Artifact escapes run root: {self.path}")
        return target

    def verify(self, run_root):
        path = self.resolve(run_root)
        try:
            digest = sha256_file(path)
        except OSError as exc:
            raise ArtifactError(f"Cannot read artifact {self.path}: {exc}") from exc
        if digest != self.sha256:
            raise ArtifactError(f"Artifact checksum mismatch: {self.path}")
        return path

    def to_dict(self):
        result = {"artifact_id": self.artifact_id, "role": self.role,
                  "path": self.path, "sha256": self.sha256}
        if self.source_record_id is not None:
            result["source_record_id"] = self.source_record_id
        return result

    @classmethod
    def from_dict(cls, value):
        try:
            return cls(**value)
        except (TypeError, ValueError) as exc:
            raise ArtifactError(f"Invalid artifact reference: {exc}") from exc

    @classmethod
    def capture(cls, run_root, path, role, *, artifact_id=None, source_record_id=None):
        """Assign once and persist the returned reference before reuse."""
        probe = cls(artifact_id or str(uuid4()), role, path, "0" * 64, source_record_id)
        target = probe.resolve(run_root)
        try:
            digest = sha256_file(target)
        except OSError as exc:
            raise ArtifactError(f"Cannot capture artifact {path}: {exc}") from exc
        return cls(probe.artifact_id, role, path, digest, source_record_id)
