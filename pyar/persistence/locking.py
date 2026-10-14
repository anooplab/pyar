"""Local POSIX single-writer guard; never replace/unlink a locked inode."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import socket
import stat
from uuid import uuid4

from .json_store import canonical_json_bytes, fsync_directory, _pairs, _nonfinite, _finite_float
from .artifacts import validate_uuid


class LockError(RuntimeError):
    """A run has an active writer or requires explicit stale-owner recovery."""


class RunLock:
    """Lock a canonical root's sibling file without creating the run itself.

    Supports local POSIX filesystems with flock semantics. Parent must exist.
    This is not a distributed lock or a guarantee about process-wide cwd.
    """
    def __init__(self, root, *, recover_stale=False):
        if type(recover_stale) is not bool:
            raise LockError("recover_stale must be an explicit boolean")
        self.root = Path(root).expanduser().resolve()
        if self.root == self.root.parent:
            raise LockError("Filesystem root cannot be a calculation root")
        self.path = self.root.parent / f".{self.root.name}.pyar-v2.lock"
        self.recover_stale = recover_stale
        self._descriptor = None
        self._owner = None

    def _write(self, record):
        payload = canonical_json_bytes(record)
        os.lseek(self._descriptor, 0, os.SEEK_SET)
        offset = 0
        while offset < len(payload):
            written = os.write(self._descriptor, payload[offset:])
            if written <= 0:
                raise OSError("Run lock metadata write made no progress")
            offset += written
        os.ftruncate(self._descriptor, len(payload))
        os.fsync(self._descriptor)

    def acquire(self):
        if self._descriptor is not None:
            raise LockError("This guard is already acquired")
        try:
            import fcntl
        except ImportError as exc:
            raise LockError("Run locking requires POSIX flock support") from exc
        if not hasattr(os, "O_NOFOLLOW"):
            raise LockError("Run locking requires no-follow file opening")
        try:
            descriptor = os.open(self.path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK, 0o600)
        except OSError as exc:
            raise LockError(f"Cannot open run lock {self.path}: {exc}") from exc
        self._descriptor = descriptor
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise LockError("Run lock must be a regular file")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise LockError(f"Run already has an active writer: {self.root}") from exc
            opened = os.fstat(descriptor)
            named = os.stat(self.path, follow_symlinks=False)
            if (opened.st_dev, opened.st_ino) != (named.st_dev, named.st_ino):
                raise LockError("Run lock inode changed; refusing unsafe acquisition")
            with os.fdopen(os.dup(descriptor), "rb") as stream:
                raw = stream.read()
            previous = None
            stale = False
            if raw:
                try:
                    previous = json.loads(raw, object_pairs_hook=_pairs, parse_constant=_nonfinite, parse_float=_finite_float)
                    if (not isinstance(previous, dict) or previous.get("format") != "pyar-run-lock"
                            or type(previous.get("lock_schema_version")) is not int
                            or previous["lock_schema_version"] != 1
                            or previous.get("run_root") != str(self.root)
                            or type(previous.get("pid")) is not int or previous["pid"] < 1
                            or not isinstance(previous.get("host"), str)
                            or not isinstance(previous.get("acquired_utc"), str)
                            or type(previous.get("active")) is not bool):
                        raise ValueError("Invalid owner record")
                    validate_uuid(previous.get("owner_id"))
                    stale = previous["active"]
                except (ValueError, UnicodeError):
                    stale = True
            if stale and not self.recover_stale:
                raise LockError("Stale/corrupt owner metadata; explicitly recover after obtaining the OS lock")
            now = datetime.now(timezone.utc).isoformat()
            owner = {"format": "pyar-run-lock", "lock_schema_version": 1,
                     "active": True, "run_root": str(self.root),
                     "owner_id": str(uuid4()), "pid": os.getpid(),
                     "host": socket.gethostname(), "acquired_utc": now}
            if isinstance(previous, dict) and "recovery" in previous:
                owner["recovery"] = previous["recovery"]
            if stale:
                owner["recovery"] = {"recovered_utc": now,
                                     "previous_metadata_sha256": hashlib.sha256(raw).hexdigest(),
                                     "previous_owner": ({key: value for key, value in previous.items()
                                                         if key != "recovery"}
                                                        if isinstance(previous, dict) else None)}
            self._owner = owner
            self._write(owner)
            fsync_directory(self.path.parent)
            return self
        except BaseException:
            # Closing releases a lock acquired by this descriptor. Do not
            # explicitly unlock a competing process's descriptor.
            os.close(descriptor)
            self._descriptor = None
            self._owner = None
            raise

    def release(self):
        if self._descriptor is None:
            return
        if self._owner["pid"] != os.getpid():
            raise LockError("A forked process cannot release its parent's guard")
        import fcntl
        try:
            self._write(dict(self._owner, active=False,
                             released_utc=datetime.now(timezone.utc).isoformat()))
        finally:
            fcntl.flock(self._descriptor, fcntl.LOCK_UN)
            os.close(self._descriptor)
            self._descriptor = None
            self._owner = None

    def __enter__(self):
        return self.acquire()

    def __exit__(self, *exception):
        self.release()
