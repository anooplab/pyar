"""Strict JSON and durable same-directory publication for new v2 records."""
import json
import math
import os
from pathlib import Path
import stat
import tempfile
from uuid import uuid4


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _nonfinite(value):
    raise ValueError(f"Nonfinite JSON number: {value}")


def _finite_float(value):
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"Nonfinite JSON number: {value}")
    return parsed


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream, object_pairs_hook=_pairs, parse_constant=_nonfinite, parse_float=_finite_float)


def canonical_json_bytes(value):
    # Reject non-string keys even when Python's encoder would silently coerce them.
    def validate(item):
        if isinstance(item, dict):
            if any(not isinstance(key, str) for key in item):
                raise ValueError("JSON object keys must be strings")
            for child in item.values():
                validate(child)
        elif isinstance(item, (list, tuple)):
            for child in item:
                validate(child)
    validate(value)
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def fsync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_json(path, value, *, overwrite=True):
    """Publish atomically; overwrite=False commits a new immutable record.

    A durability error after publication is propagated: callers must validate
    the file before recovery, rather than assuming the write never happened.
    """
    payload = canonical_json_bytes(value)
    _publish_bytes(path, payload, overwrite=overwrite, kind="JSON")


def write_bytes_atomic(path, payload):
    """Durably replace a mutable report while preserving its caller-owned bytes."""
    if not isinstance(payload, bytes):
        raise TypeError("Atomic report payload must be bytes")
    _publish_bytes(path, payload, overwrite=True, kind="report")


def _publish_bytes(path, payload, *, overwrite, kind):
    path = Path(path)
    if path.is_symlink():
        raise ValueError(f"Refusing {kind} publication through a symlink: {path}")
    missing = []
    parent = path.parent
    while not parent.exists():
        missing.append(parent)
        parent = parent.parent
    path.parent.mkdir(parents=True, exist_ok=True)
    # Persist newly created directory entries as well as the eventual file.
    for directory in reversed(missing):
        fsync_directory(directory.parent)
    temporary = None
    try:
        if kind == "JSON":
            stream_context = tempfile.NamedTemporaryFile(
                dir=path.parent, prefix=".pyar-json-", delete=False,
            )
        else:
            # Exclusive creation retains the mode of the previous text writer:
            # 0666 masked by the process umask. NamedTemporaryFile uses 0600.
            temporary = path.parent / f".pyar-report-{uuid4().hex}"
            stream_context = temporary.open("xb")
        with stream_context as stream:
            if kind == "JSON":
                temporary = Path(stream.name)
            elif path.exists():
                os.fchmod(stream.fileno(), stat.S_IMODE(path.stat().st_mode))
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            # Link is an atomic no-clobber publication; the temporary name is
            # removed below. This never links a mutable user result structure.
            os.link(temporary, path)
        temporary.unlink(missing_ok=True)
        fsync_directory(path.parent)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
