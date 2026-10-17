"""Context-local, non-mutating subprocess observations and artifact snapshots."""
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from uuid import NAMESPACE_URL, uuid5

from pyar.persistence import ArtifactRef, sha256_file
from pyar.persistence.layout import relative_path
from .models import ProcessRecord

_processes = ContextVar('pyar_optimization_process_observer', default=None)


@contextmanager
def capture_processes():
    records = []
    token = _processes.set(records)
    try:
        yield records
    finally:
        _processes.reset(token)


def record_process(command, returncode, *, exception=None):
    """No I/O, parsing or callbacks; legacy execution is unaffected when unused."""
    records = _processes.get()
    if records is not None:
        records.append(ProcessRecord(Path(command[0]).name, returncode, exception))


def safe_path(root, relative):
    root = Path(root).absolute()
    if root.is_symlink() or root.resolve() != root:
        raise ValueError('Observation root must not contain symlinks')
    path = root
    for part in relative_path(relative).parts:
        path = path / part
        if path.is_symlink():
            raise ValueError(f'Symlink artifact rejected: {relative}')
    return path


def snapshot_artifacts(root):
    """Deterministic derived refs; never assign or publish workflow artifact IDs."""
    root = Path(root).absolute()
    if root.is_symlink() or root.resolve() != root:
        raise ValueError('Observation root must not contain symlinks')
    if not root.exists():
        return ()
    refs = []
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root).as_posix()
        safe_path(root, relative)
        if path.is_file():
            digest = sha256_file(path)
            identity = str(uuid5(NAMESPACE_URL, 'pyar-ore-observation:' + relative + ':' + digest))
            refs.append(ArtifactRef(identity, 'diagnostic-evidence', relative, digest))
    return tuple(refs)
