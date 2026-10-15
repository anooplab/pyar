Run-layout persistence primitives
=================================

``pyar.persistence`` provides the building blocks used by fresh layout v2 runs.
Existing v1 trees retain their workflow owners and paths. See :doc:`run_layout`
for output discovery and compatibility. The architecture contract and baseline
inventory are in ``docs/run_layout_v2/``. The primitives do not implement a
generic scientific state machine.

Path and artifact references
----------------------------

``RunLayout(root, workflow)`` is immutable and creates no files or directories.
It accepts canonical workflow names (``grow``, ``microsolvation``, ``aggregate``,
``reaction``, ``scan-bond``, ``conformer``), without introducing CLI aliases.
``path(relative)`` rejects absolute paths, traversal, noncanonical separators,
Windows drive prefixes and symlink escapes. The run root is resolved once.

``ArtifactRef.capture(root, path, role)`` hashes the actual file bytes and assigns
an opaque UUID. Persist its dictionary once; use ``from_dict`` and ``verify``
when reusing it, rather than creating a new identity from a rank or basename.
A reference includes a role, canonical run-relative path, SHA-256 checksum and
optional source-record UUID. Missing or modified artifacts fail validation.

JSON publication
----------------

``write_json(path, value)`` serializes before creating output, using UTF-8,
sorted string keys, compact separators and one trailing LF. Nonfinite numbers
and non-string object keys are rejected. ``read_json`` rejects duplicate keys
and nonfinite values, including numeric overflow.

Publication uses a same-directory temporary file, file fsync, atomic replace,
and directory fsync. ``overwrite=False`` uses atomic no-clobber publication for
immutable records. The temporary file is removed on failure. An error after
publication (for example directory fsync failure) is propagated; recovery must
validate the published file. A successful file publication alone does not
establish scientific stage completion. Workflows must publish verified artifacts,
then completion records, then their mutable progress state.

Layout detection
----------------

``detect_version(root, workflow)`` reads only: 0 means absent/empty, 1 means a
legacy run positively recognized by the supplied ``validate_v1`` callback, and
2 means an explicit ``layout.json`` marker and matching versioned authorities.
Nonempty unrecognized output, malformed or unknown versions, mixed v1/v2
records, identity mismatches and incomplete initialization raise ``LayoutError``.
Scan recognition is delegated to its legacy request/result validator; it does
not require a fabricated ``state.json``. Scientific request and completion/gate
validation remain workflow-owned; a layout version is not permission to reuse
work. ``validate_v2`` may perform additional read-only workflow validation.

Local single-writer guard
-------------------------

``RunLock(root)`` construction is read-only. Acquiring it creates/opens the
sibling ``.<root-name>.pyar-v2.lock``; the parent must exist, and the calculation
root is not created. It uses nonblocking POSIX ``flock`` on a regular, non-symlink
file. A second process is rejected even with ``recover_stale=True``.
Never unlink or replace the lock file to recover an owner: inode identity is
part of mutual exclusion. Normal release marks the owner inactive. Process
death releases the OS lock but leaves active owner metadata; explicit
``recover_stale=True`` is then required. Recovery is recorded with the previous
metadata hash and owner. Corrupt owner metadata similarly requires explicit
recovery after acquiring the OS lock.

This contract is for local POSIX filesystems with reliable flock, hard-link,
atomic replace and fsync semantics. Network/distributed filesystem guarantees
are not provided. These APIs assume trusted run directories and one cooperating
writer; containment checks do not secure a directory concurrently modified by
an adversarial process. Locking does not make process-wide ``os.chdir`` safe
for threads. Read-only preflight must never acquire a writer lock.

Example for primitive validation
--------------------------------

.. code-block:: python

   from pyar.persistence import ArtifactRef, RunLayout, RunLock, write_json

   layout = RunLayout("grow", "grow")  # no output created
   with RunLock(layout.root):
       layout.inputs.mkdir(parents=True, exist_ok=True)
       seed = layout.path("inputs/seed.xyz")
       seed.write_text("1\ngeometry only\nHe 0 0 0\n")
       reference = ArtifactRef.capture(layout.root, "inputs/seed.xyz", "input")
       write_json(layout.path("inputs/metadata.json"), reference.to_dict(),
                  overwrite=False)
       reference.verify(layout.root)

This illustrates the primitives, not initialization of a complete v2 workflow
or an independent grow resume API. Workflow-owned loaders and scientific gates
remain necessary for safe continuation.
