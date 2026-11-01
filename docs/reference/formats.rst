Output formats, status and provenance
=====================================

Supported workflow layouts contain an explicit marker and series identity,
request, immutable stage/selection records where applicable, checkpoints,
execution artifacts, and published results. Output filenames differ by task;
use the request/result or manifest references rather than assuming every
workflow has the same directory tree.

* Layout markers identify supported 1.6 runs and their owning workflow/run ID.
* Config/profile/run-record TOML uses the supported versioned schema. Versions
  are protocol identifiers, not numbers to edit to force compatibility.
* Artifact references include run-relative paths and checksums; state records
  keep the scientific request and the authority for accepted results.
* NEB stage records preserve stage settings, physical model, inputs, artifacts,
  dependency hashes and scientific gates. They do not accept older unverified
  parameter summaries.
* Standalone result JSON declares its units: coordinates in angstrom, energy in
  hartree and gradients in hartree/bohr. Backend/optimizer intermediate files
  can use other units and must be interpreted through their documented adapter.

A calculation can be completed but scientifically unqualified. Read both status
and qualification/gates. Nonzero exit may indicate invalid input, execution or
scientific qualification failure; inspect the emitted message and result evidence.
JSON output is for structured inspection, not permission to overwrite artifacts.
Keep raw backend output and original inputs with exported tables and plots.
