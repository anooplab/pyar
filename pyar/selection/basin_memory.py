"""Versioned geometry-backed memory for basin selection.

Schema 2 stores representative XYZ geometries as the source of truth and a
small versioned pair-distance descriptor for inspection/indexing. Legacy
fingerprints are migrated as opaque data: their historical implementation
could fall back to a different quantity, so they are never used to discard or
rank current candidates without a geometry to compare.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections import Counter
from types import SimpleNamespace

import numpy as np

from pyar.selection.features import _distance_histogram

REGISTRY_FORMAT = "pyar.basin-memory"
REGISTRY_VERSION = 2
HISTOGRAM_VERSION = 1
HISTOGRAM_BINS = 64
HISTOGRAM_MAX_DISTANCE_ANGSTROM = 12.0

__all__ = [
    "BasinMemoryError",
    "REGISTRY_FORMAT",
    "REGISTRY_VERSION",
    "_apply_basin_memory",
    "_basin_novelty_scores",
    "_basin_registry_path",
    "_entry_fingerprint",
    "_fingerprint_signature",
    "_load_basin_registry",
    "_persist_basin_registry",
    "_stoichiometry_label",
    "migrate_basin_registry",
    "record_selected_basins",
]


class BasinMemoryError(ValueError):
    """A basin registry is malformed or uses an unsupported schema."""


def _stoichiometry_label(molecule):
    """Return a compact stoichiometry label for basin snapshots."""
    counts = Counter(molecule.atoms_list)
    parts = []
    if "C" in counts:
        carbon = counts.pop("C")
        parts.append("C" if carbon == 1 else f"C{carbon}")
    if "H" in counts:
        hydrogen = counts.pop("H")
        parts.append("H" if hydrogen == 1 else f"H{hydrogen}")
    for element in sorted(counts):
        count = counts[element]
        parts.append(element if count == 1 else f"{element}{count}")
    return "".join(parts) if parts else "unknown"


def _basin_registry_path(molecule, output_root="selected", group_by_stoichiometry=True):
    """Return the persistence path for basin memory."""
    if group_by_stoichiometry:
        return os.path.join(output_root, f"stoichiometry_{_stoichiometry_label(molecule)}", "basin_registry.json")
    return os.path.join(output_root, "basin_registry.json")


def _fingerprint_signature(molecule):
    """Return the historical normalized Coulomb-spectrum signature.

    Kept only for compatibility with callers and schema-1 migration. New
    basin comparisons use stored geometries and the current distance policy.
    """
    from pyar import representations

    signature = representations.fingerprint(molecule.atoms_list, np.asarray(molecule.coordinates, dtype=float))
    signature = np.asarray(np.real_if_close(signature), dtype=float).ravel()
    norm = np.linalg.norm(signature)
    if norm > 0:
        signature = signature / norm
    return signature.tolist()


def _entry_fingerprint(entry):
    """Return a historical/legacy fingerprint vector when one is present."""
    try:
        if "fingerprint" in entry:
            return np.asarray(entry["fingerprint"], dtype=float).ravel()
        values = entry["legacy_descriptors"]["pyar-fingerprint-v1"]["values"]
        return np.asarray(values, dtype=float).ravel()
    except (KeyError, TypeError, ValueError):
        return None


def _valid_geometry(geometry):
    if not isinstance(geometry, dict):
        return False
    atoms = geometry.get("atoms")
    try:
        coords = np.asarray(geometry.get("coordinates"), dtype=float)
    except (TypeError, ValueError):
        return False
    return (
        isinstance(atoms, list)
        and bool(atoms)
        and coords.shape == (len(atoms), 3)
        and np.isfinite(coords).all()
        and all(isinstance(atom, str) and atom for atom in atoms)
    )


def _migrate_entry(entry):
    """Convert an unversioned entry while preserving its fingerprint exactly."""
    if not isinstance(entry, dict):
        return None
    migrated = {key: value for key, value in entry.items() if key not in {"fingerprint", "legacy_descriptors"}}
    if _valid_geometry(entry.get("geometry")):
        migrated["geometry"] = entry["geometry"]
    fingerprint = _entry_fingerprint(entry)
    legacy = dict(entry.get("legacy_descriptors", {})) if isinstance(entry.get("legacy_descriptors"), dict) else {}
    if fingerprint is not None and np.isfinite(fingerprint).all():
        legacy["pyar-fingerprint-v1"] = {
            "name": "pyar-fingerprint",
            "version": 1,
            "values": fingerprint.tolist(),
            "source_uncertain": True,
            "migration_note": "May be Coulomb eigenvalues or the historical numerical fallback; not comparable safely.",
        }
    if legacy:
        migrated["legacy_descriptors"] = legacy
    return migrated


def _normalize_payload(payload):
    """Return a validated schema-2 payload and migration metadata."""
    if isinstance(payload, list):
        raw_entries = payload
        old_version = 1
        stoichiometry = None
    elif isinstance(payload, dict):
        registry_format = payload.get("format")
        if registry_format not in (None, REGISTRY_FORMAT):
            raise BasinMemoryError(
                f"Unrecognized basin registry format {registry_format!r}; preserving it unchanged."
            )
        old_version = payload.get("schema_version", payload.get("version", 1))
        if isinstance(old_version, bool) or not isinstance(old_version, int):
            raise BasinMemoryError(f"Invalid basin registry schema version: {old_version!r}")
        if old_version > REGISTRY_VERSION:
            raise BasinMemoryError(
                f"Registry schema {old_version} is newer than supported schema {REGISTRY_VERSION}; preserving it unchanged."
            )
        if old_version < 1:
            raise BasinMemoryError(f"Invalid basin registry schema version: {old_version!r}")
        raw_entries = payload.get("entries")
        stoichiometry = payload.get("stoichiometry")
        if not isinstance(raw_entries, list):
            raise BasinMemoryError("Basin registry 'entries' must be a list; preserving it unchanged.")
    else:
        raise BasinMemoryError("Basin registry must contain a JSON object or legacy list; preserving it unchanged.")

    if any(not isinstance(entry, dict) for entry in raw_entries):
        raise BasinMemoryError("Basin registry contains a non-object entry; preserving it unchanged.")
    entries = []
    migrated = old_version != REGISTRY_VERSION
    for raw in raw_entries:
        if "geometry" in raw and not _valid_geometry(raw["geometry"]):
            raise BasinMemoryError("Basin registry contains malformed geometry; preserving it unchanged.")
        if "fingerprint" in raw:
            fingerprint = _entry_fingerprint(raw)
            if fingerprint is None or not np.isfinite(fingerprint).all():
                raise BasinMemoryError("Basin registry contains an invalid legacy fingerprint; preserving it unchanged.")
        if old_version < REGISTRY_VERSION or "fingerprint" in raw:
            entry = _migrate_entry(raw)
            migrated = True
        else:
            entry = raw
        entries.append(entry)
    return {
        "format": REGISTRY_FORMAT,
        "schema_version": REGISTRY_VERSION,
        "stoichiometry": stoichiometry,
        "entries": entries,
    }, migrated


def _read_payload(registry_path):
    def reject_nonstandard_constant(value):
        raise ValueError(f"non-standard JSON numeric constant {value}")

    try:
        with open(registry_path, "r", encoding="utf-8") as fp:
            raw = json.load(fp, parse_constant=reject_nonstandard_constant)
    except (OSError, ValueError) as exc:
        raise BasinMemoryError(f"Could not read basin registry {registry_path}: {exc}") from exc
    return _normalize_payload(raw)


def _load_basin_registry(registry_path):
    """Load registry entries, migrating older entry shapes in memory only."""
    if not registry_path or not os.path.exists(registry_path):
        return []
    payload, migrated = _read_payload(registry_path)
    if migrated:
        from pyar.selection import clustering

        legacy_count = sum(bool(entry.get("legacy_descriptors")) for entry in payload["entries"])
        clustering.cluster_logger.info(
            "Loaded legacy basin registry %s as schema %d (%d opaque legacy descriptors retained). "
            "Run pyar-basin-memory migrate to persist the schema migration.",
            registry_path, REGISTRY_VERSION, legacy_count,
        )
    return payload["entries"]


def _write_payload(registry_path, payload):
    directory = os.path.dirname(os.path.abspath(registry_path))
    os.makedirs(directory, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=directory,
            prefix=".basin-registry-", suffix=".json", delete=False,
        ) as fp:
            json.dump(payload, fp, indent=2, sort_keys=True, allow_nan=False)
            fp.write("\n")
            fp.flush()
            os.fsync(fp.fileno())
            temporary_path = fp.name
        os.replace(temporary_path, registry_path)
    finally:
        if temporary_path and os.path.exists(temporary_path):
            os.unlink(temporary_path)


def migrate_basin_registry(registry_path, *, dry_run=False):
    """Explicitly migrate a legacy registry to schema 2 without losing vectors.

    Old Coulomb/fallback vectors are retained as opaque, provenance-tagged
    descriptors. Migration cannot recreate geometries that were never stored.
    """
    if not os.path.isfile(registry_path):
        raise BasinMemoryError(f"Basin registry does not exist: {registry_path}")
    payload, migrated = _read_payload(registry_path)
    report = {
        "path": str(registry_path),
        "schema_version": REGISTRY_VERSION,
        "changed": migrated,
        "entries": len(payload["entries"]),
        "entries_with_geometry": sum(_valid_geometry(entry.get("geometry")) for entry in payload["entries"]),
        "entries_with_opaque_legacy_descriptor": sum(bool(entry.get("legacy_descriptors")) for entry in payload["entries"]),
    }
    if migrated and not dry_run:
        _write_payload(registry_path, payload)
    return report


def _geometry_descriptor(molecule):
    """Build the versioned, rigid-motion-invariant inspection descriptor."""
    atoms = [str(atom).capitalize() for atom in molecule.atoms_list]
    species = tuple(sorted(set(atoms)))
    descriptor = _distance_histogram(
        [molecule], species, bins=HISTOGRAM_BINS,
        maximum_distance=HISTOGRAM_MAX_DISTANCE_ANGSTROM,
    )[0]
    return {
        "name": "element-pair-distance-histogram",
        "version": HISTOGRAM_VERSION,
        "parameters": {
            "species": list(species),
            "bins": HISTOGRAM_BINS,
            "maximum_distance_angstrom": HISTOGRAM_MAX_DISTANCE_ANGSTROM,
            "normalization": "per-element-pair",
        },
        "values": descriptor.tolist(),
    }


def _geometry_record(molecule):
    coordinates = np.asarray(molecule.coordinates, dtype=float)
    if not _valid_geometry({"atoms": list(molecule.atoms_list), "coordinates": coordinates.tolist()}):
        raise BasinMemoryError(f"Cannot archive invalid geometry for {getattr(molecule, 'name', '<unnamed>')}")
    return {
        "atoms": [str(atom).capitalize() for atom in molecule.atoms_list],
        "coordinates": coordinates.tolist(),
        "charge": getattr(molecule, "charge", None),
        "multiplicity": getattr(molecule, "multiplicity", None),
    }


def _geometry_id(geometry):
    canonical = json.dumps(geometry, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _geometry_identity_signature(geometry):
    """Return a rigid-motion and atom-order invariant candidate key.

    This is only an index: entries sharing it still need a complete structural
    comparison before one can be discarded. Rounding limits false negatives
    from harmless coordinate serialization noise; it never proves identity.
    """
    atoms = geometry["atoms"]
    if len(atoms) < 2:
        # With no pair distances the index cannot safely distinguish anything.
        return None
    coordinates = np.asarray(geometry["coordinates"], dtype=float)
    pairs = []
    for left in range(len(atoms)):
        for right in range(left + 1, len(atoms)):
            pair_type = tuple(sorted((atoms[left], atoms[right])))
            distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
            pairs.append((pair_type, round(distance, 4)))
    return (tuple(sorted(Counter(atoms).items())), tuple(sorted(pairs)))


def _verified_duplicate(entry, geometry, comparator):
    """Only suppress an archive entry after a complete near-zero RMSD result."""
    previous = entry.get("geometry")
    if not _valid_geometry(previous):
        return False
    if (previous.get("charge"), previous.get("multiplicity")) != (
            geometry.get("charge"), geometry.get("multiplicity")):
        return False
    from types import SimpleNamespace

    def as_molecule(record):
        return SimpleNamespace(
            atoms_list=list(record["atoms"]),
            coordinates=np.asarray(record["coordinates"], dtype=float),
        )

    try:
        result = comparator.compare(as_molecule(previous), as_molecule(geometry))
    except Exception:
        return False
    return (result.compatible and result.metadata.get("comparison_complete") is True
            and result.distance is not None and np.isfinite(result.distance)
            and result.distance <= 1.0e-5)


def _persist_basin_registry(registry_path, selected_molecules, existing_entries=None, max_entries=200):
    """Atomically persist selected geometries and versioned descriptors."""
    if not registry_path or not selected_molecules:
        return
    if isinstance(max_entries, bool) or not isinstance(max_entries, int) or max_entries < 1:
        raise ValueError("max_entries must be a positive integer")
    if existing_entries is None:
        existing_entries = _load_basin_registry(registry_path)
    entries = [_migrate_entry(entry) for entry in existing_entries if isinstance(entry, dict)]
    entries = [entry for entry in entries if entry is not None]
    from pyar.structure_comparison import GraphFirstDeduplicationComparator

    comparator = GraphFirstDeduplicationComparator(threshold=1.0e-5, atom_mode="all")
    indexed_entries = {}
    for entry in entries:
        geometry = entry.get("geometry")
        if _valid_geometry(geometry) and _geometry_identity_signature(geometry) is not None:
            indexed_entries.setdefault(_geometry_identity_signature(geometry), []).append(entry)
    compacted = []
    for entry in entries:
        geometry = entry.get("geometry")
        signature = _geometry_identity_signature(geometry) if _valid_geometry(geometry) else None
        if signature is not None and any(
                candidate is not entry and _verified_duplicate(candidate, geometry, comparator)
                for candidate in indexed_entries[signature] if candidate in compacted):
            continue
        compacted.append(entry)
    entries = compacted
    seen_ids = {entry.get("geometry_id") for entry in entries if entry.get("geometry_id")}
    indexed_entries = {}
    for entry in entries:
        geometry = entry.get("geometry")
        if _valid_geometry(geometry) and _geometry_identity_signature(geometry) is not None:
            indexed_entries.setdefault(_geometry_identity_signature(geometry), []).append(entry)
    for molecule in selected_molecules:
        geometry = _geometry_record(molecule)
        geometry_id = _geometry_id(geometry)
        signature = _geometry_identity_signature(geometry)
        if geometry_id in seen_ids or any(
                _verified_duplicate(entry, geometry, comparator)
                for entry in indexed_entries.get(signature, ())):
            continue
        seen_ids.add(geometry_id)
        energy = getattr(molecule, "energy", None)
        try:
            energy = None if energy is None else float(energy)
        except (TypeError, ValueError):
            energy = None
        if energy is not None and not np.isfinite(energy):
            energy = None
        entries.append({
            "geometry_id": geometry_id,
            "name": str(getattr(molecule, "name", "")),
            "energy": energy,
            "geometry": geometry,
            "descriptors": {"geometry": _geometry_descriptor(molecule)},
        })
        indexed_entries.setdefault(signature, []).append(entries[-1])
    directory_label = os.path.basename(os.path.dirname(registry_path))
    stoichiometry = directory_label.removeprefix("stoichiometry_") if directory_label.startswith("stoichiometry_") else None
    payload = {
        "format": REGISTRY_FORMAT,
        "schema_version": REGISTRY_VERSION,
        "stoichiometry": stoichiometry,
        "entries": entries[-max_entries:],
    }
    _write_payload(registry_path, payload)
    from pyar.selection import clustering

    clustering.cluster_logger.info(
        "Updated geometry-backed basin registry with %d selected geometries at %s",
        len(selected_molecules), registry_path,
    )


def _as_molecule(entry):
    geometry = entry.get("geometry")
    if not _valid_geometry(geometry):
        return None
    return SimpleNamespace(
        name=entry.get("name", "basin-representative"),
        atoms_list=list(geometry["atoms"]),
        coordinates=np.asarray(geometry["coordinates"], dtype=float),
        energy=entry.get("energy", 0.0),
        charge=geometry.get("charge"),
        multiplicity=geometry.get("multiplicity"),
    )


def _geometry_novelty_scores(molecules, basin_entries, *, feature, distance_metric, system_type,
                             algorithm, distance_options, diagnostics=None):
    representatives = [
        molecule for molecule in (_as_molecule(entry) for entry in basin_entries)
        if molecule is not None
    ]
    if not representatives:
        return None
    from pyar.selection.features import compute_feature_matrix, standardize_features
    from pyar.selection.structural_distances import compute_distance_matrix

    combined = representatives + list(molecules)
    # Resolve ``auto`` from the current candidate pool first. Adding historical
    # representatives must not change the system classification that governs
    # this run; the resolved feature is then recomputed over the combined pool
    # so current and archived rows share one species vocabulary and scaling.
    candidate_features = compute_feature_matrix(
        molecules, feature, system_type=system_type, algorithm=algorithm,
    )
    feature_result = compute_feature_matrix(
        combined, candidate_features.name,
        system_type=candidate_features.system_type,
        algorithm=algorithm,
    )
    if diagnostics is not None:
        diagnostics.update({
            "feature_used": feature_result.name,
            "feature_fallbacks": list(feature_result.fallbacks),
            "system_type": feature_result.system_type,
        })
    standardized = standardize_features(feature_result.values)
    distance_result = compute_distance_matrix(
        combined, distance_metric, feature_values=standardized,
        options=distance_options, allow_fallbacks=True,
    )
    if diagnostics is not None:
        diagnostics.update(distance_result.to_dict())
    values = distance_result.values[len(representatives):, :len(representatives)]
    if not np.isfinite(values).all():
        return None
    novelties = np.min(values, axis=1)
    return [
        (float(novelty), float(molecule.energy), molecule)
        for novelty, molecule in zip(novelties, molecules)
    ]


def _basin_novelty_scores(molecules, basin_entries, *, feature="auto", distance_metric="euclidean",
                          system_type="auto", algorithm="auto", distance_options=None,
                          diagnostics=None):
    """Rank candidates by distance to stored geometries under the current policy.

    Legacy vectors without geometry are intentionally ignored: they cannot be
    safely mapped to a known descriptor implementation or coordinate witness.
    """
    try:
        scores = _geometry_novelty_scores(
            molecules, basin_entries, feature=feature, distance_metric=distance_metric,
            system_type=system_type, algorithm=algorithm, distance_options=distance_options,
            diagnostics=diagnostics,
        )
    except Exception as exc:
        scores = None
        if diagnostics is not None:
            diagnostics["comparison_failure"] = f"{type(exc).__name__}: {exc}"
    return scores


def _apply_basin_memory(molecules, maximum_number_of_seeds, basin_entries, *, feature="auto",
                        distance_metric="euclidean", system_type="auto", algorithm="auto",
                        distance_options=None, diagnostics=None):
    """Reduce a large candidate pool using geometry-backed basin novelty.

    If no stored geometry exists or a complete comparable distance matrix
    cannot be built, leave the candidate pool untouched (the in-doubt-keep
    policy). Opaque legacy Coulomb fingerprints are never used to prune.
    """
    if diagnostics is not None:
        diagnostics.update({
            "schema_version": REGISTRY_VERSION,
            "input_candidates": len(molecules),
            "registry_entries": len(basin_entries),
            "geometry_representatives": sum(_valid_geometry(entry.get("geometry")) for entry in basin_entries),
            "opaque_legacy_descriptors": sum(bool(entry.get("legacy_descriptors")) for entry in basin_entries),
            "feature_requested": feature,
            "distance_requested": distance_metric,
            "distance_options": dict(distance_options or {}),
            "applied": False,
        })
    if not basin_entries or len(molecules) <= maximum_number_of_seeds:
        if diagnostics is not None:
            diagnostics["status"] = "no-registry-entries" if not basin_entries else "within-seed-budget"
        return molecules
    representatives = [entry for entry in basin_entries if _valid_geometry(entry.get("geometry"))]
    if not representatives:
        if diagnostics is not None:
            diagnostics["status"] = "legacy-only-no-geometry"
        return molecules
    pool_cap = min(len(molecules), maximum_number_of_seeds * 2)
    if len(molecules) <= pool_cap:
        if diagnostics is not None:
            diagnostics["status"] = "within-memory-pool-cap"
        return molecules
    scored = _basin_novelty_scores(
        molecules, representatives, feature=feature, distance_metric=distance_metric,
        system_type=system_type, algorithm=algorithm, distance_options=distance_options,
        diagnostics=diagnostics,
    )
    if scored is None:
        from pyar.selection import clustering

        clustering.cluster_logger.warning(
            "Basin geometry comparison failed; keeping all %d candidates without memory pruning.",
            len(molecules),
        )
        if diagnostics is not None:
            diagnostics["status"] = "comparison-unavailable-keep-all"
        return molecules
    novelty_ranked = sorted(scored, key=lambda item: (-item[0], item[1], item[2].name))
    energy_ranked = sorted(scored, key=lambda item: (item[1], item[2].name))
    selected = []
    selected_ids = set()
    for _, _, molecule in energy_ranked[:maximum_number_of_seeds]:
        selected.append(molecule)
        selected_ids.add(id(molecule))
    for _, _, molecule in novelty_ranked:
        if len(selected) >= pool_cap:
            break
        if id(molecule) not in selected_ids:
            selected.append(molecule)
            selected_ids.add(id(molecule))
    from pyar.selection import clustering

    clustering.cluster_logger.info(
        "Geometry-backed basin memory reduced candidate pool from %d to %d before clustering.",
        len(molecules), len(selected),
    )
    if diagnostics is not None:
        diagnostics.update({"status": "applied", "output_candidates": len(selected), "applied": True})
    return selected


def record_selected_basins(selected_molecules, output_root="selected"):
    """Persist final selected geometries in the per-stoichiometry registry."""
    if not selected_molecules:
        return None
    registry_path = _basin_registry_path(selected_molecules[0], output_root=output_root)
    _persist_basin_registry(registry_path, selected_molecules)
    return registry_path


def main(argv=None):
    """Explicitly migrate basin registry files without deleting legacy data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("registry", nargs="+", help="basin_registry.json file(s) to migrate")
    parser.add_argument("--dry-run", action="store_true", help="report conversion without writing")
    args = parser.parse_args(argv)
    try:
        reports = [migrate_basin_registry(path, dry_run=args.dry_run) for path in args.registry]
    except BasinMemoryError as exc:
        parser.error(str(exc))
    for report in reports:
        print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
