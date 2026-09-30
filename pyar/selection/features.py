"""Pool-consistent structural features for selection and clustering.

Feature construction is separate from cluster assignment. Descriptor rows in a
pool always share one species vocabulary and configuration; failures are
returned as diagnostics so callers can try an explicit fallback.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

FEATURES = ("mbtr", "soap", "distance-histogram")
FEATURE_FALLBACKS = {
    "mbtr": ("soap", "distance-histogram"),
    "soap": ("mbtr", "distance-histogram"),
    "distance-histogram": (),
}


class FeatureComputationError(RuntimeError):
    """Raised when no configured structural feature can be computed."""


@dataclass(frozen=True)
class FeatureMatrix:
    """A two-dimensional, finite feature matrix and its computation record."""

    name: str
    values: np.ndarray
    species: tuple[str, ...]
    fallbacks: tuple[dict[str, str], ...] = ()


def _species_for_pool(molecules):
    return tuple(sorted({str(atom).capitalize() for molecule in molecules for atom in molecule.atoms_list}))


def _distance_histogram(molecules, species, *, bins=64, maximum_distance=12.0):
    """Build a smooth-cutoff pair-distance histogram invariant to rigid motions.

    This dependency-free representation is a last-resort clustering feature.
    It preserves element-pair counts and pair-distance distributions, but it is
    not a unique molecular structure identifier.
    """
    species_pairs = tuple(
        (species[i], species[j])
        for i in range(len(species))
        for j in range(i, len(species))
    )
    pair_lookup = {pair: index for index, pair in enumerate(species_pairs)}
    output = []
    for molecule in molecules:
        atoms = [str(atom).capitalize() for atom in molecule.atoms_list]
        coordinates = np.asarray(molecule.coordinates, dtype=float)
        if not atoms or coordinates.shape != (len(atoms), 3) or not np.isfinite(coordinates).all():
            raise ValueError(f"Invalid coordinates for {getattr(molecule, 'name', '<unnamed>')}")
        composition = np.array([atoms.count(element) for element in species], dtype=float)
        if len(atoms):
            composition /= len(atoms)
        histograms = np.zeros((len(species_pairs), bins), dtype=float)
        pair_counts = np.zeros(len(species_pairs), dtype=float)
        for left in range(len(atoms)):
            for right in range(left + 1, len(atoms)):
                pair = tuple(sorted((atoms[left], atoms[right]), key=species.index))
                pair_index = pair_lookup[pair]
                distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
                if not np.isfinite(distance):
                    raise ValueError(f"Non-finite pair distance in {getattr(molecule, 'name', '<unnamed>')}")
                # Saturate at the final bin so a separated aggregate still
                # contributes information instead of disappearing entirely.
                bin_index = min(int(distance / maximum_distance * bins), bins - 1)
                histograms[pair_index, bin_index] += 1.0
                pair_counts[pair_index] += 1.0
        for index, count in enumerate(pair_counts):
            if count:
                histograms[index] /= count
        output.append(np.concatenate((composition, histograms.reshape(-1))))
    return np.asarray(output, dtype=float)


def _compute_feature_matrix(molecules, feature, species):
    if feature == "distance-histogram":
        return _distance_histogram(molecules, species)

    from pyar import representations

    descriptors = []
    for molecule in molecules:
        if feature == "mbtr":
            values = representations.mbtr_descriptor(
                molecule.atoms_list,
                molecule.coordinates,
                species=species,
                include_angles=True,
            )
        elif feature == "soap":
            values = representations.soap_structure_descriptor(
                molecule.atoms_list,
                molecule.coordinates,
                species=species,
            )
        else:
            raise ValueError(f"Unknown structural feature: {feature!r}")
        descriptors.append(np.asarray(values, dtype=float).reshape(-1))

    lengths = {len(row) for row in descriptors}
    if len(lengths) != 1:
        raise ValueError(f"Feature {feature!r} produced inconsistent vector lengths: {sorted(lengths)}")
    return np.vstack(descriptors)


def compute_feature_matrix(molecules, feature="mbtr", *, allow_fallbacks=True):
    """Compute one finite descriptor matrix, trying recorded fallbacks on error.

    The feature order is deterministic. The distance histogram is always the
    final fallback and requires only NumPy. No candidate is removed when a
    feature fails.
    """
    molecules = list(molecules)
    if not molecules:
        return FeatureMatrix(str(feature), np.empty((0, 0)), ())
    requested = str(feature).strip().lower().replace("_", "-")
    aliases = {"histogram": "distance-histogram", "pair-distance": "distance-histogram"}
    requested = aliases.get(requested, requested)
    if requested not in FEATURES:
        raise ValueError(f"Unknown feature {feature!r}. Choose one of: {', '.join(FEATURES)}")

    species = _species_for_pool(molecules)
    attempts = [requested]
    if allow_fallbacks:
        attempts.extend(FEATURE_FALLBACKS[requested])
    failures = []
    for candidate in attempts:
        try:
            matrix = _compute_feature_matrix(molecules, candidate, species)
            if matrix.ndim != 2 or matrix.shape[0] != len(molecules) or matrix.shape[1] == 0:
                raise ValueError(f"Feature returned invalid matrix shape {matrix.shape}")
            if not np.isfinite(matrix).all():
                raise ValueError("Feature contains NaN or infinite values")
            return FeatureMatrix(candidate, matrix, species, tuple(failures))
        except Exception as exc:
            failures.append({"feature": candidate, "reason": f"{type(exc).__name__}: {exc}"})
    raise FeatureComputationError(
        "Could not compute a discriminating feature matrix; "
        + "; ".join(f"{item['feature']}: {item['reason']}" for item in failures)
    )


def standardize_features(values):
    """Scale non-constant feature columns without requiring scikit-learn."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("Feature values must be a finite two-dimensional matrix")
    center = values.mean(axis=0)
    scale = values.std(axis=0)
    active = scale > 1e-14
    if not np.any(active):
        # A constant descriptor means the structures are indistinguishable in
        # this representation. Preserve that result and let the clusterer
        # return a single group; do not invent diversity from energy or order.
        return np.zeros((len(values), 1), dtype=float)
    return (values[:, active] - center[active]) / scale[active]


__all__ = [
    "FEATURES",
    "FeatureComputationError",
    "FeatureMatrix",
    "compute_feature_matrix",
    "standardize_features",
]
