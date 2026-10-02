"""Complete structural distance matrices with explicit, pool-wide fallbacks."""

from dataclasses import dataclass, field

import numpy as np
from scipy.special import logsumexp

from pyar.selection.distances import (
    DISTANCE_METRICS, FEATURE_DISTANCE_METRICS, pairwise_distances, validate_distance_matrix,
)


@dataclass(frozen=True)
class DistanceMatrix:
    """One distance scale shared by clustering and cluster-minima trimming."""

    values: np.ndarray
    requested: str
    used: str
    units: str
    parameters: dict = field(default_factory=dict)
    fallbacks: tuple[dict, ...] = ()

    def to_dict(self):
        return {"distance_requested": self.requested, "distance_used": self.used,
                "distance_units": self.units, "distance_parameters": self.parameters,
                "distance_fallbacks": list(self.fallbacks)}


def validate_distance_options(options=None):
    """Validate scientific parameters before any operational fallback."""
    options = dict(options or {})
    unknown = set(options) - {"atom_mode", "bond_scale", "max_mappings", "irmsd_timeout",
                              "soap_cutoff", "rematch_alpha", "rematch_iterations"}
    if unknown:
        raise ValueError(f"Unknown distance options: {', '.join(sorted(unknown))}")
    if "atom_mode" in options and options["atom_mode"] not in {"heavy", "all"}:
        raise ValueError("Distance atom_mode must be 'heavy' or 'all'")
    for name in ("bond_scale", "irmsd_timeout", "soap_cutoff", "rematch_alpha"):
        if name in options and (not np.isfinite(options[name]) or options[name] <= 0):
            raise ValueError(f"{name} must be finite and positive")
    for name in ("max_mappings", "rematch_iterations"):
        if name in options and (isinstance(options[name], bool)
                               or not isinstance(options[name], int) or options[name] < 1):
            raise ValueError(f"{name} must be a positive integer")
    return options


def _rematch_similarity(local_kernel, *, alpha=1.0, max_iterations=1000, tolerance=1e-8):
    """Bounded log-domain Sinkhorn balancing for entropy-regularized matching.

    This is the REMatch transport-weighted local similarity, with uniform atom
    marginals. Log-space updates avoid underflow at small positive alpha.
    Nonconvergence raises so the caller can fall back instead of hanging.
    """
    kernel = np.asarray(local_kernel, dtype=float)
    if kernel.ndim != 2 or min(kernel.shape) < 1 or not np.isfinite(kernel).all():
        raise ValueError("REMatch requires finite nonempty local similarities")
    log_kernel = (kernel - 1) / alpha
    n, m = kernel.shape
    log_u = np.full(n, -np.log(n))
    log_v = np.full(m, -np.log(m))
    for _ in range(max_iterations):
        log_v = -np.log(m) - logsumexp(log_kernel + log_u[:, None], axis=0)
        log_u = -np.log(n) - logsumexp(log_kernel + log_v[None, :], axis=1)
        weights = np.exp(log_kernel + log_u[:, None] + log_v[None, :])
        residual = max(np.max(np.abs(weights.sum(axis=0) - 1 / m)),
                       np.max(np.abs(weights.sum(axis=1) - 1 / n)))
        if residual < tolerance:
            return float(np.sum(weights * kernel))
    raise RuntimeError("REMatch transport did not converge within its iteration limit")


def _soap_rematch_matrix(molecules, options):
    from ase import Atoms
    from dscribe.descriptors import SOAP

    species = sorted({str(atom).capitalize() for molecule in molecules for atom in molecule.atoms_list})
    parameters = {"species": species, "soap_cutoff": options.get("soap_cutoff", 5.0),
                  "n_max": 6, "l_max": 4, "sigma": .5,
                  "rematch_alpha": options.get("rematch_alpha", 1.0),
                  "rematch_iterations": options.get("rematch_iterations", 1000),
                  "transport_tolerance": 1e-8, "average": "off"}
    descriptor = SOAP(species=species, periodic=False, r_cut=parameters["soap_cutoff"],
                      n_max=6, l_max=4, sigma=.5, average="off", sparse=False)
    if sum(len(molecule.atoms_list) for molecule in molecules) * descriptor.get_number_of_features() * 8 > 256 * 1024**2:
        raise MemoryError("Local SOAP descriptor pool exceeds the 256 MiB budget")
    environments = []
    for molecule in molecules:
        values = np.asarray(descriptor.create(Atoms(symbols=molecule.atoms_list,
                                                   positions=molecule.coordinates)), dtype=float)
        norms = np.linalg.norm(values, axis=1)
        if not np.isfinite(values).all() or np.any(norms <= 0):
            raise ValueError("SOAP produced nonfinite or zero local environments")
        environments.append(values / norms[:, None])

    def similarity(left, right):
        return _rematch_similarity(left @ right.T, alpha=parameters["rematch_alpha"],
                                   max_iterations=parameters["rematch_iterations"])

    self_similarity = [similarity(values, values) for values in environments]
    distances = np.zeros((len(molecules), len(molecules)))
    for left in range(len(molecules)):
        for right in range(left + 1, len(molecules)):
            kernel = similarity(environments[left], environments[right])
            kernel /= np.sqrt(self_similarity[left] * self_similarity[right])
            if not np.isfinite(kernel) or kernel > 1 + 1e-8 or kernel < -1e-8:
                raise ValueError("Normalized REMatch similarity is outside [0, 1]")
            distance = np.sqrt(max(0., 2 * (1 - min(1., kernel))))
            distances[left, right] = distances[right, left] = distance
    parameters["distance_formula"] = "sqrt(2*(1-normalized_kernel))"
    return distances, parameters


def _mapped_matrix(molecules, metric, options):
    if metric == "graph-rmsd":
        from pyar.structure_comparison import GraphFirstDeduplicationComparator

        parameters = {"atom_mode": options.get("atom_mode", "heavy"),
                      "bond_scale": options.get("bond_scale", 1.15),
                      "max_isomorphisms": options.get("max_mappings", 10000),
                      "irmsd_timeout": options.get("irmsd_timeout", 30.)}
        comparator = GraphFirstDeduplicationComparator(**parameters)
    else:
        from pyar.structure_comparison.fragment_rmsd import FragmentRMSDComparator

        parameters = {"atom_mode": options.get("atom_mode", "all"),
                      "bond_scale": options.get("bond_scale", 1.15),
                      "max_mappings": options.get("max_mappings", 10000)}
        comparator = FragmentRMSDComparator(**parameters)
        parameters.update(max_fragment_permutations=720, distance_is_upper_bound=True)
    distances = np.zeros((len(molecules), len(molecules)))
    fallbacks = 0
    for left in range(len(molecules)):
        for right in range(left + 1, len(molecules)):
            results = [comparator.compare(molecules[left], molecules[right])]
            # Exact graph enumeration is symmetric. Approximate fragment fits
            # and native iRMSD witnesses need both directions, conservatively.
            if metric == "fragment-rmsd":
                results.append(comparator.compare(molecules[right], molecules[left]))
            for result in results:
                if not result.compatible or result.distance is None:
                    reason = result.metadata.get("reason", result.metadata.get("fallback_status", "incompatible or incomplete mapping"))
                    raise ValueError(f"Pair ({left}, {right}) has no {metric} distance: {reason}")
                if result.metadata.get("fallback_method"):
                    fallbacks += 1
            distances[left, right] = distances[right, left] = max(result.distance for result in results)
    parameters["verified_irmsd_pairs"] = fallbacks
    if metric == "graph-rmsd":
        parameters["distance_is_upper_bound"] = bool(fallbacks)
    return distances, parameters


def compute_distance_matrix(molecules, metric="euclidean", *, feature_values=None,
                            allow_fallbacks=True, options=None):
    """Build a complete matrix using one backend and scale for the whole pool.

    Structural failures try SOAP/REMatch, then standardized descriptor
    Euclidean distance. Fragment matching additionally tries exact graph RMSD.
    No missing or incompatible pair is replaced by zero or a large sentinel.
    """
    molecules = list(molecules)
    metric = str(metric).strip().lower()
    if metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {metric!r}")
    options = validate_distance_options(options)
    if metric == "fragment-rmsd":
        options.setdefault("atom_mode", "all")
    for molecule in molecules:
        coords = np.asarray(molecule.coordinates, dtype=float)
        if not molecule.atoms_list or coords.shape != (len(molecule.atoms_list), 3) or not np.isfinite(coords).all():
            raise ValueError("Structural distances require finite, nonempty geometries")
    attempts = [metric]
    if allow_fallbacks and metric not in FEATURE_DISTANCE_METRICS:
        if metric == "fragment-rmsd":
            attempts.append("graph-rmsd")
        if metric != "soap-rematch":
            attempts.append("soap-rematch")
        attempts.append("euclidean")
    failures = []
    for candidate in attempts:
        try:
            if candidate in FEATURE_DISTANCE_METRICS:
                if feature_values is None:
                    raise ValueError("Descriptor distances require feature_values")
                values = pairwise_distances(feature_values() if callable(feature_values) else feature_values, candidate)
                parameters, units = {}, "standardized-feature-units"
                if candidate == "cosine":
                    units = "dimensionless"
            elif candidate == "soap-rematch":
                values, parameters = _soap_rematch_matrix(molecules, options)
                units = "dimensionless"
            else:
                values, parameters = _mapped_matrix(molecules, candidate, options)
                units = "angstrom"
            values = validate_distance_matrix(values, len(molecules))
            if candidate not in FEATURE_DISTANCE_METRICS:
                tolerance = 1e-10 if units == "angstrom" else 1e-7
                parameters["numerical_zero_tolerance"] = tolerance
                values = values.copy()
                values[values < tolerance] = 0.
            return DistanceMatrix(values, metric, candidate, units, parameters, tuple(failures))
        except Exception as exc:
            failures.append({"distance": candidate, "reason": f"{type(exc).__name__}: {exc}"})
    raise RuntimeError("No complete distance matrix could be computed: " + "; ".join(item["reason"] for item in failures))


__all__ = ["DistanceMatrix", "compute_distance_matrix", "validate_distance_options"]
