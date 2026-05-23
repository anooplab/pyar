import itertools
import json
import logging
import math
import operator
import os
from collections import Counter
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import (
    DBSCAN,
    OPTICS,
    AffinityPropagation,
    AgglomerativeClustering,
    KMeans,
    MeanShift,
    SpectralClustering,
    estimate_bandwidth,
)
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import rbf_kernel
import hdbscan
import pyar.property
import pyar.representations

cluster_logger = logging.getLogger('pyar.cluster')
_MBTR_RUNTIME_DISABLED = False
_MBTR_DISABLE_REASON = None


def _log_seed_shortfall(requested, available, context):
    """Log that fewer unique geometries exist than the requested seed count."""
    cluster_logger.info(
        "Requested %d seeds, but only %d unique geometries were available after %s.",
        requested,
        available,
        context,
    )


def _stoichiometry_label(molecule):
    """Return a compact stoichiometry label for basin snapshots."""
    counts = Counter(molecule.atoms_list)
    parts = []

    if 'C' in counts:
        carbon = counts.pop('C')
        parts.append('C' if carbon == 1 else f'C{carbon}')
    if 'H' in counts:
        hydrogen = counts.pop('H')
        parts.append('H' if hydrogen == 1 else f'H{hydrogen}')

    for element in sorted(counts):
        count = counts[element]
        parts.append(element if count == 1 else f'{element}{count}')

    return ''.join(parts) if parts else 'unknown'


def _basin_registry_path(molecule, output_root='selected', group_by_stoichiometry=True):
    """Return the persistence path for basin memory."""
    if group_by_stoichiometry:
        return os.path.join(
            output_root,
            f'stoichiometry_{_stoichiometry_label(molecule)}',
            'basin_registry.json',
        )
    return os.path.join(output_root, 'basin_registry.json')


def _fingerprint_signature(molecule):
    """Return a stable, normalized fingerprint signature for a molecule."""
    coordinates = np.asarray(molecule.coordinates, dtype=float)
    signature = pyar.representations.fingerprint(molecule.atoms_list, coordinates)
    signature = np.real_if_close(signature)
    signature = np.asarray(signature, dtype=float).ravel()
    norm = np.linalg.norm(signature)
    if norm > 0:
        signature = signature / norm
    return signature.tolist()


def _entry_fingerprint(entry):
    """Return a registry fingerprint as an array, or None if it is invalid."""
    try:
        return np.asarray(entry['fingerprint'], dtype=float).ravel()
    except (KeyError, TypeError, ValueError):
        return None


def _load_basin_registry(registry_path):
    """Load previously selected basin signatures from disk."""
    if not registry_path or not os.path.exists(registry_path):
        return []

    try:
        with open(registry_path, 'r') as fp:
            payload = json.load(fp)
    except Exception as exc:
        cluster_logger.warning("Could not load basin registry %s: %s", registry_path, exc)
        return []

    entries = payload.get('entries', []) if isinstance(payload, dict) else payload
    return [entry for entry in entries if isinstance(entry, dict) and 'fingerprint' in entry]


def _persist_basin_registry(registry_path, selected_molecules, existing_entries=None, max_entries=200):
    """Persist selected basins so later runs can avoid rediscovering them."""
    if not registry_path or not selected_molecules:
        return

    if existing_entries is None:
        existing_entries = _load_basin_registry(registry_path)

    seen_signatures = {
        tuple(np.round(fingerprint, 6))
        for fingerprint in (_entry_fingerprint(entry) for entry in existing_entries)
        if fingerprint is not None
    }
    updated_entries = list(existing_entries)
    for molecule in selected_molecules:
        signature = _fingerprint_signature(molecule)
        signature_key = tuple(np.round(np.asarray(signature, dtype=float), 6))
        if signature_key in seen_signatures:
            continue
        seen_signatures.add(signature_key)
        updated_entries.append(
            {
                'name': molecule.name,
                'energy': float(molecule.energy),
                'fingerprint': signature,
            }
        )

    if not updated_entries:
        return

    target_directory = os.path.dirname(registry_path)
    if target_directory and not os.path.exists(target_directory):
        os.makedirs(target_directory, exist_ok=True)

    payload = {
        'stoichiometry': os.path.basename(os.path.dirname(registry_path)).replace('stoichiometry_', ''),
        'entries': updated_entries[-max_entries:],
    }
    with open(registry_path, 'w') as fp:
        json.dump(payload, fp, indent=2)

    cluster_logger.info(
        "Updated basin registry with %d selected geometries at %s",
        len(selected_molecules),
        registry_path,
    )


def _basin_novelty_scores(molecules, basin_entries):
    """Rank molecules by novelty against previous basin signatures."""
    basin_signatures = [
        fingerprint
        for fingerprint in (_entry_fingerprint(entry) for entry in basin_entries)
        if fingerprint is not None
    ]
    scores = []
    for molecule in molecules:
        signature = np.asarray(_fingerprint_signature(molecule), dtype=float)
        comparable_signatures = [
            basin for basin in basin_signatures if basin.shape == signature.shape
        ]
        if comparable_signatures:
            novelty = min(np.linalg.norm(signature - basin) for basin in comparable_signatures)
        else:
            novelty = np.inf
        scores.append((novelty, float(molecule.energy), molecule))
    return scores


def _apply_basin_memory(molecules, maximum_number_of_seeds, basin_entries):
    """Reduce the MBTR candidate pool using previous basin discoveries."""
    if not basin_entries or len(molecules) <= maximum_number_of_seeds:
        return molecules

    pool_cap = min(len(molecules), maximum_number_of_seeds * 2)
    if len(molecules) <= pool_cap:
        return molecules

    scored_molecules = _basin_novelty_scores(molecules, basin_entries)
    novelty_ranked = sorted(scored_molecules, key=lambda item: (-item[0], item[1], item[2].name))
    energy_ranked = sorted(scored_molecules, key=lambda item: (item[1], item[2].name))

    selected = []
    selected_ids = set()
    for _, _, molecule in energy_ranked[:maximum_number_of_seeds]:
        selected.append(molecule)
        selected_ids.add(id(molecule))
    for _, _, molecule in novelty_ranked:
        if len(selected) >= pool_cap:
            break
        if id(molecule) in selected_ids:
            continue
        selected.append(molecule)
        selected_ids.add(id(molecule))

    cluster_logger.info(
        "Basin memory reduced candidate pool from %d to %d before MBTR selection.",
        len(molecules),
        len(selected),
    )
    return selected


def _prefer_connected_structures(molecules):
    """Prefer geometries that remain connected under a covalent-radius graph."""
    if len(molecules) < 2:
        return molecules

    try:
        from pyar import trial_generation
    except Exception as exc:
        cluster_logger.debug(
            "Connectivity filter unavailable, keeping original candidate pool: %s",
            exc,
        )
        return molecules

    connected = []
    disconnected = []
    for molecule in molecules:
        try:
            if trial_generation.broken(molecule):
                disconnected.append(molecule)
            else:
                connected.append(molecule)
        except Exception as exc:
            cluster_logger.debug(
                "Connectivity check failed for %s, keeping it in the connected pool: %s",
                getattr(molecule, 'name', 'unknown'),
                exc,
            )
            connected.append(molecule)

    if not connected:
        cluster_logger.warning(
            "All candidate geometries are disconnected; continuing with the original pool."
        )
        return molecules

    if disconnected:
        cluster_logger.info(
            "Discarded %d disconnected candidate geometries before clustering.",
            len(disconnected),
        )

    return connected


def _finalize_selection(selected_molecules, basin_registry_path, existing_entries=None):
    """Persist basin memory and return the selected molecules."""
    if basin_registry_path:
        _persist_basin_registry(basin_registry_path, selected_molecules, existing_entries=existing_entries)
    return selected_molecules


def record_selected_basins(selected_molecules, output_root='selected'):
    """Persist final selected geometries in the per-stoichiometry basin registry."""
    if not selected_molecules:
        return None
    registry_path = _basin_registry_path(selected_molecules[0], output_root=output_root)
    _persist_basin_registry(registry_path, selected_molecules)
    return registry_path


def _limit_seed_count(molecules, maximum_number_of_seeds, reason):
    """Enforce a hard upper bound on selected molecules."""
    if len(molecules) <= maximum_number_of_seeds:
        if len(molecules) < maximum_number_of_seeds:
            _log_seed_shortfall(maximum_number_of_seeds, len(molecules), reason)
        return molecules
    cluster_logger.warning(
        "Seed count exceeded limit after %s: keeping lowest-energy %d of %d",
        reason,
        maximum_number_of_seeds,
        len(molecules),
    )
    limited = sorted(molecules, key=lambda m: float(m.energy))[:maximum_number_of_seeds]
    print_energy_table(limited)
    return limited


def _max_min_diversity_select(features, molecules, maximum_number_of_seeds, initial_selected_indices=None):
    """Select diverse seeds via greedy max-min (farthest-point) selection.

    When initial_selected_indices is provided, those points are treated as fixed
    anchors and only the additional selected molecules are returned.
    """
    if initial_selected_indices is None:
        initial_selected_indices = []

    if len(molecules) <= maximum_number_of_seeds and not initial_selected_indices:
        return molecules

    selected_indices = list(initial_selected_indices)
    returned_indices = []
    remaining_indices = set(range(len(molecules))) - set(selected_indices)

    if not selected_indices:
        energies = np.array([float(m.energy) for m in molecules])
        start_index = int(np.argmin(energies))
        selected_indices.append(start_index)
        returned_indices.append(start_index)
        remaining_indices.remove(start_index)

    while len(selected_indices) < maximum_number_of_seeds and remaining_indices:
        best_idx = None
        best_score = -np.inf
        for idx in remaining_indices:
            min_distance = min(
                np.linalg.norm(features[idx] - features[sidx]) for sidx in selected_indices
            )
            if min_distance > best_score and not np.isclose(min_distance, best_score):
                best_score = min_distance
                best_idx = idx
            elif np.isclose(min_distance, best_score) and best_idx is not None:
                if float(molecules[idx].energy) < float(molecules[best_idx].energy):
                    best_idx = idx
        selected_indices.append(best_idx)
        returned_indices.append(best_idx)
        remaining_indices.remove(best_idx)

    selected = [molecules[i] for i in returned_indices]
    if initial_selected_indices:
        cluster_logger.info(
            "Added %d filler seeds by max-min strategy.",
            len(selected),
        )
    else:
        cluster_logger.info(
            "Selected %d diverse seeds by max-min strategy.",
            len(selected),
        )
    print_energy_table(selected)
    return selected


def remove_similar(list_of_molecules):
    ordered_molecules = sorted(list_of_molecules, key=lambda molecule: (float(molecule.energy), molecule.name))
    final_list = []
    rmsd_threshold = _adaptive_duplicate_rmsd_threshold(ordered_molecules)
    cluster_logger.debug('Number of molecules before similarity elimination,  {}'.format(len(ordered_molecules)))
    for candidate in ordered_molecules:
        duplicate = False
        for kept in final_list:
            if len(candidate.atoms_list) < 2 or len(kept.atoms_list) < 2:
                continue
            if not _structure_is_similar(candidate, kept):
                continue
            if _rmsd_after_alignment(candidate, kept) < rmsd_threshold:
                duplicate = True
                cluster_logger.debug(
                    'Removing {} as a near-duplicate of {}'.format(candidate.name, kept.name)
                )
                break
        if not duplicate:
            final_list.append(candidate)
    cluster_logger.debug('Number of molecules after similarity elimination,  {}'.format(len(final_list)))
    print_energy_table(final_list)
    return final_list


def _adaptive_duplicate_rmsd_threshold(molecules):
    """Estimate a duplicate RMSD threshold from the current candidate pool."""
    if len(molecules) < 2:
        return 0.10

    sampled_rmsd = []
    sample_limit = 40
    for left_index, left in enumerate(molecules):
        if len(sampled_rmsd) >= sample_limit:
            break
        for right in molecules[left_index + 1:]:
            if len(sampled_rmsd) >= sample_limit:
                break
            if len(left.atoms_list) != len(right.atoms_list):
                continue
            if not _structure_is_similar(left, right):
                continue
            try:
                sampled_rmsd.append(_rmsd_after_alignment(left, right))
            except Exception:
                continue

    if not sampled_rmsd:
        return 0.10

    finite_rmsd = [value for value in sampled_rmsd if np.isfinite(value)]
    if not finite_rmsd:
        return 0.10

    baseline = float(np.percentile(finite_rmsd, 5))
    threshold = max(0.05, baseline * 0.75)
    return float(np.clip(threshold, 0.05, 0.15))


def _structure_is_similar(candidate, kept):
    """Quick prefilter using fingerprint distance."""
    fingerprint_distance = calc_fingerprint_distance(candidate, kept)
    return abs(fingerprint_distance) < 1.0


def _kabsch_rotation(reference, mobile):
    """Return the proper rotation matrix aligning ``mobile`` onto ``reference``."""
    covariance = mobile.T @ reference
    left_vectors, _, right_vectors = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(left_vectors @ right_vectors) < 0:
        correction[-1, -1] = -1.0
    return left_vectors @ correction @ right_vectors


def _kabsch_rmsd(reference, mobile):
    """Return RMSD after optimal translation and proper rotation."""
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    rotation = _kabsch_rotation(reference_centered, mobile_centered)
    aligned = mobile_centered @ rotation
    diff = aligned - reference_centered
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _equivalent_atom_groups(atoms):
    """Return indices grouped by element, in deterministic element order."""
    return {
        element: [index for index, atom in enumerate(atoms) if atom == element]
        for element in sorted(set(atoms))
    }


def _exact_element_orders(reference_atoms, mobile_atoms, maximum_permutations=720):
    """Yield element-preserving mobile atom orders when exact matching is tractable."""
    reference_groups = _equivalent_atom_groups(reference_atoms)
    mobile_groups = _equivalent_atom_groups(mobile_atoms)
    permutation_count = math.prod(
        math.factorial(len(indices)) for indices in mobile_groups.values()
    )
    if permutation_count > maximum_permutations:
        return None

    elements = list(reference_groups)
    permutations = [
        list(itertools.permutations(mobile_groups[element])) for element in elements
    ]

    def orders():
        for element_orders in itertools.product(*permutations):
            mobile_order = [None] * len(reference_atoms)
            for element, element_order in zip(elements, element_orders):
                for reference_index, mobile_index in zip(
                    reference_groups[element],
                    element_order,
                ):
                    mobile_order[reference_index] = mobile_index
            yield mobile_order

    return orders()


def _assigned_element_order(reference_atoms, reference, mobile_atoms, mobile):
    """Match mobile atoms to reference atoms with element-constrained assignment."""
    order = [None] * len(reference_atoms)
    reference_groups = _equivalent_atom_groups(reference_atoms)
    mobile_groups = _equivalent_atom_groups(mobile_atoms)
    for element, reference_indices in reference_groups.items():
        mobile_indices = mobile_groups[element]
        distance_matrix = np.sum(
            (
                reference[np.asarray(reference_indices)][:, None, :]
                - mobile[np.asarray(mobile_indices)][None, :, :]
            )
            ** 2,
            axis=2,
        )
        row_indices, column_indices = linear_sum_assignment(distance_matrix)
        for row_index, column_index in zip(row_indices, column_indices):
            order[reference_indices[row_index]] = mobile_indices[column_index]
    return order


def _iterative_assigned_rmsd(reference_atoms, reference, mobile_atoms, mobile):
    """Estimate permutation-aware RMSD for systems too large for exact matching."""
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    order = _assigned_element_order(
        reference_atoms,
        reference_centered,
        mobile_atoms,
        mobile_centered,
    )
    for _ in range(8):
        rotation = _kabsch_rotation(reference_centered, mobile_centered[order])
        aligned_mobile = mobile_centered @ rotation
        next_order = _assigned_element_order(
            reference_atoms,
            reference_centered,
            mobile_atoms,
            aligned_mobile,
        )
        if next_order == order:
            break
        order = next_order
    return _kabsch_rmsd(reference, mobile[order])


def _rmsd_after_alignment(candidate, kept):
    """Return translation-, rotation-, and element-order-invariant RMSD."""
    candidate_coords = np.asarray(candidate.coordinates, dtype=float)
    kept_coords = np.asarray(kept.coordinates, dtype=float)
    if (
        candidate_coords.shape != kept_coords.shape
        or Counter(candidate.atoms_list) != Counter(kept.atoms_list)
    ):
        return float('inf')

    candidate_orders = _exact_element_orders(kept.atoms_list, candidate.atoms_list)
    if candidate_orders is not None:
        return min(
            _kabsch_rmsd(kept_coords, candidate_coords[order])
            for order in candidate_orders
        )

    return _iterative_assigned_rmsd(
        kept.atoms_list,
        kept_coords,
        candidate.atoms_list,
        candidate_coords,
    )

def calc_fingerprint_distance(a, b):
    """Calculate the distance between two fingerprints"""
    return np.linalg.norm(
        pyar.representations.fingerprint(a.atoms_list, a.coordinates)
        - pyar.representations.fingerprint(b.atoms_list, b.coordinates)
    )


def choose_geometries(
    list_of_molecules,
    maximum_number_of_seeds=12,
    persist_basin_memory=True,
    apply_basin_memory=True,
    algorithm=None,
    group_basin_by_stoichiometry=True,
):
    global _MBTR_RUNTIME_DISABLED, _MBTR_DISABLE_REASON
    if len(list_of_molecules) < 2:
        _log_seed_shortfall(maximum_number_of_seeds, len(list_of_molecules), "selection")
        basin_registry_path = _basin_registry_path(
            list_of_molecules[0],
            group_by_stoichiometry=group_basin_by_stoichiometry,
        ) if list_of_molecules and os.path.isdir('selected') else None
        write_path = basin_registry_path if persist_basin_memory else None
        return _finalize_selection(list_of_molecules, write_path)

    if len(list_of_molecules) <= maximum_number_of_seeds:
        cluster_logger.info('Not enough data for clustering. Removing similar geometries from the list')
        basin_registry_path = _basin_registry_path(
            list_of_molecules[0],
            group_by_stoichiometry=group_basin_by_stoichiometry,
        ) if os.path.isdir('selected') else None
        write_path = basin_registry_path if persist_basin_memory else None
        selected = _limit_seed_count(
            remove_similar(list_of_molecules),
            maximum_number_of_seeds,
            reason="similarity pruning",
        )
        return _finalize_selection(selected, write_path)

    # Read selection algorithm from environment variable, default to the
    # hybrid cluster-plus-fill behavior when no explicit algorithm is provided.
    if algorithm is None:
        algorithm = os.environ.get('PYAR_CLUSTERING_ALGORITHM', 'hybrid')
    algorithm = algorithm.lower()
    cluster_logger.info(f'Seed selection on {len(list_of_molecules)} geometries using {algorithm}')

    pruned_molecules = remove_similar(list_of_molecules)
    basin_registry_path = _basin_registry_path(
        pruned_molecules[0],
        group_by_stoichiometry=group_basin_by_stoichiometry,
    ) if pruned_molecules and os.path.isdir('selected') else None
    write_path = basin_registry_path if persist_basin_memory else None
    basin_entries = _load_basin_registry(basin_registry_path) if apply_basin_memory else []
    if basin_entries:
        pruned_molecules = _apply_basin_memory(pruned_molecules, maximum_number_of_seeds, basin_entries)
    pruned_molecules = _prefer_connected_structures(pruned_molecules)

    if len(pruned_molecules) <= maximum_number_of_seeds:
        cluster_logger.info(
            "Similarity pruning reduced seed pool to %d; skipping diversity selection.",
            len(pruned_molecules),
        )
        return _finalize_selection(pruned_molecules, write_path, existing_entries=basin_entries)

    if _MBTR_RUNTIME_DISABLED:
        cluster_logger.info(
            "MBTR clustering disabled for this run (%s). Falling back to similarity pruning.",
            _MBTR_DISABLE_REASON,
        )
        selected = _limit_seed_count(
            pruned_molecules,
            maximum_number_of_seeds,
            reason="MBTR-disabled fallback",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    # Use MBTR for feature representation. Some DScribe/ASE combinations can
    # fail during conversion, so disable MBTR for the current run after first
    # failure and use similarity pruning thereafter.
    try:
        dt = np.array([pyar.representations.mbtr_descriptor(m.atoms_list, m.coordinates) for m in pruned_molecules])
    except Exception as exc:
        _MBTR_RUNTIME_DISABLED = True
        _MBTR_DISABLE_REASON = str(exc)
        cluster_logger.warning(
            "MBTR clustering unavailable, disabling MBTR for this run and falling back to similarity pruning: %s",
            exc,
        )
        selected = _limit_seed_count(
            pruned_molecules,
            maximum_number_of_seeds,
            reason="MBTR failure fallback",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    # Scale descriptor features before diversity/clustering selection.
    scaler = StandardScaler()
    dt_scaled = scaler.fit_transform(dt)

    if os.environ.get("PYAR_SAVE_MBTR_FEATURES") == "1":
        pd.DataFrame(dt_scaled).to_csv("mbtr_features.csv")

    if algorithm in {"maxmin", "max-min", "max_min"}:
        selected = _limit_seed_count(
            _max_min_diversity_select(dt_scaled, pruned_molecules, maximum_number_of_seeds),
            maximum_number_of_seeds,
            reason="max-min selection",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    try:
        cluster_algorithm = "hdbscan" if algorithm == "hybrid" else algorithm
        labels = generate_labels(dt_scaled, cluster_algorithm, maximum_number_of_seeds)
    except Exception as exc:
        cluster_logger.warning(
            "Clustering algorithm %s failed (%s); falling back to max-min selection.",
            algorithm,
            exc,
        )
        cluster_logger.debug("Clustering failure traceback", exc_info=True)
        selected = _limit_seed_count(
            _max_min_diversity_select(
                dt_scaled,
                pruned_molecules,
                maximum_number_of_seeds,
            ),
            maximum_number_of_seeds,
            reason=f"{algorithm} failure fallback",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    best_from_each_cluster = select_best_from_each_cluster(labels, pruned_molecules)

    if len(best_from_each_cluster) > maximum_number_of_seeds:
        cluster_logger.info(
            "Cluster selection returned %d minima; trimming to %d with max-min.",
            len(best_from_each_cluster),
            maximum_number_of_seeds,
        )
        selected_ids = {id(molecule) for molecule in best_from_each_cluster}
        cluster_feature_indices = [
            index for index, molecule in enumerate(pruned_molecules)
            if id(molecule) in selected_ids
        ]
        cluster_subset_features = dt_scaled[cluster_feature_indices]
        cluster_subset_molecules = [pruned_molecules[index] for index in cluster_feature_indices]
        trimmed = _max_min_diversity_select(
            cluster_subset_features,
            cluster_subset_molecules,
            maximum_number_of_seeds,
        )
        selected = _limit_seed_count(
            trimmed,
            maximum_number_of_seeds,
            reason="cluster minima trimming",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    if len(best_from_each_cluster) == maximum_number_of_seeds:
        selected = _limit_seed_count(
            best_from_each_cluster,
            maximum_number_of_seeds,
            reason="cluster selection",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    selected_ids = {id(m) for m in best_from_each_cluster}
    if len(selected_ids) == len(pruned_molecules):
        return _finalize_selection(best_from_each_cluster, write_path, existing_entries=basin_entries)

    cluster_logger.info(
        "Cluster selection returned %d seeds; filling remaining %d with max-min.",
        len(best_from_each_cluster),
        maximum_number_of_seeds - len(best_from_each_cluster),
    )
    selected_indices = [
        index for index, molecule in enumerate(pruned_molecules)
        if id(molecule) in selected_ids
    ]
    filler = _max_min_diversity_select(
        dt_scaled,
        pruned_molecules,
        maximum_number_of_seeds,
        initial_selected_indices=selected_indices,
    )
    combined = best_from_each_cluster + filler
    selected = _limit_seed_count(
        combined,
        maximum_number_of_seeds,
        reason="hybrid cluster-fill selection",
    )
    return _finalize_selection(selected, write_path, existing_entries=basin_entries)

def generate_labels(dt, algorithm='hdbscan', maximum_number_of_seeds=8):
    if algorithm == 'kmeans':
        return kmeans_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'dbscan':
        return dbscan_clustering(dt)
    elif algorithm == 'optics':
        return optics_clustering(dt)
    elif algorithm in {'affinity', 'affinity_propagation'}:
        return affinity_propagation_clustering(dt)
    elif algorithm in {'meanshift', 'mean_shift'}:
        return mean_shift_clustering(dt)
    elif algorithm in {'agglomerative', 'ward'}:
        return agglomerative_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'spectral':
        return spectral_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'hdbscan':
        return hdbscan_clustering(dt)
    elif algorithm == 'gaussian_mixture':
        return gaussian_mixture_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'rbf_kernel':
        return rbf_kernel_clustering(dt)
    else:
        cluster_logger.warning(f"Unknown algorithm: {algorithm}. Using HDBSCAN.")
        return hdbscan_clustering(dt)

def kmeans_clustering(dt, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    return kmeans.fit_predict(dt)

def dbscan_clustering(dt):
    eps, min_samples = determine_dbscan_params(dt)
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    return dbscan.fit_predict(dt)


def optics_clustering(dt):
    clusterer = OPTICS(min_samples=2, xi=0.05, min_cluster_size=2)
    return clusterer.fit_predict(dt)

def hdbscan_clustering(dt):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1)
    return clusterer.fit_predict(dt)


def affinity_propagation_clustering(dt):
    clusterer = AffinityPropagation(random_state=42)
    return clusterer.fit_predict(dt)


def mean_shift_clustering(dt):
    bandwidth = estimate_bandwidth(dt, quantile=0.2, n_samples=min(len(dt), 500))
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        bandwidth = None
    clusterer = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    return clusterer.fit_predict(dt)


def spectral_clustering(dt, n_clusters):
    n_clusters = max(2, min(n_clusters, len(dt)))
    n_neighbors = max(1, min(10, len(dt) - 1))
    clusterer = SpectralClustering(
        n_clusters=n_clusters,
        random_state=42,
        assign_labels='kmeans',
        affinity='nearest_neighbors',
        n_neighbors=n_neighbors,
    )
    return clusterer.fit_predict(dt)


def agglomerative_clustering(dt, n_clusters):
    n_clusters = max(2, min(n_clusters, len(dt)))
    clusterer = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
    return clusterer.fit_predict(dt)

def gaussian_mixture_clustering(dt, n_components):
    gm = GaussianMixture(n_components=n_components, random_state=42)
    return gm.fit_predict(dt)

def rbf_kernel_clustering(dt, threshold=0.99):
    similarities = rbf_kernel(dt)
    n_samples = similarities.shape[0]
    labels = np.zeros(n_samples, dtype=int)
    current_label = 0
    for i in range(n_samples):
        if labels[i] == 0:
            current_label += 1
            labels[i] = current_label
            for j in range(i+1, n_samples):
                if similarities[i, j] > threshold:
                    labels[j] = current_label
    return labels

def determine_dbscan_params(dt):
    # Simple heuristic for DBSCAN parameters
    distances = np.sort(np.sum((dt[:, None, :] - dt[None, :, :]) ** 2, axis=-1), axis=1)
    eps = np.median(distances[:, 1])
    min_samples = 2
    return eps, min_samples

# def select_best_from_each_cluster(labels, list_of_molecules):
#     unique_labels = np.unique(labels)
#     cluster_logger.info(f"The distribution of file in each cluster: {np.bincount(labels)}")
#     best_from_each_cluster = []
#     for this_label in unique_labels:
#         if this_label != -1:  # -1 is the noise label in some clustering algorithms
#             molecules_in_this_group = [m for label, m in zip(labels, list_of_molecules) if label == this_label]
#             best_from_each_cluster.append(get_the_best_molecule(molecules_in_this_group))
#     cluster_logger.info("Lowest energy structures from each cluster")
#     print_energy_table(best_from_each_cluster)
#     return best_from_each_cluster

def select_best_from_each_cluster(labels, list_of_molecules):
    labels = np.array(labels)  # Ensure labels is a numpy array
    unique_labels = np.unique(labels)
    noise_molecules = []
    
    # Handle the case where there are negative labels (noise points)
    if np.any(labels < 0):
        cluster_logger.info("Clustering algorithm identified noise points.")
        noise_molecules = [m for l, m in zip(labels, list_of_molecules) if l == -1]
        positive_labels = labels[labels >= 0]
        if len(positive_labels) > 0:
            cluster_logger.info(f"The distribution of files in each cluster (excluding noise): {np.bincount(positive_labels)}")
        else:
            cluster_logger.info("No valid clusters found.")
    else:
        cluster_logger.info(f"The distribution of files in each cluster: {np.bincount(labels)}")

    best_from_each_cluster = []
    for label in unique_labels:
        if label != -1:  # Exclude noise points (label -1)
            molecules_in_this_group = [m for l, m in zip(labels, list_of_molecules) if l == label]
            if molecules_in_this_group:
                best_from_each_cluster.append(get_the_best_molecule(molecules_in_this_group))

    if noise_molecules:
        best_noise = get_the_best_molecule(noise_molecules)
        cluster_logger.info(
            "Including best noise-point representative: %s (energy %.6f)",
            best_noise.name,
            float(best_noise.energy),
        )
        best_from_each_cluster.append(best_noise)

    cluster_logger.info("Lowest energy structures from each cluster:")
    print_energy_table(best_from_each_cluster)
    return best_from_each_cluster

def get_the_best_molecule(list_of_molecules):
    return min(list_of_molecules, key=lambda m: m.energy)

def print_energy_table(molecules, stream=None, title=None):
    """Report energies with relative values against the global minimum."""
    e_dict = {i.name: float(i.energy) for i in molecules}
    if not e_dict:
        return

    ref = min(e_dict.values())
    lines = []
    if title:
        lines.append(title)
    lines.append(f"{'Name':>35}  {'Energy':>12}  {'R. E. (kcal/mol)':>18}")
    for name, energy in sorted(e_dict.items(), key=operator.itemgetter(1)):
        lines.append(f"{name:>35}  {energy:12.6f}  {(energy - ref) * 627.51:18.2f}")
    lines.append(f"Global minimum: {min(e_dict, key=e_dict.get)} ({ref:12.6f} Eh)")
    lines.append("")

    if stream is not None:
        for line in lines:
            print(line, file=stream)
    else:
        for line in lines:
            cluster_logger.info(line)



def read_energy_from_xyz_file(xyz_file):
    import re
    try:
        with open(xyz_file, 'r') as fr:
            lines = fr.readlines()
            second_line = lines[1].strip()
            energy_matches = re.findall(r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', second_line)
            if energy_matches:
                energy = float(energy_matches[-1])
            else:
                raise ValueError("No numeric energy found")
    except (IndexError, ValueError):
        with open(xyz_file, 'r') as fr:
            comments_line = fr.readlines()[1].rstrip()
            energy = float(re.split(r':|=|\s+', comments_line)[1])

    return energy


def plot_energy_histogram(molecules):
    energies = [i.energy for i in molecules]
    ref = min(energies)
    relative_energies = [(float(energy) - ref) for energy in energies]
    histogram_bin = np.linspace(0, max(relative_energies), 10)
    import matplotlib.pyplot as plt
    plt.hist(relative_energies, histogram_bin)
    plt.xlabel('Energy')
    plt.ylabel('Population')
    plt.title('Histogram of energies')

    
def main():
    pass

if __name__ == "__main__":
    main()
