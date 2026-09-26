"""Reaction workflow orchestration for PyAR.

This module owns the AFIR reaction-search pipeline:

* validate the reaction request and any restart state
* create or resume the ``reaction/`` run directory
* build the numeric gamma schedule
* generate trial orientations for each gamma cycle
* optimize each orientation with the selected backend
* perform unbiased relaxation when a bonded candidate is found
* deduplicate and persist unique products
* persist ``reaction/state.json`` plus product and trace artifacts
* emit trace-analysis summaries for successful paths
* return a structured :class:`~pyar.workflow_results.ReactionResult`

The module is the canonical implementation behind the ``pyar-react``
command-line entry point.
"""

import logging
import json
import math
import os
import shutil
import sys

import numpy as np
import pyar.scan
from pyar import file_manager
from pyar.data import defualt_parameters
from pyar.backend_capabilities import backend_supports_native_optimization
from pyar.selection import clustering
from pyar.optimiser import is_cycle_exceeded, is_success, is_usable, optimise
from pyar.sampling import trial_generator as trial_generation
from pyar.workflow_results import ReactionResult
from pyar import reaction_analysis
from pyar.reaction_identity import (
    molecule_identity_from_xyz,
    reaction_product_changed,
    separated_reactant_identity,
    same_molecular_identity,
    write_disconnected_reference,
)
from pyar.reaction_trace import infer_bonds
from pyar.release import _inside_covalent_radius_sum
from pyar.state.reaction import ReactionRunState, ReactionStateError, read_legacy_checkpoint
from pyar.workflows._growth import (
    sampling_configuration,
    workflow_run_directory,
    workflow_state_path,
)

reactor_logger = logging.getLogger('pyar.workflows.reaction')

saved_product_identities = {}


def _build_reaction_result(
    workdir,
    status,
    product_dir,
    run_state,
    gamma_list,
    orientations,
    sampling,
):
    """Package the current reaction outcome as a structured result.

    The returned :class:`~pyar.workflow_results.ReactionResult` captures the
    final workflow state, the reaction run directory, the persisted restart
    state, and the remaining orientation set so callers can inspect or report
    the result without parsing the on-disk state tree.
    """
    return ReactionResult(
        workflow="reaction",
        status=status,
        run_directory=workflow_run_directory(workdir, "reaction"),
        state_path=workflow_state_path(workflow_run_directory(workdir, "reaction")),
        selected_paths=tuple(product["path"] for product in run_state.data.get("products", [])),
        metadata={
            "gamma_schedule": tuple(float(value) for value in gamma_list),
            "products": tuple(run_state.data.get("products", [])),
            "product_directory": product_dir,
            "remaining_orientations": len(orientations),
            "sampling": run_state.data.get("sampling", sampling),
        },
    )


def print_header(gamma_max, gamma_min, hm_orientations, software):
    """Log the fixed header used at the start of a reaction run.

    The header mirrors the historical reactor logging style and gives the
    caller a single place to confirm the gamma range, orientation count, and
    backend before the workflow begins mutating the working directory.
    """
    reactor_logger.info("==================== PyAR Reaction Workflow ====================")
    reactor_logger.info(f"Gamma range: {gamma_min} to {gamma_max}")
    reactor_logger.info(f"Orientations: {hm_orientations}")
    reactor_logger.info(f"Software: {software}")
    reactor_logger.info("===============================================================")


def with_gamma(qc_params, gamma):
    """Return a copy of ``qc_params`` with a specific AFIR gamma applied.

    The helper also enables trace recording only for geomeTRIC-backed runs
    with a non-zero gamma, because those are the only configurations that
    produce the path trace used by the analysis tooling.
    """
    updated_qc_params = dict(qc_params)
    updated_qc_params['gamma'] = gamma
    trace_enabled = (
        updated_qc_params.get("geometry_optimizer") == "geometric"
        and float(gamma) != 0.0
    )
    updated_qc_params["trace_enabled"] = trace_enabled
    updated_qc_params["reaction_trace"] = trace_enabled
    return updated_qc_params


def without_afir_bias(qc_params):
    """Return parameters for unbiased physical relaxation.

    This is the relaxation step applied after a bonded AFIR candidate has
    been identified. It preserves the physical backend configuration while
    forcing ``gamma=0.0`` and the backend-native optimizer so the candidate
    is re-optimized on the unbiased physical objective. Bias-controller and
    trace settings are deliberately removed rather than merely disabled.
    """
    unbiased = dict(qc_params)
    unbiased.setdefault("opt_cycles", defualt_parameters.values["opt_cycles"])
    unbiased.setdefault("opt_threshold", defualt_parameters.values["opt_threshold"])
    unbiased['gamma'] = 0.0
    unbiased['geometry_optimizer'] = 'native'
    unbiased['trace_enabled'] = False
    unbiased['reaction_trace'] = False
    for key in (
        'bias_controller', 'bias_controller_restart', 'bias_alpha_min',
        'bias_alpha_margin', 'bias_alpha_smoothing', 'bias_alpha_epsilon',
        'bias_scheduled_alpha', 'bias_potential', 'softmin_beta',
        'release_job_name', 'release_retry_attempt', 'release_retry_limit',
        'release_margin_factor', 'release_alpha_critical_max', 'release_persistence',
        'release_distance_fraction',
        'release_bond_order_growth',
    ):
        unbiased.pop(key, None)
    return unbiased


def build_gamma_schedule(gamma_min, gamma_max, steps=10):
    """Build the numeric AFIR gamma schedule used by the reaction workflow.

    The schedule is inclusive and monotonic. A single-valued schedule is
    returned when the limits are equal; otherwise the workflow uses evenly
    spaced values between the endpoints.
    """
    if not np.isfinite(gamma_min) or not np.isfinite(gamma_max):
        raise ValueError("AFIR gamma limits must be finite numbers")
    if gamma_min < 0.0 or gamma_max < 0.0:
        raise ValueError("AFIR gamma limits must be non-negative")
    if gamma_max < gamma_min:
        raise ValueError("AFIR maximum gamma must be greater than or equal to minimum gamma")
    if gamma_min == gamma_max:
        return np.asarray([float(gamma_min)])
    return np.linspace(gamma_min, gamma_max, num=steps, dtype=float)


def format_gamma_id(gamma):
    """Format a gamma value for directory names and readable job labels.

    The formatter keeps directory names stable and lexicographically useful
    by zero-padding the integral part and encoding the decimal separator as
    ``p``.
    """
    value = float(gamma)
    sign = "m" if value < 0.0 else ""
    magnitude = f"{abs(value):.12g}"
    if "e" not in magnitude:
        integer, separator, fraction = magnitude.partition(".")
        magnitude = integer.zfill(4)
        if separator:
            magnitude = f"{magnitude}p{fraction}"
    else:
        magnitude = magnitude.replace(".", "p").replace("-", "m").replace("+", "")
    return f"{sign}{magnitude}"


def _molecule_signature(molecule):
    """Return stable input geometry metadata used to validate restarts.

    The signature is deliberately small and deterministic so restart
    validation can reject changed reactants before the workflow mutates the
    reaction directory.
    """
    return {
        "atoms": list(molecule.atoms_list),
        "coordinates": np.asarray(molecule.coordinates, dtype=float).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
    }


def build_reaction_request(reactant_a, reactant_b, gamma_list, hm_orientations,
                           qc_params, site, proximity_factor):
    """Build the scientific request persisted with reaction restart state.

    The request captures the scientifically relevant inputs that must remain
    fixed across restarts: the gamma schedule, the backend configuration, the
    selected site constraint, the proximity factor, and a signature of both
    reactants.
    """
    backend_parameters = dict(qc_params)
    backend_parameters.pop("gamma", None)
    if backend_parameters.get("geometry_optimizer") == "geometric":
        # Prevent resuming survivors ranked with the legacy biased objective.
        backend_parameters["reaction_energy_convention"] = "physical-v1"
    return {
        "gamma_schedule": [float(value) for value in gamma_list],
        "bias_mode": "adaptive" if qc_params.get("bias_controller") == "adaptive" else "scheduled",
        "orientations": int(hm_orientations),
        "backend_parameters": backend_parameters,
        "site": None if site is None else list(site),
        "proximity_factor": float(proximity_factor),
        "reactants": [
            _molecule_signature(reactant_a),
            _molecule_signature(reactant_b),
        ],
    }


def _restore_product_registry(run_state):
    """Restore saved product identities so resumed runs deduplicate correctly."""
    saved_product_identities.clear()
    for job_name, (inchi, smiles) in run_state.saved_product_identities().items():
        saved_product_identities[job_name] = {"inchi": inchi, "smiles": smiles}


def _ensure_reactant_identity(run_state, reactant_a, reactant_b):
    """Persist the original separated-reactant identity used for product gating."""
    if "reactant_identity" not in run_state.data:
        run_state.data["reactant_identity"] = separated_reactant_identity(
            reactant_a,
            reactant_b,
        )
        run_state.save()
    return run_state.data["reactant_identity"]


def _is_known_product(identity):
    """Return whether the product's canonical molecular identity is known.

    The registry stores canonical InChI/SMILES pairs for already accepted
    products, so this check prevents emitting duplicate products when a run
    discovers the same structure through a different orientation or gamma.
    """
    return any(
        same_molecular_identity(identity, recorded_identity)
        for recorded_identity in saved_product_identities.values()
    )


def _should_continue_after_product(optimized_molecules):
    """Return whether the reaction loop should keep advancing gamma values."""
    return len(optimized_molecules) == 0 and bool(saved_product_identities)


def relax_without_afir_bias(molecule, qc_params):
    """Relax a bonded AFIR candidate on the unbiased physical objective.

    The molecule is written to temporary XYZ snapshots before and after the
    relaxation so callers can inspect the pre- and post-relaxation geometries
    if the optimization succeeds.
    """
    software = qc_params.get("software")
    if software is not None and not backend_supports_native_optimization(software):
        raise ValueError(
            f"Backend '{software}' does not advertise native optimization; "
            "unbiased reaction relaxation cannot be performed safely."
        )
    original_name = molecule.name
    release_attempt = int(qc_params.get("release_retry_attempt", 0))
    molecule.mol_to_xyz('trial_relax.xyz')
    # Each retry starts from a newly biased geometry. Give it a distinct native
    # optimization job name so the optimizer cannot return a cached relaxation
    # from the preceding attempt.
    molecule.name = 'relax' if release_attempt == 0 else f'relax_attempt_{release_attempt}'
    try:
        status = optimise(molecule, without_afir_bias(qc_params))
    finally:
        molecule.name = original_name
    if is_success(status):
        molecule.mol_to_xyz('result_relax.xyz')
    return status


def _release_evidence(qc_params):
    """Read the latest accepted-step release evidence from the optimizer state."""
    if qc_params.get("bias_controller") != "adaptive":
        return None
    job_name = qc_params.get("release_job_name")
    state_path = os.path.join(f"job_{job_name}", "pyar_geometric_state.json") if job_name else "pyar_geometric_state.json"
    if not os.path.isfile(state_path):
        return None
    try:
        with open(state_path) as fp:
            return json.load(fp).get("release")
    except (OSError, ValueError):
        reactor_logger.warning("Could not read adaptive release evidence from %s", state_path)
        return None


def _write_release_probe(pre_release, molecule, qc_params, status, evidence):
    """Persist the pre-release geometry and native free-relaxation outcome."""
    attempt = int(qc_params.get("release_retry_attempt", 0))
    pre_release_path = f"pre_release_attempt_{attempt}.xyz"
    pre_release.mol_to_xyz(pre_release_path)
    optimized = is_success(status) and molecule.coordinates is not None
    coordinates = np.asarray(molecule.coordinates, dtype=float) if optimized else np.empty((0, 3))
    post_bonds = infer_bonds(molecule.atoms_list, coordinates) if optimized else set()
    forming_pairs = {
        tuple(pair) for pair in (evidence or {}).get("forming_pairs", [])
    }
    survived = bool(forming_pairs) and forming_pairs.issubset(post_bonds) and all(
        _inside_covalent_radius_sum(pair, molecule.atoms_list, coordinates)
        for pair in forming_pairs
    )
    post_distances = [
        float(np.linalg.norm(coordinates[left] - coordinates[right]))
        for left, right in sorted(forming_pairs)
    ] if optimized else []
    probe = {
        "pre_release_geometry": pre_release_path,
        "native_optimizer_used": True,
        "release_attempt": attempt,
        "bias_alpha_margin": qc_params.get("bias_alpha_margin"),
        "free_relax_status": bool(is_success(status)),
        "post_relax_connectivity": [list(pair) for pair in sorted(post_bonds)],
        "post_relax_distance": post_distances,
        "release_survived": bool(survived),
        "release_state": "RELEASE_SURVIVED" if optimized and survived else "RELEASE_FAILED",
    }
    with open("release_probe.json", "w") as fp:
        json.dump(probe, fp, indent=2, sort_keys=True)
    with open("release_attempts.jsonl", "a") as fp:
        json.dump(probe, fp, sort_keys=True)
        fp.write("\n")
    state_path = os.path.join(f"job_{pre_release.name}", "pyar_geometric_state.json")
    if os.path.isfile(state_path):
        try:
            with open(state_path) as fp:
                state = json.load(fp)
            release = dict(state.get("release") or {})
            release["state"] = probe["release_state"]
            release["reason"] = (
                "native_free_relaxation_survived" if optimized and survived
                else "native_free_relaxation_failed_or_dissociated"
            )
            release["release_attempt"] = probe["release_attempt"]
            release["bias_alpha_margin"] = probe["bias_alpha_margin"]
            state["release"] = release
            temporary_path = state_path + ".tmp"
            with open(temporary_path, "w") as fp:
                json.dump(state, fp, indent=2, sort_keys=True)
            os.replace(temporary_path, state_path)
        except (OSError, ValueError):
            reactor_logger.warning("Could not update adaptive release state after free relaxation")
    return survived


def _set_release_outcome(job_name, state_name, reason):
    """Persist the chemical identity decision after a successful free probe."""
    state_path = os.path.join(f"job_{job_name}", "pyar_geometric_state.json")
    try:
        with open(state_path) as fp:
            state = json.load(fp)
        release = dict(state.get("release") or {})
        release.update(state=state_name, reason=reason)
        state["release"] = release
        temporary_path = state_path + ".tmp"
        with open(temporary_path, "w") as fp:
            json.dump(state, fp, indent=2, sort_keys=True)
        os.replace(temporary_path, state_path)
    except (OSError, ValueError):
        reactor_logger.warning("Could not persist chemical release outcome for %s", job_name)


def _write_release_retry_outcome(molecule, attempt, status, evidence, reason):
    """Add the final biased-retry result to the probe summary, preserving probe history."""
    geometry_path = f"release_retry_geometry_attempt_{int(attempt)}.xyz"
    molecule.mol_to_xyz(geometry_path)
    try:
        with open("release_probe.json") as fp:
            summary = json.load(fp)
    except (OSError, ValueError):
        summary = {}
    summary.update(
        retry_outcome="retry_not_candidate",
        retry_attempt=int(attempt),
        retry_optimization_status=str(status),
        retry_geometry=geometry_path,
        retry_reason=str(reason),
        retry_evidence=evidence,
    )
    with open("release_probe.json", "w") as fp:
        json.dump(summary, fp, indent=2, sort_keys=True)


def _prepare_release_retry(pre_release, qc_params, margin, attempt):
    """Restore the accepted checkpoint with a larger adaptive margin."""
    checkpoint_path = os.path.join(f"job_{pre_release.name}", "pyar_bias_controller_state.json")
    if not os.path.isfile(checkpoint_path):
        reactor_logger.warning("Cannot escalate release: controller checkpoint is missing")
        return None
    try:
        with open(checkpoint_path) as fp:
            checkpoint = json.load(fp)
    except (OSError, ValueError) as exc:
        reactor_logger.warning("Cannot escalate release: invalid controller checkpoint: %s", exc)
        return None
    configuration = checkpoint.get("configuration", {})
    checkpoint_qc = dict(configuration.get("qc_params", {}))
    checkpoint_qc["bias_alpha_margin"] = float(margin)
    configuration["qc_params"] = checkpoint_qc
    controller = checkpoint.get("controller", {})
    controller_configuration = dict(controller.get("configuration", {}))
    controller_configuration["safety_margin"] = float(margin)
    controller["configuration"] = controller_configuration
    checkpoint["configuration"] = configuration
    checkpoint["controller"] = controller
    # A failed unbiased probe is evidence that this contact is not yet a
    # releasable product. Keep the accumulated bias load, but restart contact
    # persistence at the saved geometry so the retry can continue driving
    # rather than immediately re-raising the same candidate.
    checkpoint.pop("release_tracker", None)
    checkpoint["positions_angstrom"] = np.asarray(pre_release.coordinates, dtype=float).tolist()
    try:
        temporary_path = checkpoint_path + ".tmp"
        with open(temporary_path, "w") as fp:
            json.dump(checkpoint, fp, indent=2, sort_keys=True)
        os.replace(temporary_path, checkpoint_path)
    except OSError as exc:
        reactor_logger.warning("Cannot write release retry checkpoint: %s", exc)
        return None
    retry_params = dict(qc_params)
    retry_params.update(
        bias_alpha_margin=float(margin),
        bias_controller_restart=True,
        trace_mode="append",
        release_retry_attempt=int(attempt),
        adaptive_max_segments=8,
        adaptive_suppress_release_candidate=True,
    )
    return retry_params


def _adaptive_release_probe(pre_release, molecule, qc_params, evidence):
    """Probe release, escalating only when the saved evidence supports retry."""
    current_pre_release = pre_release
    current_molecule = molecule
    current_params = dict(qc_params)
    current_evidence = evidence
    attempts = int(current_params.get("release_retry_attempt", 0))
    initial_attempts = attempts
    retry_limit = int(current_params.get("release_retry_limit", 2))
    margin_factor = float(current_params.get("release_margin_factor", 2.0))
    if retry_limit < 0 or not math.isfinite(margin_factor) or margin_factor <= 1.0:
        raise ValueError("release retry limit must be non-negative and margin factor must exceed 1")

    while True:
        _set_release_outcome(
            current_pre_release.name, "FREE_RELAX_PROBE", "native_free_relaxation_in_progress"
        )
        status = relax_without_afir_bias(current_molecule, current_params)
        survived = _write_release_probe(
            current_pre_release, current_molecule, current_params, status, current_evidence
        )
        if is_success(status) and survived:
            return status, current_molecule, current_pre_release, current_evidence, survived, attempts - initial_attempts

        persistence = int((current_evidence or {}).get("persistence_counter", 0))
        pair_keys = {
            str(tuple(sorted(pair))) for pair in (current_evidence or {}).get("forming_pairs", [])
        }
        retry_supported = bool(
            current_evidence
            and current_evidence.get("state") == "CANDIDATE"
            and pair_keys
            and persistence >= int(current_params.get("release_persistence", 3))
        )
        if attempts >= retry_limit or not retry_supported:
            return status, current_molecule, current_pre_release, current_evidence, survived, attempts - initial_attempts

        attempts += 1
        old_margin = float(current_params.get("bias_alpha_margin") or 0.001)
        new_margin = old_margin * margin_factor
        retry_params = _prepare_release_retry(
            current_pre_release, current_params, new_margin, attempts
        )
        if retry_params is None:
            return status, current_molecule, current_pre_release, current_evidence, survived, attempts - initial_attempts
        reactor_logger.info(
            "Release attempt %d failed; retrying from pre-release geometry with bias-alpha-margin=%s",
            attempts,
            new_margin,
        )
        retry_molecule = current_pre_release.copy()
        retry_status = optimise(retry_molecule, retry_params)
        if not is_usable(retry_status) or not retry_molecule.is_bonded():
            retry_evidence = _release_evidence(retry_params)
            _write_release_retry_outcome(
                retry_molecule, attempts, retry_status, retry_evidence,
                "retry_optimization_failed_or_geometry_unbonded",
            )
            return retry_status, retry_molecule, current_pre_release, retry_evidence, False, attempts - initial_attempts
        retry_evidence = _release_evidence(retry_params)
        if not retry_evidence or retry_evidence.get("state") != "CANDIDATE":
            _write_release_retry_outcome(
                retry_molecule, attempts, retry_status, retry_evidence,
                "release_evidence_not_candidate",
            )
            return retry_status, retry_molecule, current_pre_release, retry_evidence, False, attempts - initial_attempts
        current_pre_release = retry_molecule.copy()
        current_molecule = retry_molecule
        current_params = retry_params
        current_evidence = retry_evidence


def _orientation_final_coordinate_path(gamma_directory, orientation_directory, job_name):
    """Return the existing coordinate artifact that represents this orientation."""
    candidate_paths = (
        os.path.join(orientation_directory, "result_relax.xyz"),
        os.path.join(orientation_directory, "job_relax", "result_relax.xyz"),
        os.path.join(
            orientation_directory,
            f"job_{job_name}",
            f"result_{job_name}.xyz",
        ),
    )
    for candidate_path in candidate_paths:
        if os.path.exists(os.path.join(gamma_directory, candidate_path)):
            return candidate_path
    return "unavailable"


def initialize_reaction_run(reactant_a, reactant_b, gamma_min, gamma_max, hm_orientations,
                            qc_params, site, proximity_factor):
    """Prepare a reaction run and return the mutable workflow state.

    This function owns the restart contract for the reaction workflow. It
    either resumes an existing reaction state, migrates a legacy checkpoint,
    or creates a new ``reaction/`` tree with all trial geometries staged for
    the first gamma cycle.
    """
    current_workdir = os.getcwd()
    adaptive_mode = qc_params.get("bias_controller") == "adaptive"
    if gamma_max is None:
        raise ValueError("Adaptive reaction bias requires a finite --bias-max value")
    if adaptive_mode:
        if gamma_min is not None:
            reactor_logger.warning(
                "Ignoring --bias-min in adaptive mode; --bias-max is the single bias ceiling."
            )
        if not np.isfinite(gamma_max) or gamma_max < 0.0:
            raise ValueError("AFIR gamma maximum must be a finite non-negative number")
        requested_gamma_list = np.asarray([float(gamma_max)])
    else:
        if gamma_min is None:
            raise ValueError("Fixed or scheduled reaction bias requires --bias-min")
        requested_gamma_list = build_gamma_schedule(gamma_min, gamma_max)
    request = build_reaction_request(
        reactant_a,
        reactant_b,
        requested_gamma_list,
        hm_orientations,
        qc_params,
        site,
        proximity_factor,
    )
    sampling = sampling_configuration(
        number_of_orientations=int(hm_orientations),
        use_angles=(site is None and getattr(reactant_b, "number_of_atoms", 0) > 1),
    )
    run_state = ReactionRunState.load(current_workdir, request)
    if run_state is None:
        legacy_checkpoint = read_legacy_checkpoint(current_workdir)
        if legacy_checkpoint is not None:
            run_state = ReactionRunState.migrate_legacy(
                current_workdir, legacy_checkpoint, request, sampling=sampling
            )
            reactor_logger.warning(
                "Imported legacy jobs.pkl into reaction/state.json; "
                "legacy product deduplication history is unavailable."
            )

    if run_state is not None:
        reactor_logger.info('Reaction state detected: resuming reaction workflow')
        _ensure_reactant_identity(run_state, reactant_a, reactant_b)
        gamma_list = run_state.remaining_gamma_schedule()
        orientations_to_optimize = run_state.pending_molecules()
        os.chdir('reaction')
        cwd = os.getcwd()
        product_dir = f'{cwd}/products'
        _restore_product_registry(run_state)
        return current_workdir, cwd, run_state, gamma_list, orientations_to_optimize, product_dir

    if os.path.isdir('reaction'):
        raise ReactionStateError(
            "Existing reaction directory has no resumable state; "
            "preserve it and start in a new directory, or remove it explicitly."
        )
    os.makedirs('reaction')
    os.chdir('reaction')
    cwd = os.getcwd()

    reactor_logger.info('Starting reaction workflow')
    reactor_logger.info(
        "Reaction config: orientations=%s gamma_min=%s gamma_max=%s site=%s proximity_factor=%s",
        hm_orientations, gamma_min, gamma_max, site, proximity_factor,
    )
    reactor_logger.debug("Reaction qc_params=%s", qc_params)
    reactor_logger.debug(f'Current working directory: {cwd}')

    software = qc_params['software']
    print_header(gamma_max, gamma_min, hm_orientations, software)
    product_dir = f'{cwd}/products'
    reactor_logger.debug(f'Product directory: {product_dir}')
    file_manager.make_directories(product_dir)
    file_manager.make_directories('trial_geometries')
    os.chdir('trial_geometries')
    if site is None:
        all_orientations = trial_generation.create_trial_geometries(
            'geom',
            reactant_a,
            reactant_b,
            hm_orientations,
            site,
        )
    else:
        all_orientations = pyar.scan.generate_guess_for_bonding(
            'geom',
            reactant_a,
            reactant_b,
            site[0],
            site[1],
            hm_orientations,
            d_scale=proximity_factor,
        )

    os.chdir(cwd)

    gamma_list = requested_gamma_list
    orientations_to_optimize = all_orientations[:]
    run_state = ReactionRunState.create(
        current_workdir,
        request,
        orientations_to_optimize,
        (reactant_a, reactant_b),
        sampling=sampling,
    )
    _ensure_reactant_identity(run_state, reactant_a, reactant_b)
    _restore_product_registry(run_state)
    return current_workdir, cwd, run_state, gamma_list, orientations_to_optimize, product_dir


def react(reactant_a, reactant_b, gamma_min, gamma_max, hm_orientations, qc_params,
          site, proximity_factor):
    """Run the reaction-search workflow for two reactants.

    The workflow iterates over the gamma schedule, optimizes each orientation,
    records products and trace summaries, and returns a structured result that
    summarizes the final reaction state.
    """
    sampling = sampling_configuration(
        number_of_orientations=int(hm_orientations),
        use_angles=(site is None and getattr(reactant_b, "number_of_atoms", 0) > 1),
    )
    workdir, cwd, run_state, gamma_list, orientations_to_optimize, product_dir = initialize_reaction_run(
        reactant_a,
        reactant_b,
        gamma_min,
        gamma_max,
        hm_orientations,
        qc_params,
        site,
        proximity_factor,
    )

    for gamma in gamma_list:
        adaptive_mode = qc_params.get("bias_controller") == "adaptive"
        gamma_id = "adaptive" if adaptive_mode else format_gamma_id(gamma)
        reactor_logger.info(f"Gamma cycle path: reaction/gamma_{gamma_id}")
        reactor_logger.info(f'Gamma cycle start: {gamma_id}')
        gamma_home = f'{cwd}/gamma_{gamma_id}'
        if not os.path.exists(gamma_home):
            file_manager.make_directories(gamma_home)
        os.chdir(gamma_home)

        gamma_qc_params = with_gamma(qc_params, gamma)
        optimized_molecules = optimize_all(
            gamma_id,
            orientations_to_optimize,
            run_state,
            product_dir,
            gamma_qc_params,
        )

        reactor_logger.info(
            f"Gamma cycle optimized geometries: {len(optimized_molecules)}")
        if len(optimized_molecules) == 0:
            run_state.complete_cycle(gamma, [])
            if _should_continue_after_product(optimized_molecules):
                reactor_logger.info(
                    "Reaction search found %d unique product(s); continuing "
                    "through the remaining gamma schedule.",
                    len(saved_product_identities),
                )
                orientations_to_optimize = []
                continue
            else:
                reactor_logger.info("No orientations left for next gamma cycle.")
                run_state.finish("completed_no_candidates")
                result_status = "completed_no_candidates"
            os.chdir(workdir)
            return _build_reaction_result(
                workdir,
                result_status,
                product_dir,
                run_state,
                gamma_list,
                orientations_to_optimize,
                sampling,
            )
        if len(optimized_molecules) == 1:
            orientations_to_optimize = optimized_molecules[:]
        else:
            orientations_to_optimize = clustering.remove_similar(
                optimized_molecules)
        run_state.complete_cycle(gamma, orientations_to_optimize)
        reactor_logger.info(f"Products found so far: {len(saved_product_identities)}")
        reactor_logger.info(f"Next cycle candidate geometries: {len(orientations_to_optimize)}")

        reactor_logger.debug("the keys of the molecules for next gamma cycle")
        for this_orientation in orientations_to_optimize:
            reactor_logger.debug(f"{this_orientation.name}")

    os.chdir(workdir)
    terminal_status = "completed_products_found" if saved_product_identities else "completed"
    run_state.finish(terminal_status)
    reactor_logger.info("Reaction workflow completed. State retained in reaction/state.json.")
    return _build_reaction_result(
        workdir,
        terminal_status,
        product_dir,
        run_state,
        gamma_list,
        orientations_to_optimize,
        sampling,
    )


def optimize_all(gamma_id, orientations, run_state, product_dir, qc_param):
    """Optimize all trial geometries for one gamma cycle.

    Each orientation is written to its own job directory, optimized with the
    current gamma, and then either retained for the next gamma value or
    promoted to a unique product if it survives the unbiased relaxation step.
    """
    gamma = qc_param['gamma']
    cwd = os.getcwd()
    table_of_optimized_molecules = (
        run_state.current_survivor_molecules()
        if run_state is not None
        else []
    )
    reactant_identity = (
        getattr(run_state, "data", {}).get("reactant_identity")
        if run_state is not None
        else None
    )
    pending_orientations = list(orientations)

    def record_orientation_completion(job_name, status):
        """Persist the processed job and any orientations still pending."""
        pending_orientations.pop(0)
        if run_state is not None:
            run_state.record_job(
                job_name,
                gamma,
                status,
                pending_orientations,
                table_of_optimized_molecules,
            )
        reactor_logger.info(
            "Orientation completed! status=%s final_coordinate=%s",
            status,
            _orientation_final_coordinate_path(
                cwd,
                orientations_home,
                job_name,
            ),
        )

    for this_molecule in orientations:
        job_key = this_molecule.name
        reactor_logger.info(f'   Orientation: {job_key}')
        o_key = f"_{job_key[-8:]}"
        orientations_home = f'orientation{o_key}'
        file_manager.make_directories(orientations_home)
        os.chdir(orientations_home)
        job_name = gamma_id + o_key
        this_molecule.name = job_name
        reactor_logger.info(f'Optimizing {this_molecule.name}')
        start_xyz_file_name = f'trial_{this_molecule.name}.xyz'
        this_molecule.mol_to_xyz(start_xyz_file_name)
        reference_xyz_file_name = f'reactants_{this_molecule.name}.xyz'
        write_disconnected_reference(this_molecule, reference_xyz_file_name)
        start_identity = molecule_identity_from_xyz(reference_xyz_file_name)
        orientation_qc_param = dict(qc_param)
        orientation_qc_param["release_job_name"] = job_name
        status = optimise(this_molecule, orientation_qc_param)
        this_molecule.name = job_name
        reactor_logger.info('Optimization step completed')
        if is_usable(status):
            before_relax = this_molecule.copy()
            reactor_logger.info("Energy E({}): {:12.7f}".format(job_name, this_molecule.energy))

            if this_molecule.is_bonded():
                evidence = _release_evidence(orientation_qc_param)
                release_candidate = (
                    evidence is not None and evidence.get("state") == "CANDIDATE"
                )
                if qc_param.get("bias_controller") == "adaptive" and not release_candidate:
                    reactor_logger.info(
                        "%s has contact but no conservative release candidate; retaining accepted geometry",
                        job_name,
                    )
                    table_of_optimized_molecules.append(before_relax)
                    record_orientation_completion(job_name, "contact_not_release_ready")
                    os.chdir(cwd)
                    continue
                reactor_logger.info("Close contacts detected; running unbiased relaxation (gamma=0.0)")
                release_survived = True
                release_attempts = 0
                if qc_param.get("bias_controller") == "adaptive":
                    (
                        status,
                        this_molecule,
                        before_relax,
                        evidence,
                        release_survived,
                        release_attempts,
                    ) = _adaptive_release_probe(
                        before_relax, this_molecule, orientation_qc_param, evidence
                    )
                else:
                    status = relax_without_afir_bias(this_molecule, qc_param)
                if is_success(status):
                    if qc_param.get("bias_controller") == "adaptive" and not release_survived:
                        reactor_logger.info("Native free relaxation lost the forming connectivity")
                        table_of_optimized_molecules.append(before_relax)
                        record_orientation_completion(job_name, "release_failed")
                        os.chdir(cwd)
                        continue
                    current_identity = molecule_identity_from_xyz('result_relax.xyz')
                    current_inchi = current_identity["inchi"]
                    current_smile = current_identity["smiles"]

                    reactor_logger.info('Relaxation completed')
                    reactor_logger.info("Checking product formation using SMILES and InChI")

                    reactor_logger.info(f"Start SMILE: {start_identity['smiles']} Current SMILE: {current_smile}")
                    reactor_logger.info(f"Start InChi: {start_identity['inchi']} Current InChi: {current_inchi}")

                    reference_identity = reactant_identity or start_identity
                    if reactant_identity is not None:
                        reactor_logger.info(
                            f"Reactant SMILE: {reactant_identity['smiles']} Current SMILE: {current_smile}"
                        )
                        reactor_logger.info(
                            f"Reactant InChi: {reactant_identity['inchi']} Current InChi: {current_inchi}"
                        )

                    # Product validity is determined against the original
                    # separated reactants, not a distorted higher-gamma
                    # survivor that may serialize differently.
                    if not reaction_product_changed(reference_identity, current_identity):
                        if qc_param.get("bias_controller") == "adaptive":
                            _set_release_outcome(
                                job_name, "RELEASE_FAILED", "relaxed_identity_matches_reactants"
                            )
                        table_of_optimized_molecules.append(before_relax)
                        reactor_logger.info(f'{job_name} kept for higher-gamma optimization')
                    else:
                        if qc_param.get("bias_controller") == "adaptive":
                            _set_release_outcome(
                                job_name, "PRODUCT_CONFIRMED", "free_relaxation_and_identity_changed"
                            )
                        reactor_logger.info(
                            "Relaxed identity changed from the starting structure."
                        )

                        reactor_logger.info("Checking whether product is new")
                        if _is_known_product(current_identity):
                            reactor_logger.info("Product matches an existing product; discarded")
                            product_status = "duplicate_product"
                        else:
                            reactor_logger.info("New product detected; saving")
                            saved_product_identities[job_name] = current_identity
                            shutil.copy('result_relax.xyz', f'{product_dir}/{job_name}.xyz')
                            trace_directory = os.path.join(os.getcwd(), f'job_{job_name}')
                            trace_summary = None
                            try:
                                trace_summary = reaction_analysis.analyse_reaction_trace(trace_directory)
                                reactor_logger.info(
                                    "Reaction trace analyzed for %s: %s",
                                    job_name,
                                    trace_summary["candidate_ts_directory"] if trace_summary else "no trace records",
                                )
                                if trace_summary is not None:
                                    reactor_logger.info(
                                        "Trace candidates for %s: highest_backend=%s pre_product=%s first_topology_change=%s highest_total=%s",
                                        job_name,
                                        trace_summary.get("highest_backend_energy_index"),
                                        trace_summary.get("pre_product_index"),
                                        trace_summary.get("first_topology_change_index"),
                                        trace_summary.get("highest_total_energy_index"),
                                    )
                            except Exception:
                                reactor_logger.exception(
                                    "Reaction trace analysis failed for %s", job_name
                                )
                            if run_state is not None:
                                run_state.record_product(
                                    job_name,
                                    gamma,
                                    current_inchi,
                                    current_smile,
                                    f'{product_dir}/{job_name}.xyz',
                                    trace_summary=trace_summary,
                                )
                            product_status = "new_product"
                        os.chdir(cwd)
                        record_orientation_completion(job_name, product_status)
                        continue
                elif is_cycle_exceeded(status):
                    table_of_optimized_molecules.append(before_relax)
                    reactor_logger.info(f'{job_name} kept for higher-gamma optimization')
                elif qc_param.get("bias_controller") == "adaptive":
                    table_of_optimized_molecules.append(before_relax)
                    reactor_logger.info(
                        "%s native free relaxation failed; preserving pre-release geometry",
                        job_name,
                    )

            else:
                table_of_optimized_molecules.append(this_molecule)
                reactor_logger.info('No close contacts found')
                reactor_logger.info(f'{job_name} kept for higher-gamma optimization')

        record_orientation_completion(job_name, status)
        os.chdir(cwd)
        sys.stdout.flush()
    return table_of_optimized_molecules


def main():
    """Compatibility entry point retained for older imports and scripts."""
    return None


if __name__ == "__main__":
    main()
