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
import os
import shutil
import sys

import numpy as np
import pyar.scan
from pyar import file_manager
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
    forcing ``gamma=0.0`` so the candidate can be re-optimized without the
    AFIR bias.
    """
    return with_gamma(qc_params, 0.0)


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
    return {
        "gamma_schedule": [float(value) for value in gamma_list],
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
    original_name = molecule.name
    molecule.mol_to_xyz('trial_relax.xyz')
    molecule.name = 'relax'
    try:
        status = optimise(molecule, without_afir_bias(qc_params))
    finally:
        molecule.name = original_name
    if is_success(status):
        molecule.mol_to_xyz('result_relax.xyz')
    return status


def initialize_reaction_run(reactant_a, reactant_b, gamma_min, gamma_max, hm_orientations,
                            qc_params, site, proximity_factor):
    """Prepare a reaction run and return the mutable workflow state.

    This function owns the restart contract for the reaction workflow. It
    either resumes an existing reaction state, migrates a legacy checkpoint,
    or creates a new ``reaction/`` tree with all trial geometries staged for
    the first gamma cycle.
    """
    current_workdir = os.getcwd()
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
        gamma_id = format_gamma_id(gamma)
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
        status = optimise(this_molecule, qc_param)
        this_molecule.name = job_name
        reactor_logger.info('Optimization step completed')
        if is_usable(status):
            before_relax = this_molecule.copy()
            reactor_logger.info("Energy E({}): {:12.7f}".format(job_name, this_molecule.energy))

            if this_molecule.is_bonded():
                reactor_logger.info("Close contacts detected; running unbiased relaxation (gamma=0.0)")
                status = relax_without_afir_bias(this_molecule, qc_param)
                if is_success(status):
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
                        table_of_optimized_molecules.append(before_relax)
                        reactor_logger.info(f'{job_name} kept for higher-gamma optimization')
                    else:
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
                                        "Trace candidates for %s: highest_backend=%s pre_product=%s max_bond_change=%s highest_total=%s",
                                        job_name,
                                        trace_summary.get("highest_backend_energy_index"),
                                        trace_summary.get("pre_product_index"),
                                        trace_summary.get("max_bond_change_index"),
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
