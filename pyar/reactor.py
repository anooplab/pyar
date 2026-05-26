"""Reaction workflow orchestration for PyAR."""

import logging
import os
import shutil
import sys

import numpy as np
import pyar.scan
from pyar import trial_generation, file_manager
from pyar.data_analysis import clustering
from pyar.optimiser import is_cycle_exceeded, is_success, is_usable, optimise
from pyar.workflow_results import ReactionResult
from pyar.reaction_identity import (
    molecule_identity_from_xyz,
    same_molecular_identity,
    write_disconnected_reference,
)
from pyar.state.reaction import ReactionRunState, ReactionStateError, read_legacy_checkpoint

reactor_logger = logging.getLogger('pyar.reactor')

saved_product_identities = {}


def _build_reaction_result(workdir, status, product_dir, run_state, gamma_list, orientations):
    """Package the current reaction outcome as a structured result."""
    return ReactionResult(
        workflow="reaction",
        status=status,
        run_directory=str(os.path.join(workdir, "reaction")),
        state_path=str(os.path.join(workdir, "reaction", "state.json")),
        selected_paths=tuple(product["path"] for product in run_state.data.get("products", [])),
        metadata={
            "gamma_schedule": tuple(float(value) for value in gamma_list),
            "products": tuple(run_state.data.get("products", [])),
            "product_directory": product_dir,
            "remaining_orientations": len(orientations),
        },
    )


def print_header(gamma_max, gamma_min, hm_orientations, software):
    """Log a concise header for a reactor run."""
    reactor_logger.info("==================== PyAR Reaction Workflow ====================")
    reactor_logger.info(f"Gamma range: {gamma_min} to {gamma_max}")
    reactor_logger.info(f"Orientations: {hm_orientations}")
    reactor_logger.info(f"Software: {software}")
    reactor_logger.info("===============================================================")


def with_gamma(qc_params, gamma):
    """Return a copy of ``qc_params`` with the current gamma value applied."""
    updated_qc_params = dict(qc_params)
    updated_qc_params['gamma'] = gamma
    return updated_qc_params


def without_afir_bias(qc_params):
    """Return parameters for unbiased physical relaxation of a reaction candidate."""
    return with_gamma(qc_params, 0.0)


def build_gamma_schedule(gamma_min, gamma_max, steps=10):
    """Build the numeric AFIR gamma schedule used by the reaction workflow."""
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
    """Format a gamma value for directory names and readable job labels."""
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
    """Return stable input geometry metadata used to validate restarts."""
    return {
        "atoms": list(molecule.atoms_list),
        "coordinates": np.asarray(molecule.coordinates, dtype=float).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
    }


def build_reaction_request(reactant_a, reactant_b, gamma_list, hm_orientations,
                           qc_params, site, proximity_factor):
    """Build the scientific request persisted with reaction restart state."""
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


def _is_known_product(identity):
    """Return whether the product's canonical molecular identity is known."""
    return any(
        same_molecular_identity(identity, recorded_identity)
        for recorded_identity in saved_product_identities.values()
    )


def relax_without_afir_bias(molecule, qc_params):
    """Relax a bonded AFIR candidate on the unbiased physical objective."""
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
    """Prepare a reaction run and return the mutable workflow state."""
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
    run_state = ReactionRunState.load(current_workdir, request)
    if run_state is None:
        legacy_checkpoint = read_legacy_checkpoint(current_workdir)
        if legacy_checkpoint is not None:
            run_state = ReactionRunState.migrate_legacy(
                current_workdir, legacy_checkpoint, request
            )
            reactor_logger.warning(
                "Imported legacy jobs.pkl into reaction/state.json; "
                "legacy product deduplication history is unavailable."
            )

    if run_state is not None:
        reactor_logger.info('Reaction state detected: resuming reaction workflow')
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
    )
    _restore_product_registry(run_state)
    return current_workdir, cwd, run_state, gamma_list, orientations_to_optimize, product_dir


def react(reactant_a, reactant_b, gamma_min, gamma_max, hm_orientations, qc_params,
          site, proximity_factor):
    """Run the reaction-search workflow for two reactants."""
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
            if saved_product_identities:
                reactor_logger.info(
                    "Reaction search completed with %d unique product(s); "
                    "no candidates require higher-gamma optimization.",
                    len(saved_product_identities),
                )
                run_state.finish("completed_products_found")
                result_status = "completed_products_found"
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
    )


def optimize_all(gamma_id, orientations, run_state, product_dir, qc_param):
    """Optimize all trial geometries for one gamma cycle."""
    gamma = qc_param['gamma']
    cwd = os.getcwd()
    table_of_optimized_molecules = (
        run_state.current_survivor_molecules()
        if run_state is not None
        else []
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

                    if same_molecular_identity(start_identity, current_identity):
                        table_of_optimized_molecules.append(before_relax)
                        reactor_logger.info(f'{job_name} kept for higher-gamma optimization')

                    else:
                        reactor_logger.info("Geometry differs from starting structure.")

                        reactor_logger.info("Checking whether product is new")
                        if _is_known_product(current_identity):
                            reactor_logger.info("Product matches an existing product; discarded")
                            product_status = "duplicate_product"

                        else:
                            reactor_logger.info("New product detected; saving")
                            saved_product_identities[job_name] = current_identity
                            shutil.copy('result_relax.xyz', f'{product_dir}/{job_name}.xyz')
                            if run_state is not None:
                                run_state.record_product(
                                    job_name,
                                    gamma,
                                    current_inchi,
                                    current_smile,
                                    f'{product_dir}/{job_name}.xyz',
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
    """Entry point placeholder retained for compatibility."""
    return None


if __name__ == "__main__":
    main()
