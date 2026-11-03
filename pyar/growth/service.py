"""Shared growth-workflow helpers for aggregate and solvation paths."""

from __future__ import annotations

import copy
import json
import logging
import os
import random
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from pyar.growth.errors import GrowthStopped
from pyar.growth.policy import (GrowthPolicy, GrowthRun, SeedRun, SeedOptimization,
    OptimizationStages, FinalSelectionRun, FinalGroupSelection)

import numpy as np

from pyar import file_manager
from pyar.core.molecule import Molecule, atomic_data
from pyar.io import xyz as molecule_io
from pyar.selection import clustering
from pyar.selection import reports as selection_reports
from pyar.sampling import trial_generator as trial_generation
from pyar.optimiser import is_cycle_exceeded, is_success, optimise
from pyar.backend_errors import BackendExecutionError
from pyar.optimization_service import optimize_workflow_once
from pyar.structure_comparison.coordinate_graph import analyze_growth_transition
from pyar.workflow_results import AggregateResult, ReactionResult

aggregator_logger = logging.getLogger("pyar.workflows.aggregate")

_WORKFLOW_RESULTS = {
    "aggregate": AggregateResult,
    "reaction": ReactionResult,
}

_CONNECTIVITY_POLICIES = {"off", "prefer", "strict"}



def _write_growth_structure_analysis(input_parts, outputs, filename):
    """Write best-effort coordinate-only input/product diagnostics."""
    if not outputs:
        return
    try:
        report = {
            "stage": "pre-similarity-selection",
            "method": "coordinate-only-adjacency",
            "structures": [
                {
                    "name": str(getattr(output, "name", "<unnamed>")),
                    **analyze_growth_transition(input_parts, output),
                }
                for output in outputs
            ],
        }
        with open(filename, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
    except (KeyError, ValueError, TypeError) as exc:
        # This diagnostic must never block structure generation or selection.
        aggregator_logger.warning("Coordinate-only structural analysis skipped: %s", exc)


def _optimize_growth_batch(
    molecules, qc_params, workflow_role, failure_filename, *, directory=None
):
    """Attempt each independent candidate and persist backend failures."""
    directory = Path(directory or Path.cwd()).resolve()
    statuses = []
    failures = []
    failure_path = Path(str(directory / failure_filename))
    if failure_path.exists():
        try:
            failures = json.loads(failure_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            failures = []
    for molecule in molecules:
        try:
            status, _diagnostic, _diagnostic_path = optimize_workflow_once(
                molecule,
                qc_params,
                workflow_role=workflow_role,
                optimizer=optimise,
                directory=directory,
            )
        except BackendExecutionError as exc:
            statuses.append(False)
            failures.append({"molecule": molecule.name, "workflow_role": workflow_role, "error": str(exc)})
            failure_path.write_text(json.dumps(failures, indent=2) + "\n", encoding="utf-8")
            aggregator_logger.error("Backend failed for %s (%s): %s", molecule.name, workflow_role, exc)
        else:
            statuses.append(status)
    return statuses


def _select_growth_seed_pool(optimized_molecules, selected_from_restart, run):
    """Select fresh candidates and combine them with the resumed selected pool."""
    selection_diagnostics = {}
    selected_seeds = clustering.choose_geometries(
        optimized_molecules,
        maximum_number_of_seeds=run.policy.maximum_number_of_seeds,
        persist_basin_memory=not run.stages.staged,
        group_basin_by_stoichiometry=False,
        connectivity_policy=run.connectivity_policy,
        feature=run.policy.selection.feature,
        algorithm=run.policy.selection.algorithm,
        distance_metric=run.policy.selection.distance,
        system_type=run.policy.selection.system_type,
        diagnostics=selection_diagnostics,
        directory=run.directory,
    )
    _write_growth_selection_diagnostics(
        selection_diagnostics,
        str(run.directory / os.path.join("selected", "selection_diagnostics.json")),
    )
    for molecule in selected_seeds:
        molecule.connectivity_policy_hint = run.connectivity_policy
    selected_seeds = selected_from_restart + selected_seeds
    if len(selected_seeds) > run.policy.maximum_number_of_seeds:
        aggregator_logger.warning(
            "Restart plus current selection produced %d seeds; keeping lowest-energy %d.",
            len(selected_seeds),
            run.policy.maximum_number_of_seeds,
        )
        selected_seeds = sorted(
            selected_seeds, key=lambda molecule: float(molecule.energy)
        )[: run.policy.maximum_number_of_seeds]
    return selected_seeds


def _write_growth_selection_diagnostics(diagnostics, filename):
    """Persist selection and basin-memory provenance without gating workflow progress."""
    if not diagnostics:
        return
    target = Path(filename)
    temporary_path = None
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=target.parent,
            prefix=f".{target.name}-", suffix=".tmp", delete=False,
        ) as stream:
            json.dump(diagnostics, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
            temporary_path = stream.name
        os.replace(temporary_path, target)
    except (OSError, TypeError, ValueError) as exc:
        aggregator_logger.warning("Selection diagnostics could not be written to %s: %s", target, exc)
    finally:
        if temporary_path and os.path.exists(temporary_path):
            os.unlink(temporary_path)


def workflow_run_directory(root_directory, workflow_name):
    """Return the absolute workflow run directory for the current invocation."""
    return str(Path(root_directory).resolve() / workflow_name)


def workflow_state_path(run_directory):
    """Return the canonical JSON state file path for a workflow run."""
    return str(Path(run_directory) / "state.json")


# Public metadata builder: retain keyword controls; there is no forwarding chain.
def sampling_configuration(
    *,
    direction_method="fibonacci",
    rotation_method="halton",
    number_of_orientations=None,
    sequence_offset=0,
    oversample_factor=8,
    seed=None,
    use_angles=None,
):
    """Return the sampling metadata recorded in workflow state and results."""
    configuration = {
        "direction_method": direction_method,
        "rotation_method": rotation_method,
        "sequence_offset": int(sequence_offset),
        "oversample_factor": int(oversample_factor),
    }
    if number_of_orientations is not None:
        configuration["number_of_orientations"] = int(number_of_orientations)
    if seed is not None:
        configuration["seed"] = int(seed)
    if use_angles is not None:
        configuration["use_angles"] = bool(use_angles)
    return configuration


def stopped_workflow_result(workflow_name, root_directory):
    """Return the standard structured result for a workflow stop request."""
    result_cls = _WORKFLOW_RESULTS[workflow_name]
    run_directory = workflow_run_directory(root_directory, workflow_name)
    return result_cls(
        workflow=workflow_name,
        status="stopped",
        run_directory=run_directory,
        state_path=workflow_state_path(run_directory),
    )


def _format_path(path):
    return "".join(m.name for m in path)


def _stoichiometry_label(molecule):
    """Return a compact chemical formula label for a molecule."""
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


def _atom_count(molecule):
    """Return the atom count for a molecule-like object."""
    number_of_atoms = getattr(molecule, "number_of_atoms", None)
    if number_of_atoms is not None:
        return int(number_of_atoms)
    return len(getattr(molecule, "atoms_list", ()))


def _connectivity_hint(molecule):
    """Return any explicit connectivity preference recorded on a molecule."""
    for attribute in (
        "connectivity_policy_hint",
        "connectivity_policy",
        "growth_mode",
        "workflow_kind",
        "origin",
        "source_kind",
    ):
        value = getattr(molecule, attribute, None)
        if value in _CONNECTIVITY_POLICIES:
            return value
    return None


def _is_formula_generated(molecule):
    return getattr(molecule, "connectivity_policy_hint", None) == "prefer" or getattr(molecule, "source_kind", None) == "formula"


def _is_radical_like(molecule):
    """Return True for simple radical-like metadata cues."""
    multiplicity = getattr(molecule, "multiplicity", None)
    if multiplicity is not None:
        try:
            if int(multiplicity) % 2 == 0:
                return True
        except (TypeError, ValueError):
            pass

    atomic_numbers = getattr(molecule, "atomic_number", None)
    if atomic_numbers and len(atomic_numbers) == 1:
        try:
            charge = int(getattr(molecule, "charge", 0) or 0)
            electron_count = sum(int(number) for number in atomic_numbers) - charge
        except (TypeError, ValueError):
            return False
        return electron_count % 2 == 1

    return False


def resolve_connectivity_policy_for_growth(seed, monomer, qc_params=None):
    """Resolve how growth workflows should treat covalent connectivity."""
    qc_params = qc_params or {}
    requested_policy = qc_params.get("connectivity_policy")
    if requested_policy is not None:
        requested_policy = str(requested_policy).lower()
        if requested_policy == "auto":
            requested_policy = None
        else:
            if requested_policy not in _CONNECTIVITY_POLICIES:
                raise ValueError(
                    f"Unknown connectivity policy: {requested_policy!r}. "
                    "Expected one of 'off', 'prefer', 'strict', or 'auto'."
                )
            return requested_policy

    for molecule in (seed, monomer):
        hint = _connectivity_hint(molecule)
        if hint is not None:
            return hint

    seed_atoms = _atom_count(seed)
    monomer_atoms = _atom_count(monomer)

    if _is_formula_generated(seed) or _is_formula_generated(monomer):
        return "prefer"

    if _is_radical_like(seed) or _is_radical_like(monomer):
        return "prefer"

    if seed_atoms == 1 and monomer_atoms == 1:
        return "prefer"

    if monomer_atoms == 1 and (
        getattr(seed, "growth_kind", None) == "atomic_cluster"
        or getattr(seed, "aggregate_kind", None) == "atomic_cluster"
        or getattr(seed, "workflow_kind", None) == "atomic_cluster"
    ):
        return "prefer"

    if seed_atoms == 1 and monomer_atoms > 1:
        return "off"

    if seed_atoms > 1 or monomer_atoms > 1:
        return "off"

    return "off"


def resolve_connectivity_policy_for_aggregate(molecules):
    """Resolve the aggregate-level connectivity policy for final clustering."""
    if any(_is_formula_generated(molecule) for molecule in molecules):
        return "prefer"
    if any(_is_radical_like(molecule) for molecule in molecules):
        return "prefer"
    if molecules and all(_atom_count(molecule) == 1 for molecule in molecules):
        return "prefer"
    return "off"


def _snapshot_selected_geometries(
    selected_seeds,
    output_root="selected",
    summary_lines=None,
    group_by_stoichiometry=False,
):
    """Persist the selected geometries to disk."""
    if not selected_seeds:
        return None

    snapshot_dir = output_root
    if group_by_stoichiometry:
        stoichiometry = _stoichiometry_label(selected_seeds[0])
        snapshot_dir = os.path.join(output_root, f"stoichiometry_{stoichiometry}")
    else:
        stoichiometry = _stoichiometry_label(selected_seeds[0])
    if not os.path.exists(snapshot_dir):
        file_manager.make_directories(snapshot_dir)
    else:
        for old_result in Path(snapshot_dir).glob("result_*.xyz"):
            old_result.unlink()

    for molecule in selected_seeds:
        target_file = os.path.join(snapshot_dir, f"result_{molecule.name}.xyz")
        molecule_io.write_xyz(molecule, target_file)
        molecule.relative_path = os.path.relpath(target_file, output_root)

    if all(m.energy is not None for m in selected_seeds):
        selection_reports.print_energy_table(selected_seeds, title="Selected energy table:")

    if summary_lines:
        summary_path = os.path.join(snapshot_dir, "README.txt")
        with open(summary_path, "w") as fp:
            fp.write(f"Selected geometries for stoichiometry {stoichiometry}\n")
            fp.write("=" * 60 + "\n")
            for line in summary_lines:
                fp.write(f"{line}\n")

    aggregator_logger.info(
        "Saved %d selected geometries to %s",
        len(selected_seeds),
        snapshot_dir,
    )
    return snapshot_dir


def _resolve_orientation_count(hm_orientations):
    """Return the integer orientation count requested by the caller."""
    if hm_orientations == "auto":
        return 8
    return int(hm_orientations)


def _seed_specific_orientation_count(seed, requested_count):
    """Return the trial count for a seed, honoring monoatomic special cases."""
    return 1 if len(seed) == 1 else requested_count


def _discover_selected_result_files(aggregate_root):
    """Return result files collected from all pathway-level selected folders."""
    root_path = Path(aggregate_root)
    result_files = set()
    for ag_dir in root_path.glob("ag_*"):
        if not ag_dir.is_dir():
            continue
        result_files.update(ag_dir.rglob("selected/result_*.xyz"))
    return sorted(result_files)


def _selection_name_from_path(result_file, aggregate_root):
    """Build a filesystem-safe name from a collected result path."""
    result_path = Path(result_file)
    try:
        relative_path = result_path.relative_to(Path(aggregate_root))
    except ValueError:
        relative_path = result_path
    return str(relative_path.with_suffix("")).replace(os.sep, "_")


def _finalize_selected_geometries(run: FinalSelectionRun, policy: GrowthPolicy):
    """Group original pathway evidence, select each composition and publish it."""
    result_files = (
        _discover_selected_result_files(run.directory)
        if run.source_records is None
        else sorted(run.source_records, key=lambda path: Path(run.source_records[path]))
    )
    if not result_files:
        aggregator_logger.info("No pathway-level selected results were found for final clustering.")
        return []
    groups = _group_final_candidates(run, result_files)
    selected_root = os.path.join(run.directory, "selected")
    if not os.path.exists(selected_root):
        file_manager.make_directories(selected_root)
    selected = []
    for stoichiometry, molecules in sorted(groups.items()):
        aggregator_logger.info(
            "Final clustering for stoichiometry %s from %d pathway results",
            stoichiometry,
            len(molecules),
        )
        group = _select_final_group(run, policy, stoichiometry, molecules)
        selected.extend(group.selected)
        _publish_final_group(run, policy, group)
    return selected


def _group_final_candidates(run, result_files):
    grouped = defaultdict(list)
    for each_file in result_files:
        molecule = Molecule.from_xyz(str(each_file))
        molecule.energy = (
            selection_reports.read_energy_from_xyz_file(str(each_file))
            if run.energy_ranked
            else None
        )
        molecule.name = (
            _selection_name_from_path(each_file, run.directory)
            if run.source_records is None
            else _selection_name_from_path(run.source_records[each_file], ".")
        )
        grouped[_stoichiometry_label(molecule)].append(molecule)
    return grouped

def _select_final_group(run, policy, stoichiometry, molecules):
    selection_diagnostics = {}
    if run.energy_ranked:
        selected = clustering.choose_geometries(
            molecules,
            maximum_number_of_seeds=policy.maximum_number_of_seeds,
            apply_basin_memory=False,
            algorithm=policy.selection.algorithm,
            feature=policy.selection.feature,
            distance_metric=policy.selection.distance,
            connectivity_policy=policy.connectivity_policy or "off",
            system_type=policy.selection.system_type,
            diagnostics=selection_diagnostics,
            directory=run.directory,
        )
    else:
        selected = select_geometry_only(molecules, policy, directory=run.directory)
        selection_diagnostics.update(selection_mode="geometry-diversity", energy_ranked=False,
                                     input_count=len(molecules), selected_count=len(selected))
    return FinalGroupSelection(
        stoichiometry, molecules, selected, selection_diagnostics
    )

def _publish_final_group(run, policy, group):
    selected_root = os.path.join(run.directory, "selected")
    snapshot_directory = _snapshot_selected_geometries(
        group.selected,
        output_root=selected_root,
        summary_lines=[
            f"Stoichiometry: {group.stoichiometry}",
            "Selection mode: final cross-path clustering",
            f"Source geometries: {len(group.sources)}",
            f"Selected geometries: {len(group.selected)}",
            f"Algorithm: {policy.selection.algorithm}",
            f"Connectivity policy: {policy.connectivity_policy or 'off'}",
        ],
        group_by_stoichiometry=True,
    )
    if snapshot_directory:
        _write_growth_selection_diagnostics(
            group.diagnostics,
            os.path.join(snapshot_directory, "selection_diagnostics.json"),
        )



def _supports_two_layer_optimization(qc_params):
    """Return True when the backend supports staged loose/normal optimization."""
    if "_two_layer_optimization" in qc_params:
        return bool(qc_params.get("_two_layer_optimization"))
    return qc_params.get("software") in {"xtb", "xtb-aiqm1", "xtb-aimnet2"}


def random_permutation(iterable, r=None, *, rng=None):
    """Random selection from itertools.permutations(iterable, r)."""
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple((random if rng is None else rng).sample(pool, r))


def read_old_path():
    """Read restart pathway markers from a prior log, if one exists."""
    if not os.path.exists("pyar.log"):
        return None

    path = []
    with open("pyar.log") as log_file:
        lines = log_file.readlines()
    needed_elements = 3 if any("DEBUG" in line for line in lines) else 2
    for line in lines:
        split_line = line.split(":")
        if len(split_line) == needed_elements and split_line[-2].strip().isnumeric():
            path.append(split_line[-1].rstrip())
    return path or None


def old_path_to_new_path(monomers_to_be_added, old_path):
    monomers_by_name = {}
    for monomer in monomers_to_be_added:
        monomers_by_name.setdefault(monomer.name, monomer)
    complete_pathways = []
    for each in old_path:
        label = str(each).strip()
        try:
            complete_pathways.append(tuple(monomers_by_name[name] for name in label))
        except KeyError as exc:
            raise ValueError(
                f"Legacy aggregation path {label!r} contains unknown fragment {exc.args[0]!r}."
            ) from exc
    return complete_pathways


def select_pathways(monomers_to_be_added, number_of_pathways, *, seed=0):
    """Return up to ``number_of_pathways`` distinct addition orders, reproducibly.

    A seeded generator and first-seen ordering make fresh runs repeatable; the
    attempt budget lets duplicate draws be retried without looping forever when
    fewer distinct orders exist than were requested.
    """
    rng = random.Random(seed)
    complete_pathways = []
    seen = set()
    for _ in range(max(0, int(number_of_pathways)) * 20):
        if len(complete_pathways) >= number_of_pathways:
            break
        new_permutation = random_permutation(monomers_to_be_added, rng=rng)
        key = tuple(monomer.name for monomer in new_permutation)
        if key not in seen:
            seen.add(key)
            complete_pathways.append(new_permutation)
    return complete_pathways


def read_orientations(molecule_id, noo, *, directory=None):
    directory = Path(directory or Path.cwd()).resolve()
    orientations = []
    for i in range(noo):
        xyz_file = f"trial_{i:03d}_{molecule_id}.xyz"
        new_orientation = Molecule.from_xyz(str(directory / xyz_file))
        new_orientation.name = f"{i:03d}_{molecule_id}"
        orientations.append(new_orientation)
    return orientations


def _raise_if_growth_stopped(directory):
    if check_stop_signal(directory=directory):
        aggregator_logger.info("Function: add_one")
        raise GrowthStopped(f"Growth stopped in {directory}")


def _resolve_growth_connectivity(policy, seeds, monomer):
    connectivity = (
        resolve_connectivity_policy_for_growth(seeds[0], monomer, policy.qc_params)
        if policy.connectivity_policy is None and seeds
        else policy.connectivity_policy or "off"
    )
    connectivity = connectivity.lower()
    if connectivity == "auto":
        connectivity = (
            resolve_connectivity_policy_for_growth(seeds[0], monomer, policy.qc_params)
            if seeds
            else "off"
        )
    if connectivity not in _CONNECTIVITY_POLICIES:
        raise ValueError(
            f"Unknown connectivity policy: {connectivity!r}. Expected one of 'off', 'prefer', 'strict', or 'auto'."
        )
    return connectivity


def _optimization_stages(qc_params):
    staged = _supports_two_layer_optimization(qc_params)
    loose = dict(qc_params)
    if staged:
        loose["opt_threshold"] = "loose"
    refine = dict(qc_params)
    refine["opt_threshold"] = "normal"
    return OptimizationStages(loose, refine, staged)


def _prepare_growth_run(aggregate_id, seeds, monomer, policy, directory):
    directory = Path(directory or Path.cwd()).resolve()
    _raise_if_growth_stopped(directory)
    aggregator_logger.info(f"  Seed pool for {aggregate_id}: {len(seeds)} molecules")
    connectivity = _resolve_growth_connectivity(policy, seeds, monomer)
    run = GrowthRun(
        aggregate_id,
        seeds,
        monomer,
        policy,
        directory,
        connectivity,
        _optimization_stages(policy.qc_params),
    )
    return run


def add_one(aggregate_id, seeds, monomer, policy: GrowthPolicy, *, directory=None):
    """Orchestrate one addition using owned settings and stage results."""
    run = _prepare_growth_run(aggregate_id, seeds, monomer, policy, directory)
    if not policy.qc_params.get("software"):
        return _grow_geometry_only(run)
    candidates = _optimize_growth_seeds(run)
    resumed = _resume_growth_selection(run, candidates)
    selected = _select_growth_seed_pool(candidates, resumed, run)
    return _refine_and_select_growth(run, selected)


def _prepare_seed_run(run, index, molecule):
    _raise_if_growth_stopped(run.directory)
    aggregator_logger.info("   Processing seed %03d/%03d", index + 1, len(run.seeds))
    seed = SeedRun(run, index, molecule)
    if not seed.directory.exists():
        file_manager.make_directories(seed.directory)
    molecule.mol_to_xyz(str(seed.directory / "seed.xyz"))
    run.monomer.mol_to_xyz(str(seed.directory / "monomer.xyz"))
    return seed


def _optimize_growth_seeds(run):
    aggregator_logger.info(
        "Optimization layers: %s",
        "two-stage (loose then normal)" if run.stages.staged else "single-stage",
    )
    accepted = []
    for index, molecule in enumerate(run.seeds):
        seed = _prepare_seed_run(run, index, molecule)
        accepted.extend(_retry_cycle_limited(seed, list(generate_orientations(seed))))
    return accepted


def _optimize_orientation_round(seed, pending, round_index):
    run = seed.run
    statuses = _optimize_growth_batch(
        pending,
        run.stages.loose,
        "growth_preselection",
        "backend_failures.json",
        directory=seed.directory,
    )
    accepted = [
        molecule for molecule, status in zip(pending, statuses) if is_success(status)
    ]
    retry = [
        molecule
        for molecule, status in zip(pending, statuses)
        if is_cycle_exceeded(status)
        and (run.connectivity_policy == "off" or not trial_generation.broken(molecule))
    ]
    _write_growth_structure_analysis(
        [seed.molecule, run.monomer],
        retry,
        str(seed.directory / f"structural_analysis_round_{round_index + 1:02d}.json"),
    )
    return SeedOptimization(accepted, clustering.remove_similar(retry), statuses)


def _retry_cycle_limited(seed, pending):
    accepted = []
    statuses = [False for _ in pending]
    for round_index in range(10):
        if not pending:
            aggregator_logger.info("    All trial molecules processed for this seed")
            break
        aggregator_logger.info(
            f"    Optimization round {round_index + 1:d}: pending={len(pending):d}"
        )
        outcome = _optimize_orientation_round(seed, pending, round_index)
        accepted.extend(outcome.accepted)
        pending, statuses = outcome.retry, outcome.statuses
    else:
        # Retain the original terminal logging alignment; changing it is not a refactor.
        aggregator_logger.info("    Molecules still unconverged after 10 rounds:")
        for molecule, status in zip(pending, statuses):
            if is_cycle_exceeded(status) and (
                seed.run.connectivity_policy == "off"
                or not trial_generation.broken(molecule)
            ):
                aggregator_logger.info("      %s", molecule.name)
    return accepted


def _resume_growth_selection(run, candidates):
    selected = run.directory / "selected"
    resumed = []
    if selected.exists():
        # The established helper removes adopted members from the fresh pool.
        resumed = check_for_the_finished_jobs_on_restart(candidates, str(run.directory))
        if resumed:
            aggregator_logger.info(
                "Found %d already refined selected geometries from restart.",
                len(resumed),
            )
    else:
        file_manager.make_directories(selected)
    aggregator_logger.info("  Selecting optimized pool")
    if run.seeds:
        _write_growth_structure_analysis(
            [run.seeds[0], run.monomer],
            candidates,
            str(run.directory / "structural_analysis_preselection.json"),
        )
    return resumed


def _refine_growth_seed(run, molecule):
    original = copy.deepcopy(molecule)
    directory = run.directory / "selected"
    status = _optimize_growth_batch(
        [molecule],
        run.stages.refine,
        "growth_refinement",
        "backend_failures.json",
        directory=directory,
    )[0]
    if is_success(status):
        shutil.copy(
            str(directory / f"job_{molecule.name}/result_{molecule.name}.xyz"),
            str(directory / "."),
        )
        return molecule, None
    return None, original


def _refine_and_select_growth(run, pool):
    if not run.stages.staged:
        return _publish_single_stage_growth(run, pool)
    aggregator_logger.info("Refining selected geometries with normal threshold")
    refined, loose = [], []
    for molecule in pool:
        accepted, fallback = _refine_growth_seed(run, molecule)
        if accepted is not None:
            refined.append(accepted)
        else:
            loose.append(fallback)
    if refined:
        aggregator_logger.info("Selection result: %d refined molecules", len(refined))
        selected, stage = refined, "Stage 2: normal-threshold refinement"
    else:
        aggregator_logger.info(
            "Selection result: no refined molecules, returning loose set (%d)",
            len(loose),
        )
        selected, stage = (
            loose,
            "Stage 2: normal-threshold refinement failed for all selected seeds",
        )
    _snapshot_selected_geometries(
        selected,
        output_root=str(run.directory / "selected"),
        summary_lines=[
            "Selection mode: two-stage",
            "Stage 1: loose preselection",
            stage,
            f"Backend: {run.policy.qc_params.get('software')}",
            f"Selected geometries: {len(selected)}",
            f"Connectivity policy: {run.connectivity_policy}",
        ],
        group_by_stoichiometry=False,
    )
    return selected


def _publish_single_stage_growth(run, selected):
    aggregator_logger.info(
        "Backend does not support staged optimization; returning selected seeds after a single pass."
    )
    _snapshot_selected_geometries(
        selected,
        output_root=str(run.directory / "selected"),
        summary_lines=[
            "Selection mode: single-stage",
            f"Backend: {run.policy.qc_params.get('software')}",
            f"Selected geometries: {len(selected)}",
            f"Connectivity policy: {run.connectivity_policy}",
        ],
        group_by_stoichiometry=False,
    )
    return selected


def _grow_geometry_only(run):
    # Structural diversity uses no fabricated energies; consumption order is fixed.
    orientations = []
    for index, molecule in enumerate(run.seeds):
        seed = _prepare_seed_run(run, index, molecule)
        orientations.extend(generate_orientations(seed))
    aggregator_logger.info(
        "Generated %d trial geometries without QM optimization", len(orientations)
    )
    if not run.policy.bound_geometry_only:
        return orientations
    from dataclasses import replace

    selected = select_geometry_only(
        orientations,
        replace(run.policy, connectivity_policy=run.connectivity_policy),
        directory=run.directory,
    )
    _snapshot_selected_geometries(
        selected,
        output_root=str(run.directory / "selected"),
        summary_lines=["Selection mode: geometry-diversity; no energies available"],
    )
    return selected


def check_for_the_finished_jobs_on_restart(list_of_optimized_molecules, cwd):
    selected = Path(cwd) / "selected"
    selected_names = set(path.name for path in selected.iterdir())
    result, remaining = [], []
    for molecule in list_of_optimized_molecules:
        if (
            f"job_{molecule.name}" in selected_names
            and f"result_{molecule.name}.xyz" in selected_names
        ):
            result.append(molecule)
        else:
            remaining.append(molecule)
    list_of_optimized_molecules[:] = remaining
    return result


def generate_orientations(seed: SeedRun):
    count = _seed_specific_orientation_count(
        seed.molecule, seed.run.policy.orientations
    )
    directory, molecule_id = seed.directory, seed.molecule_id
    aggregator_logger.debug("Making orientations")
    if not all(
        (directory / f"trial_{i:03d}_{molecule_id}.xyz").exists() for i in range(count)
    ):
        yield from trial_generation.create_trial_geometries(
            molecule_id,
            seed.molecule,
            seed.run.monomer,
            count,
            seed.run.policy.site,
            directory=directory,
        )
        aggregator_logger.debug("Orientations are made.")
    else:
        yield from read_orientations(molecule_id, count, directory=directory)


def update_id(aid, the_monomer):
    """Increment one monomer count inside an aggregate identifier."""
    from collections import deque

    parts_of_aid = deque(aid.split("_"))
    prefix = parts_of_aid.popleft()
    d = {}
    new_id = prefix
    while parts_of_aid:
        cid = parts_of_aid.popleft()
        n = parts_of_aid.popleft()
        d[cid] = int(n)
        if cid == the_monomer:
            d[cid] += 1
        new_id += f"_{cid}_{d[cid]:03d}"

    return new_id


def check_stop_signal(*, directory=None):
    directory = Path(directory or Path.cwd()).resolve()
    if os.path.exists(str(directory / "stop")) or os.path.exists(
        str(directory / "STOP")
    ):
        aggregator_logger.info(f"Found stop file, in {str(directory)}")
        return 1


_FORMULA_SYMBOLS = set(atomic_data.atomic_number)


def _normalize_symbol(symbol):
    """Normalize an element token to standard capitalization."""
    if len(symbol) == 1:
        return symbol.upper()
    return symbol[0].upper() + symbol[1:].lower()


def _read_formula_symbol(formula, position):
    """Read one element symbol from a formula string."""
    char = formula[position]
    if not char.isalpha():
        raise ValueError(f"Invalid chemical formula: {formula!r}")

    next_char = formula[position + 1] if position + 1 < len(formula) else ""
    if char.isupper() and next_char.islower():
        symbol = _normalize_symbol(formula[position : position + 2])
        if symbol in _FORMULA_SYMBOLS:
            return symbol, position + 2
        raise ValueError(f"Unknown element symbol: {formula[position:position + 2]}")

    if char.islower() and next_char.islower():
        symbol = _normalize_symbol(formula[position : position + 2])
        if symbol in _FORMULA_SYMBOLS:
            return symbol, position + 2

    symbol = _normalize_symbol(char)
    if symbol in _FORMULA_SYMBOLS:
        return symbol, position + 1

    if next_char and next_char.islower():
        raise ValueError(f"Unknown element symbol: {formula[position:position + 2]}")
    raise ValueError(f"Unknown element symbol: {char}")


def _parse_formula(formula):
    """Expand a chemical formula string into a list of atom symbols."""
    formula = formula.strip()
    atoms_list = []
    position = 0
    while position < len(formula):
        if formula[position].isspace():
            position += 1
            continue

        symbol, position = _read_formula_symbol(formula, position)

        digit_start = position
        while position < len(formula) and formula[position].isdigit():
            position += 1
        count = int(formula[digit_start:position]) if position > digit_start else 1
        atoms_list.extend([symbol] * count)

    if not atoms_list:
        raise ValueError(f"Invalid chemical formula: {formula!r}")
    return atoms_list


def expand_formula_to_aggregate_inputs(formula):
    """Convert a formula into aggregate input symbols and multiplicities."""
    atoms_list = _parse_formula(formula)
    fragments = []
    aggregate_sizes = []
    for atom in atoms_list:
        if fragments and fragments[-1] == atom:
            aggregate_sizes[-1] += 1
        else:
            fragments.append(atom)
            aggregate_sizes.append(1)
    return fragments, aggregate_sizes


def _atom_radius(atom):
    """Return a packing radius for an element symbol."""
    radius = atomic_data.vdw_radius.get(atom)
    if radius is None or np.isnan(radius):
        radius = atomic_data.covalent_radius[atom]
    return float(radius)


def _estimate_formula_box_size(atoms_list, box_size=None):
    """Estimate a coordinate box large enough for an initial formula geometry."""
    if box_size is not None:
        if box_size <= 0:
            raise ValueError("box_size must be positive")
        return float(box_size)

    max_radius = max(_atom_radius(atom) for atom in atoms_list)
    volume_scale = len(atoms_list) ** (1.0 / 3.0)
    return max(2.0 * max_radius + 2.0, 2.5 * volume_scale * max_radius + 2.0)


def _pack_formula_coordinates(atoms_list, box_size, rng=None, minimum_distance_factor=0.85, max_attempts=500):
    """Pack atoms into a box while avoiding obvious overlaps."""
    rng = np.random.default_rng() if rng is None else rng
    half_box = box_size / 2.0
    coordinates = np.zeros((len(atoms_list), 3), dtype=float)
    placed_indices = []
    placement_order = sorted(
        range(len(atoms_list)),
        key=lambda index: (-_atom_radius(atoms_list[index]), index),
    )

    for index in placement_order:
        atom = atoms_list[index]
        attempts = 0
        while True:
            candidate = rng.uniform(-half_box, half_box, size=3)
            if all(
                np.linalg.norm(candidate - coordinates[other_index])
                >= minimum_distance_factor
                * (_atom_radius(atom) + _atom_radius(atoms_list[other_index]))
                for other_index in placed_indices
            ):
                coordinates[index] = candidate
                placed_indices.append(index)
                break

            attempts += 1
            if attempts >= max_attempts:
                box_size *= 1.15
                half_box = box_size / 2.0
                attempts = 0
                aggregator_logger.debug(
                    "Expanded formula packing box to %s for atom %s", box_size, atom
                )

    return coordinates


def generate_molecule_from_formula(formula, box_size=None, rng=None):
    """Generate a molecule object from a chemical formula."""
    atoms_list = _parse_formula(formula)
    box_size = _estimate_formula_box_size(atoms_list, box_size=box_size)
    coordinates = _pack_formula_coordinates(atoms_list, box_size, rng=rng)
    molecule = Molecule(atoms_list, coordinates, name=formula, title=formula)
    molecule.connectivity_policy_hint = "prefer"
    molecule.source_kind = "formula"
    return molecule


def select_geometry_only(molecules, policy: GrowthPolicy, *, directory=None):
    """Bound an energy-free pool using graph-first pruning and max-min diversity."""
    from pyar.selection.deduplication import remove_similar, _prefer_connected_structures
    from pyar.selection.features import compute_feature_matrix, standardize_features
    from pyar.selection.diversity import _max_min_diversity_select
    from pyar.selection.structural_distances import compute_distance_matrix

    valid = [m for m in molecules if np.isfinite(m.coordinates).all()]
    unique = remove_similar(valid, energy_ranked=False)
    if (policy.connectivity_policy or "off") != "off":
        unique = _prefer_connected_structures(
            unique, policy=policy.connectivity_policy or "off"
        )
    diagnostics = {"selection_mode": "geometry-diversity", "energy_ranked": False,
                   "input_count": len(molecules), "unique_count": len(unique)}
    if len(unique) > policy.maximum_number_of_seeds:
        features = compute_feature_matrix(
            unique,
            feature=policy.selection.feature,
            system_type=policy.selection.system_type,
        )
        values = standardize_features(features.values)
        distances = compute_distance_matrix(
            unique, policy.selection.distance, feature_values=values
        )
        unique = _max_min_diversity_select(
            values,
            unique,
            policy.maximum_number_of_seeds,
            distance_metric=policy.selection.distance,
            distance_matrix=distances.values,
            energy_ranked=False,
        )
        diagnostics.update(distances.to_dict())
        diagnostics.update(feature_used=features.name, feature_fallbacks=list(features.fallbacks))
    diagnostics["selected_count"] = len(unique)
    _write_growth_selection_diagnostics(diagnostics, str(Path(directory or Path.cwd()) / "selection_diagnostics.json"))
    return unique
