"""Aggregate workflow orchestration for PyAR.

This module owns the aggregation and cluster-generation workflow used by
``pyar-cli -a``. It is
responsible for:

* validating the aggregate request and restart state
* creating or resuming the ``aggregates/`` run directory
* selecting pathways and orientation sets
* generating trial candidates and running backend optimization
* selecting survivors from pathway-level and final clustering stages
* persisting ``aggregates/state.json`` plus selected-geometry snapshots
* returning a structured :class:`~pyar.workflow_results.AggregateResult`
* running the backend optimization steps for each pathway
* collecting the selected geometries and workflow result metadata

The public entry points are :func:`aggregate` and
:func:`aggregate_from_formulas`.
"""

from __future__ import annotations

import copy
import json
import os
import string
from collections import OrderedDict
from pathlib import Path
from pyar.structure_comparison.coordinate_graph import analyze_coordinate_structure

from pyar.aggregation import AggregateRequest
from pyar.state.aggregate import AggregateRunState, AggregateStateError
from pyar import file_manager
from pyar.workflow_results import AggregateResult
from pyar.growth.service import (
    add_one,
    aggregator_logger,
    check_stop_signal,
    expand_formula_to_aggregate_inputs,
    generate_molecule_from_formula,
    old_path_to_new_path,
    read_old_path,
    select_pathways,
    resolve_connectivity_policy_for_aggregate,
    update_id,
    _finalize_selected_geometries,
    _format_path,
    sampling_configuration,
    _resolve_orientation_count,
    workflow_run_directory,
    workflow_state_path,
    stopped_workflow_result,
    _working_directory,
)

__all__ = [
    "aggregate",
    "aggregate_from_formulas",
    "expand_formula_to_aggregate_inputs",
    "generate_molecule_from_formula",
]


def _pathway_from_label(monomers_to_be_added, label):
    """Restore one persisted pathway label to its input molecule objects.

    Pathways are persisted as labels so they can survive restarts and legacy
    checkpoint migration. This helper maps a stored label back to the concrete
    molecule objects used in the current run.
    """
    by_name = {}
    for molecule in monomers_to_be_added:
        by_name.setdefault(molecule.name, molecule)
    try:
        return tuple(by_name[name] for name in label)
    except KeyError as exc:
        raise AggregateStateError(
            f"Saved aggregation pathway {label!r} contains unknown fragment {exc.args[0]!r}."
        ) from exc


def _selected_result_paths(pattern):
    """Return result paths relative to the aggregate workflow directory."""
    return sorted(str(path) for path in Path(".").glob(pattern))


def _aggregate_v1(
    molecules,
    aggregate_sizes,
    hm_orientations,
    qc_params,
    maximum_number_of_seeds,
    first_pathway,
    number_of_pathways,
    site,
    connectivity_policy="auto",
    selection_feature="auto",
    selection_algorithm="auto",
    selection_distance="euclidean",
    selection_system_type="auto",
):
    """Run an aggregate or cluster-generation workflow.

    The function creates or resumes the aggregation state, iterates over the
    selected pathways, runs the per-orientation optimization for each path,
    and returns an :class:`~pyar.workflow_results.AggregateResult` describing
    the final run state.
    """
    if check_stop_signal():
        aggregator_logger.info("Function: aggregate")
        return stopped_workflow_result("aggregate", os.getcwd())

    aggregate_request = AggregateRequest.from_options(
        molecules,
        aggregate_sizes,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        first_pathway,
        number_of_pathways,
        site,
        connectivity_policy,
        selection_feature,
        selection_algorithm,
        selection_distance,
        selection_system_type,
    )
    molecules = list(aggregate_request.molecules)
    aggregate_sizes = list(aggregate_request.aggregate_sizes)
    qc_params = dict(aggregate_request.backend_parameters)
    maximum_number_of_seeds = aggregate_request.maximum_number_of_seeds
    first_pathway = aggregate_request.first_pathway
    number_of_pathways = aggregate_request.number_of_pathways
    site = None if aggregate_request.site is None else list(aggregate_request.site)
    connectivity_policy = aggregate_request.connectivity_policy
    selection_feature = aggregate_request.selection_feature
    selection_algorithm = aggregate_request.selection_algorithm
    selection_distance = aggregate_request.selection_distance
    selection_system_type = aggregate_request.selection_system_type

    number_of_orientations = _resolve_orientation_count(hm_orientations)
    sampling = sampling_configuration(
        number_of_orientations=number_of_orientations,
        use_angles=None,
    )

    root_directory = os.getcwd()
    parent_folder = "aggregates"
    request = aggregate_request.to_state_dict()
    run_state = AggregateRunState.load(root_directory, request)
    old_path = read_old_path()
    parent_folder_exists = os.path.isdir(parent_folder)
    restart = run_state is not None or bool(old_path and parent_folder_exists)

    if not parent_folder_exists:
        file_manager.make_directories(parent_folder)
        if old_path:
            aggregator_logger.info(
                "Ignoring restart markers because %s was missing; starting a new aggregation.",
                parent_folder,
            )
    elif run_state is None and not old_path and any(Path(parent_folder).iterdir()):
        raise AggregateStateError(
            "Existing aggregates directory has no resumable state; "
            "preserve it and start in a new directory, or remove it explicitly."
        )

    with _working_directory(parent_folder):
        starting_directory = os.getcwd()

        # Capture what the workflow was given before aggregation changes the
        # structures. This is a coordinate-only diagnostic, not identity
        # assignment and not an input to candidate selection.
        structural_analysis_directory = Path("structural_analysis")
        structural_analysis_directory.mkdir(exist_ok=True)
        try:
            input_analysis = {
                "stage": "workflow-input",
                "workflow": "aggregate",
                "method": "coordinate-only-adjacency",
                "structures": [
                    analyze_coordinate_structure(molecule) for molecule in molecules
                ],
                "limitations": [
                    "XYZ coordinates do not provide bond orders, charge, spin, or fragment intent.",
                    "Adjacency is geometric and should not be interpreted as definitive chemical identity.",
                ],
            }
        except (KeyError, ValueError, TypeError) as exc:
            # Diagnostic analysis is informative only and must not gate a run.
            aggregator_logger.warning("Input structural analysis skipped: %s", exc)
            input_analysis = {
                "stage": "workflow-input",
                "workflow": "aggregate",
                "method": "coordinate-only-adjacency",
                "status": "skipped",
                "reason": f"{type(exc).__name__}: {exc}",
                "structures": [],
            }
        with (structural_analysis_directory / "input.json").open("w", encoding="utf-8") as stream:
            json.dump(input_analysis, stream, indent=2)
            stream.write("\n")

        if restart:
            aggregator_logger.info(f"Restarting aggregation in {starting_directory}")
        else:
            aggregator_logger.info(f"Starting aggregation in {starting_directory}")
        aggregator_logger.info(
            "Aggregate config: fragments=%d sizes=%s orientations=%s software=%s max_seeds=%s pathways=%s first_pathway=%s",
            len(molecules),
            aggregate_sizes,
            hm_orientations,
            qc_params.get("software"),
            maximum_number_of_seeds,
            number_of_pathways,
            first_pathway,
        )
        aggregator_logger.debug("Aggregate site=%s restart=%s", site, restart)

        seed_names = string.ascii_lowercase
        ag_id = "ag"
        requested_connectivity_policy = connectivity_policy
        if requested_connectivity_policy == "auto":
            final_connectivity_policy = resolve_connectivity_policy_for_aggregate(molecules)
        else:
            final_connectivity_policy = requested_connectivity_policy
        if final_connectivity_policy == "prefer":
            connectivity_reason = "atomic/formula growth"
        elif final_connectivity_policy == "strict":
            connectivity_reason = "explicit CLI override"
        elif requested_connectivity_policy == "auto":
            connectivity_reason = "molecular aggregate / noncovalent complex"
        else:
            connectivity_reason = "explicit CLI override"
        aggregator_logger.info(
            "Connectivity policy: %s (%s)",
            final_connectivity_policy,
            connectivity_reason,
        )

        monomers_to_be_added = []
        for seed_molecule, seed_name, size_of_this_seed in zip(
            molecules, seed_names, aggregate_sizes
        ):
            seed_molecule.name = seed_name
            monomers_to_be_added.extend(seed_molecule for _ in range(size_of_this_seed))
            ag_id += f"_{seed_name}_000"

        aggregator_logger.info("Aggregate identifier: %s", ag_id)
        aggregator_logger.info("Aggregate output root: %s/%s", starting_directory, ag_id)

        if run_state is not None:
            aggregator_logger.info(
                "Aggregation state detected: resuming %d remaining pathways",
                len(run_state.remaining_pathway_labels()),
            )
        elif old_path and parent_folder_exists:
            pathways_to_calculate = old_path_to_new_path(monomers_to_be_added, old_path)
            run_state = AggregateRunState.create(
                root_directory,
                request,
                [_format_path(path) for path in pathways_to_calculate],
                legacy_import="pyar.log",
                sampling=sampling,
            )
            aggregator_logger.warning(
                "Imported legacy pyar.log pathway markers into aggregates/state.json."
            )
        else:
            if len(molecules) == 1:
                pathways_to_calculate = [monomers_to_be_added]
            else:
                pathways_to_calculate = select_pathways(
                    monomers_to_be_added,
                    number_of_pathways,
                )
            aggregator_logger.info("Selected build paths:")
            for i, path in enumerate(pathways_to_calculate):
                aggregator_logger.info("  %03d: %s", i, _format_path(path))
            run_state = AggregateRunState.create(
                root_directory,
                request,
                [_format_path(path) for path in pathways_to_calculate],
                sampling=sampling,
            )

        pathways_to_calculate = [
            (label, _pathway_from_label(monomers_to_be_added, label))
            for label in run_state.remaining_pathway_labels()
        ]

        seed_storage = OrderedDict()
        initial_storage = copy.deepcopy(seed_storage)
        initial_aggregate_id = ag_id
        outside_counter = first_pathway + run_state.completed_pathway_count()

        for label, pathway in pathways_to_calculate:
            aggregator_logger.info("Path start: %s", label)
            for this_monomer in pathway:
                if len(seed_storage) < 1:
                    ag_id = update_id(ag_id, this_monomer.name)
                    seed_storage[ag_id] = [this_monomer]
                    continue
                this_seed = seed_storage[ag_id]
                ag_id = update_id(ag_id, this_monomer.name)
                ag_home = "{}_{:03d}".format(ag_id, outside_counter)
                if not os.path.exists(ag_home):
                    file_manager.make_directories(ag_home)
                with _working_directory(ag_home):
                    aggregator_logger.info("Aggregate cycle path: %s", ag_home)
                    seed_storage[ag_id] = add_one(
                        ag_id,
                        this_seed,
                        this_monomer,
                        number_of_orientations,
                        qc_params,
                        maximum_number_of_seeds,
                        site,
                        connectivity_policy=final_connectivity_policy,
                        selection_feature=selection_feature,
                        selection_algorithm=selection_algorithm,
                        selection_distance=selection_distance,
                        selection_system_type=selection_system_type,
                    )
                if len(seed_storage[ag_id]) == 0:
                    aggregator_logger.info(
                        f"No molecules were found from {ag_id} to continue this pathway."
                    )
                    aggregator_logger.info("Stopping this path.")
                    break
                seed_storage.popitem(last=False)
            run_state.complete_pathway(
                label,
                _selected_result_paths(f"ag_*_{outside_counter:03d}/selected/result_*.xyz"),
            )
            outside_counter += 1
            seed_storage = copy.copy(initial_storage)
            ag_id = initial_aggregate_id

        final_selected = _finalize_selected_geometries(
            aggregate_root=".",
            energy_ranked=bool(qc_params.get("software")),
            maximum_number_of_seeds=maximum_number_of_seeds,
            feature=selection_feature,
            algorithm=selection_algorithm,
            distance_metric=selection_distance,
            connectivity_policy=final_connectivity_policy,
            system_type=selection_system_type,
        )
        if final_selected:
            aggregator_logger.info(
                "Final cross-path selected geometries written: %d",
                len(final_selected),
            )
        run_state.finish(_selected_result_paths("selected/**/result_*.xyz"))
        result = AggregateResult(
            workflow="aggregate",
            status=run_state.data["status"],
            run_directory=workflow_run_directory(root_directory, "aggregates"),
            state_path=workflow_state_path(workflow_run_directory(root_directory, "aggregates")),
            selected_paths=tuple(run_state.data.get("final_selected_results", [])),
            metadata={
                "pathways": tuple(run_state.data.get("pathways", [])),
                "completed_pathways": tuple(
                    record["label"] for record in run_state.data.get("completed_pathways", [])
                ),
                "legacy_import": run_state.data.get("legacy_import"),
                "sampling": run_state.data.get("sampling", sampling),
                "connectivity_policy": final_connectivity_policy,
            },
        )

        if hm_orientations == "auto" and number_of_orientations <= 256:
            number_of_orientations += 8
        aggregator_logger.info("Aggregation workflow completed.")
        return result


def aggregate_from_formulas(
    formulas,
    aggregate_sizes,
    hm_orientations,
    qc_params,
    maximum_number_of_seeds,
    first_pathway,
    number_of_pathways,
    site,
    connectivity_policy="auto",
    selection_feature="auto",
    selection_algorithm="auto",
    selection_distance="euclidean",
    selection_system_type="auto",
):
    """Generate initial molecules from formulas and run the aggregate workflow.

    This is the convenience entry point used by ``--formula`` aggregation
    runs. It converts the provided formulas into molecule objects and then
    delegates to :func:`aggregate`.
    """
    molecules = [generate_molecule_from_formula(formula) for formula in formulas]
    return aggregate(
        molecules,
        aggregate_sizes,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        first_pathway,
        number_of_pathways,
        site,
        connectivity_policy=connectivity_policy,
        selection_feature=selection_feature,
        selection_algorithm=selection_algorithm,
        selection_distance=selection_distance,
        selection_system_type=selection_system_type,
    )


def aggregate(molecules, aggregate_sizes, hm_orientations, qc_params, maximum_number_of_seeds,
              first_pathway, number_of_pathways, site, connectivity_policy="auto",
              selection_feature="auto", selection_algorithm="auto", selection_distance="euclidean",
              selection_system_type="auto", *, recover_stale=False):
    """Use explicit v2 records for fresh runs, retaining v1 and legacy imports."""
    from pyar.persistence import detect_version, RunLock
    if check_stop_signal():
        return stopped_workflow_result('aggregate', os.getcwd())
    root = Path.cwd().resolve()
    args = (molecules, aggregate_sizes, hm_orientations, qc_params, maximum_number_of_seeds,
            first_pathway, number_of_pathways, site, connectivity_policy, selection_feature,
            selection_algorithm, selection_distance, selection_system_type)
    request = AggregateRequest.from_options(*args)
    payload = request.to_state_dict()
    directory = root / 'aggregates'
    def legacy(path):
        return AggregateRunState.load(root, payload) is not None
    # Established pyar.log imports belong exclusively to the v1 implementation.
    if not (directory / 'layout.json').exists() and directory.is_dir() and read_old_path():
        return _aggregate_v1(*args)
    version = detect_version(directory, 'aggregate', validate_v1=legacy)
    if version == 1:
        return _aggregate_v1(*args)
    with RunLock(directory, recover_stale=recover_stale):
        version = detect_version(directory, 'aggregate', validate_v1=legacy)
        if version == 1:
            return _aggregate_v1(*args)
        return _aggregate_v2(root, request, hm_orientations)


def _aggregate_v2(root, request, hm_orientations):
    from pyar.aggregation.layout_v2 import AggregateRunStateV2, verify_refs
    from pyar.persistence import RunLayout, write_json
    directory = root / 'aggregates'
    layout = RunLayout(directory, 'aggregate')
    payload = request.to_state_dict()
    state = AggregateRunStateV2.load(root, payload) if layout.marker.exists() else None
    molecules = list(request.molecules)
    number_of_orientations = _resolve_orientation_count(hm_orientations)
    sampling = sampling_configuration(number_of_orientations=number_of_orientations, use_angles=None)
    qc_params = dict(request.backend_parameters)
    site = None if request.site is None else list(request.site)
    final_connectivity_policy = (resolve_connectivity_policy_for_aggregate(molecules)
                                 if request.connectivity_policy == 'auto' else request.connectivity_policy)
    monomers_to_be_added = []
    initial_id = 'ag'
    for molecule, label, size in zip(molecules, string.ascii_lowercase, request.aggregate_sizes):
        molecule.name = label
        monomers_to_be_added.extend(molecule for _ in range(size))
        initial_id += f'_{label}_000'
    if state is None:
        pathways = ([monomers_to_be_added] if len(molecules) == 1 else
                    select_pathways(monomers_to_be_added, request.number_of_pathways))
        state = AggregateRunStateV2.create(root, payload, molecules,
                                          [_format_path(path) for path in pathways], sampling)
    for index in range(len(state.data['completed_pathways']), len(state.data['pathways'])):
        label = state.data['pathways'][index]
        pathway = _pathway_from_label(monomers_to_be_added, label)
        ag_id = initial_id
        seeds = None
        for step, monomer in enumerate(pathway):
            ag_id = update_id(ag_id, monomer.name)
            if state.selection_path(index, step).exists():
                if str(step) not in state.data['steps'].get(str(index), {}):
                    seeds = state.adopt_step(index, step)
                else:
                    seeds = verify_refs(state.selection(index, step)['selected'], directory)
                if not seeds:
                    break
                continue
            if step == 0:
                seeds = [monomer]
                state.commit_step(index, step, seeds, None, None)
                continue
            legacy_home = '{}_{:03d}'.format(ag_id, request.first_pathway + index)
            native = layout.path(f'pathways/pathway_{index:03d}/steps/step_{step:03d}/attempts/attempt_000/execution/{legacy_home}')
            native.mkdir(parents=True, exist_ok=True)
            with _working_directory(native):
                selected = add_one(ag_id, seeds, monomer, number_of_orientations,
                    qc_params, request.maximum_number_of_seeds, site,
                    connectivity_policy=final_connectivity_policy,
                    selection_feature=request.selection_feature, selection_algorithm=request.selection_algorithm,
                    selection_distance=request.selection_distance, selection_system_type=request.selection_system_type)
            state.commit_step(index, step, selected, native, legacy_home)
            seeds = selected
            if not seeds:
                break
        state.complete_pathway(index, label)
    if layout.path('final_selection.json').exists():
        state.adopt_final()
    else:
        execution = layout.path('selection/attempts/attempt_000/execution')
        execution.mkdir(parents=True, exist_ok=True)
        with _working_directory(execution):
            selected = _finalize_selected_geometries(aggregate_root=str(execution),
                maximum_number_of_seeds=request.maximum_number_of_seeds,
                algorithm=request.selection_algorithm, feature=request.selection_feature,
                distance_metric=request.selection_distance, system_type=request.selection_system_type,
                connectivity_policy=final_connectivity_policy,
                energy_ranked=bool(request.backend_parameters.get('software')),
                source_records=state.native_sources())
        state.finish(selected)
    result = AggregateResult(workflow='aggregate', status=state.data['status'],
        run_directory=str(directory), state_path=str(state.state_file),
        selected_paths=tuple(state.data['final_selected_results']),
        metadata=dict(pathways=tuple(state.data['pathways']),
                      completed_pathways=tuple(item['label'] for item in state.data['completed_pathways']),
                      legacy_import=None, sampling=state.data['sampling'],
                      connectivity_policy=final_connectivity_policy, layout_version=2))
    write_json(layout.path('summary.json'), result.to_dict())
    return result
