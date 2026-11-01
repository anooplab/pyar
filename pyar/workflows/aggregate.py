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

import os
import string
from pathlib import Path
from pyar.growth.errors import GrowthStopped

from pyar.aggregation import AggregateRequest
from pyar.state.aggregate import AggregateStateError
from pyar.workflow_results import AggregateResult
from pyar.growth.service import (
    add_one,
    check_stop_signal,
    expand_formula_to_aggregate_inputs,
    generate_molecule_from_formula,
    select_pathways,
    resolve_connectivity_policy_for_aggregate,
    update_id,
    _finalize_selected_geometries,
    _format_path,
    sampling_configuration,
    _resolve_orientation_count,
    stopped_workflow_result,
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
    """Create or resume supported immutable aggregation records."""
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
    detect_version(directory, 'aggregate')
    with RunLock(directory, recover_stale=recover_stale):
        detect_version(directory, 'aggregate')
        try:
            return _aggregate_v2(root, request, hm_orientations)
        except GrowthStopped:
            # Keep the pending checkpoint resumable; never publish a partial selection.
            return stopped_workflow_result("aggregate", root)


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
    for index in range(
        len(state.data["completed_pathways"]), len(state.data["pathways"])
    ):
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
            source_key = '{}_{:03d}'.format(ag_id, request.first_pathway + index)
            native = layout.path(f'pathways/pathway_{index:03d}/steps/step_{step:03d}/attempts/attempt_000/execution/{source_key}')
            native.mkdir(parents=True, exist_ok=True)
            selected = add_one(
                ag_id,
                seeds,
                monomer,
                number_of_orientations,
                qc_params,
                request.maximum_number_of_seeds,
                site,
                connectivity_policy=final_connectivity_policy,
                selection_feature=request.selection_feature,
                selection_algorithm=request.selection_algorithm,
                selection_distance=request.selection_distance,
                selection_system_type=request.selection_system_type,
                directory=native,
            )
            state.commit_step(index, step, selected, native, source_key)
            seeds = selected
            if not seeds:
                break
        state.complete_pathway(index, label)
    if layout.path("final_selection.json").exists():
        state.adopt_final()
    else:
        execution = layout.path('selection/attempts/attempt_000/execution')
        execution.mkdir(parents=True, exist_ok=True)
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
                      sampling=state.data['sampling'],
                      connectivity_policy=final_connectivity_policy, layout_version=2))
    write_json(layout.path('summary.json'), result.to_dict())
    return result
