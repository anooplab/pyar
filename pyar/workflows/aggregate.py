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
from pyar.growth.policy import GrowthPolicy, FinalSelectionRun
from dataclasses import dataclass, replace

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




@dataclass(frozen=True)
class AggregateRun:
    request: AggregateRequest
    root: Path
    layout: object
    state: object
    monomers: list
    initial_id: str
    policy: GrowthPolicy

    @property
    def directory(self):
        return self.root / 'aggregates'


def _open_aggregate_run(root, request, payload):
    from pyar.aggregation.layout_v2 import AggregateRunStateV2
    from pyar.persistence import RunLayout

    directory = root / 'aggregates'
    layout = RunLayout(directory, 'aggregate')
    state = AggregateRunStateV2.load(root, payload) if layout.marker.exists() else None
    molecules = list(request.molecules)
    number_of_orientations = _resolve_orientation_count(request.orientations)
    sampling = sampling_configuration(number_of_orientations=number_of_orientations, use_angles=None)
    final_connectivity_policy = (resolve_connectivity_policy_for_aggregate(molecules)
                                 if request.connectivity_policy == 'auto' else request.connectivity_policy)
    policy = replace(
        GrowthPolicy.from_aggregate_request(request, number_of_orientations),
        connectivity_policy=final_connectivity_policy,
    )
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
    return AggregateRun(
        request, root, layout, state, monomers_to_be_added, initial_id, policy
    )


@dataclass(frozen=True)
class AggregateStep:
    pathway: int
    index: int
    aggregate_id: str
    monomer: object
    seeds: list | None


def _run_aggregate_step(run, step):
    from pyar.aggregation.layout_v2 import verify_refs

    if run.state.selection_path(step.pathway, step.index).exists():
        if str(step.index) not in run.state.data["steps"].get(str(step.pathway), {}):
            return run.state.adopt_step(step.pathway, step.index)
        return verify_refs(
            run.state.selection(step.pathway, step.index)["selected"], run.directory
        )
    if step.index == 0:
        seeds = [step.monomer]
        run.state.commit_step(step.pathway, step.index, seeds, None, None)
        return seeds
    source_key = "{}_{:03d}".format(
        step.aggregate_id, run.request.first_pathway + step.pathway
    )
    native = run.layout.path(
        f"pathways/pathway_{step.pathway:03d}/steps/step_{step.index:03d}/attempts/attempt_000/execution/{source_key}"
    )
    native.mkdir(parents=True, exist_ok=True)
    selected = add_one(
        step.aggregate_id, step.seeds, step.monomer, run.policy, directory=native
    )
    run.state.commit_step(step.pathway, step.index, selected, native, source_key)
    return selected


def _run_aggregate_pathway(run: AggregateRun, index):
    label = run.state.data["pathways"][index]
    pathway = _pathway_from_label(run.monomers, label)
    aggregate_id, seeds = run.initial_id, None
    for position, monomer in enumerate(pathway):
        aggregate_id = update_id(aggregate_id, monomer.name)
        seeds = _run_aggregate_step(
            run, AggregateStep(index, position, aggregate_id, monomer, seeds)
        )
        if not seeds:
            break
    run.state.complete_pathway(index, label)


def _commit_final_aggregation(run: AggregateRun):
    if run.layout.path("final_selection.json").exists():
        run.state.adopt_final()
    else:
        execution = run.layout.path("selection/attempts/attempt_000/execution")
        execution.mkdir(parents=True, exist_ok=True)
        selected = _finalize_selected_geometries(
            FinalSelectionRun(
                execution,
                run.state.native_sources(),
                bool(run.request.backend_parameters.get("software")),
            ),
            run.policy,
        )
        run.state.finish(selected)

def _publish_aggregate_result(run: AggregateRun):
    from pyar.persistence import write_json

    result = AggregateResult(
        workflow="aggregate",
        status=run.state.data["status"],
        run_directory=str(run.directory),
        state_path=str(run.state.state_file),
        selected_paths=tuple(run.state.data["final_selected_results"]),
        metadata=dict(
            pathways=tuple(run.state.data["pathways"]),
            completed_pathways=tuple(
                (item["label"] for item in run.state.data["completed_pathways"])
            ),
            sampling=run.state.data["sampling"],
            connectivity_policy=run.policy.connectivity_policy,
            layout_version=2,
        ),
    )
    write_json(run.layout.path("summary.json"), result.to_dict())
    return result

def run_aggregate(request: AggregateRequest, *, recover_stale=False):
    """Execute a resolved aggregation request with existing pathway/restart contracts."""
    from pyar.persistence import detect_version, RunLock

    root = Path.cwd().resolve()
    if check_stop_signal():
        return stopped_workflow_result("aggregate", root)
    payload = request.to_state_dict()
    directory = root / 'aggregates'
    detect_version(directory, 'aggregate')
    with RunLock(directory, recover_stale=recover_stale):
        detect_version(directory, 'aggregate')
        run = _open_aggregate_run(root, request, payload)
        try:
            for index in range(
                len(run.state.data["completed_pathways"]),
                len(run.state.data["pathways"]),
            ):
                _run_aggregate_pathway(run, index)
            _commit_final_aggregation(run)
            return _publish_aggregate_result(run)
        except GrowthStopped:
            return stopped_workflow_result("aggregate", root)


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
    """Build formula fragments and execute one resolved aggregation request."""
    molecules = [generate_molecule_from_formula(formula) for formula in formulas]
    if check_stop_signal():
        return stopped_workflow_result('aggregate', os.getcwd())
    options = dict(locals())
    options.pop("formulas")
    return run_aggregate(AggregateRequest.from_options(**options))


def aggregate(
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
    *,
    recover_stale=False,
):
    """Keep the keyword API and convert scientific options only once."""
    if check_stop_signal():
        return stopped_workflow_result('aggregate', os.getcwd())
    options = dict(locals())
    options.pop("recover_stale")
    return run_aggregate(
        AggregateRequest.from_options(**options), recover_stale=recover_stale
    )
