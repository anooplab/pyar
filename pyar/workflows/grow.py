"""Fixed-addend growth with explicit request, step and publication ownership."""
from dataclasses import dataclass
from pathlib import Path

from pyar.growth.errors import GrowthStopped
from pyar.growth.policy import GrowthPolicy
from pyar.growth.request import GrowRequest
from pyar.growth.service import add_one, check_stop_signal
from pyar.workflow_results import GrowResult


@dataclass(frozen=True)
class GrowRun:
    request: GrowRequest
    state: object
    layout: object

    @property
    def directory(self):
        return self.state.directory


@dataclass(frozen=True)
class GrowStepResult:
    seeds: list
    status: str = 'running'


def grow(request: GrowRequest, *, output="grow", recover_stale=False):
    """Keep the supported growth interface; execute its resolved request."""
    return run_grow(request, output=output, recover_stale=recover_stale)


def run_grow(request: GrowRequest, *, output="grow", recover_stale=False):
    """Create/resume immutable growth steps using an already validated request."""
    from pyar.persistence import detect_version, RunLock

    directory = Path(output).resolve()
    payload = request.to_state_dict()
    detect_version(directory, "grow")
    directory.parent.mkdir(parents=True, exist_ok=True)
    with RunLock(directory, recover_stale=recover_stale):
        run = _open_grow_run(request, directory, payload)
        return _publish_grow_result(run, _advance_grow_steps(run))


def _open_grow_run(request, directory, payload):
    from pyar.growth.layout_v2 import GrowRunStateV2
    from pyar.persistence import detect_version, RunLayout

    version = detect_version(directory, 'grow')
    state = (GrowRunStateV2.create(directory, payload, request.seed, request.monomer)
             if version == 0 else GrowRunStateV2.load(directory, payload))
    return GrowRun(request, state, RunLayout(directory, "grow"))


def _advance_grow_steps(run):
    from pyar.growth.layout_v2 import verify_refs

    seeds = verify_refs(run.state.data["current_seeds"], run.directory)
    status = run.state.data["status"]
    if status not in {"running", "stopped"}:
        return status
    run.state.data["status"] = "running"
    for step in range(run.state.data["next_step"], run.request.count + 1):
        outcome = _run_grow_step(run, seeds, step)
        seeds = outcome.seeds
        if outcome.status != "running":
            return outcome.status
    return "completed"


def _growth_step_policy(run, seeds):
    from pyar.molecule_merge import combine_multiplicity

    request = run.request
    site = None if request.site is None else [[request.site[0]], [len(seeds[0]) + request.site[1]]]
    qc = dict(request.backend_parameters)
    qc.update(charge=seeds[0].charge + request.monomer.charge,
              multiplicity=combine_multiplicity(seeds[0].multiplicity, request.monomer.multiplicity),
              scftype='rhf' if seeds[0].scftype == request.monomer.scftype == 'rhf' else 'uhf')
    return GrowthPolicy.from_grow_request(request, qc, site)


def _stop_grow_step(run, seeds):
    run.state.data["status"] = "stopped"
    run.state.save()
    return GrowStepResult(seeds, "stopped")


def _run_grow_step(run, seeds, step):
    from pyar.growth.layout_v2 import verify_refs

    if run.state.selection_path(step).exists():
        run.state.accept_step(step)
        seeds = verify_refs(run.state.data["current_seeds"], run.directory)
        return GrowStepResult(seeds, "running" if seeds else "no_candidates")
    if check_stop_signal():
        return _stop_grow_step(run, seeds)
    execution = run.layout.path(f"steps/step_{step:03d}/attempts/attempt_000/execution")
    execution.mkdir(parents=True, exist_ok=True)
    try:
        selected = add_one(
            f"grow_{step:03d}",
            seeds,
            run.request.monomer,
            _growth_step_policy(run, seeds),
            directory=execution,
        )
    except GrowthStopped:
        return _stop_grow_step(run, seeds)
    if selected is None or selected is StopIteration:
        return _stop_grow_step(run, seeds)
    _commit_grow_step(run, step, selected, execution)
    return GrowStepResult(selected, "running" if selected else "no_candidates")


def _commit_grow_step(run, step, selected, execution):
    from pyar.growth.layout_v2 import durable_pool
    from pyar.persistence import read_json, write_json

    if len(selected) > run.request.maximum_number_of_seeds:
        raise RuntimeError('Shared growth engine exceeded the survivor budget')
    stage = run.layout.step(step)
    refs = durable_pool(selected, stage / "selected", run.directory)
    mode = (
        "energy-ranked"
        if run.request.backend_parameters.get("software")
        else "geometry-diversity"
    )
    source = execution / 'selected/selection_diagnostics.json'
    if not source.exists():
        source = execution / 'selection_diagnostics.json'
    diagnostics = (
        read_json(source)
        if source.exists()
        else dict(
            selection_mode=mode,
            selected_count=len(selected),
            maximum_number_of_seeds=run.request.maximum_number_of_seeds,
        )
    )
    write_json(stage / "selected/selection_diagnostics.json", diagnostics)
    run.state.commit_step(step, refs, mode)


def _grow_result_metadata(run, payload):
    return dict(
        requested_additions=run.request.count,
        completed_additions=len(run.state.data["completed_steps"]),
        orientations_per_seed=run.request.number_of_orientations,
        maximum_number_of_seeds=run.request.maximum_number_of_seeds,
        backend=run.request.backend_parameters.get("software"),
        connectivity_policy=run.request.connectivity_policy,
        selection_policy={k: payload[k] for k in payload if k.startswith("selection_")},
        stages=run.state.data["completed_steps"],
        sampling=payload["sampling"],
        layout_version=2,
    )


def _publish_grow_result(run, status):
    from pyar.persistence import write_json

    if status == "completed":
        run.state.finish()
    paths = (
        tuple(str(run.directory / ref["path"]) for ref in run.state.data["final_refs"])
        if status == "completed"
        else ()
    )
    result = GrowResult(
        workflow="grow",
        status=status,
        run_directory=str(run.directory),
        state_path=str(run.state.path),
        selected_paths=paths,
        metadata=_grow_result_metadata(run, run.request.to_state_dict()),
    )
    if status == "completed":
        write_json(run.directory / "summary.json", result.to_dict())
        write_json(
            run.layout.path("results/manifest.json"),
            dict(
                format="pyar-results-manifest",
                manifest_schema_version=1,
                workflow="grow",
                run_id=run.state.data["run_id"],
                artifacts=run.state.data["final_refs"],
            ),
        )
    return result
