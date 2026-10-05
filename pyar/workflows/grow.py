"""Fixed-seed, repeated-addend growth using the shared addition engine."""
from pathlib import Path
from pyar.growth.request import GrowRequest
from pyar.growth.reporting import snapshot_pool, restore_pool
from pyar.growth.service import add_one, _working_directory, check_stop_signal
from pyar.state.grow import GrowRunState, write_json
from pyar.workflow_results import GrowResult


def grow(request: GrowRequest, *, output="grow"):
    """Retain every completed growth stage and resume a compatible request."""
    directory = Path(output).resolve()
    payload = request.to_state_dict()
    state = GrowRunState.load(directory, payload)
    if state is None:
        directory.mkdir(parents=True, exist_ok=True)
        refs = snapshot_pool([request.seed], directory / "step_000" / "selected", directory)
        state = GrowRunState(directory, {"version": 1, "workflow": "grow", "request": payload,
                                        "status": "running", "next_step": 1,
                                        "completed_steps": [], "current_seeds": refs})
        state.save()
        write_json(directory / "request.json", payload)
    seeds = restore_pool(state.data["current_seeds"], directory)
    status = state.data["status"]
    if status == "running":
        for step in range(state.data["next_step"], request.count + 1):
            if check_stop_signal():
                status = "stopped"
                break
            stage = directory / f"step_{step:03d}"
            stage.mkdir(exist_ok=True)
            site = None if request.site is None else [[request.site[0]], [len(seeds[0]) + request.site[1]]]
            from pyar.molecule_merge import combine_multiplicity
            qc = dict(request.backend_parameters)
            qc.update(charge=seeds[0].charge + request.monomer.charge,
                      multiplicity=combine_multiplicity(seeds[0].multiplicity, request.monomer.multiplicity),
                      scftype="rhf" if seeds[0].scftype == request.monomer.scftype == "rhf" else "uhf")
            with _working_directory(stage):
                selected = add_one(f"grow_{step:03d}", seeds, request.monomer,
                                   request.number_of_orientations, qc,
                                   request.maximum_number_of_seeds, site,
                                   connectivity_policy=request.connectivity_policy,
                                   selection_feature=request.selection_feature,
                                   selection_algorithm=request.selection_algorithm,
                                   selection_distance=request.selection_distance,
                                   selection_system_type=request.selection_system_type)
            if selected is None or selected is StopIteration:
                status = "stopped"
                break
            if len(selected) > request.maximum_number_of_seeds:
                raise RuntimeError("Shared growth engine exceeded the survivor budget")
            refs = snapshot_pool(selected, stage / "selected", directory)
            diagnostics_path = stage / "selected" / "selection_diagnostics.json"
            if not diagnostics_path.exists():
                source = stage / "selection_diagnostics.json"
                if source.exists():
                    import json
                    diagnostics = json.loads(source.read_text())
                else:
                    diagnostics = {"selection_mode": "energy-ranked" if request.backend_parameters.get("software") else "geometry-diversity",
                                   "selected_count": len(selected), "maximum_number_of_seeds": request.maximum_number_of_seeds}
                write_json(diagnostics_path, diagnostics)
            state.data["completed_steps"].append({"step": step, "selected_count": len(selected),
                                                   "selected_paths": [r["path"] for r in refs],
                                                   "selection_mode": "energy-ranked" if request.backend_parameters.get("software") else "geometry-diversity"})
            state.data.update(current_seeds=refs, next_step=step + 1)
            seeds = selected
            if not seeds:
                status = "no_candidates"
                state.data["status"] = status
                state.save()
                break
            state.save()
        else:
            status = "completed"
    paths = ()
    if status == "completed":
        refs = snapshot_pool(seeds, directory / "final", directory)
        paths = tuple(str(directory / ref["path"]) for ref in refs)
        state.data.update(status=status, final_paths=list(paths))
        state.save()
    result = GrowResult(workflow="grow", status=status, run_directory=str(directory),
                        state_path=str(state.path), selected_paths=paths,
                        metadata={"requested_additions": request.count,
                                  "completed_additions": len(state.data["completed_steps"]),
                                  "orientations_per_seed": request.number_of_orientations,
                                  "maximum_number_of_seeds": request.maximum_number_of_seeds,
                                  "backend": request.backend_parameters.get("software"),
                                  "connectivity_policy": request.connectivity_policy,
                                  "selection_policy": {k: payload[k] for k in payload if k.startswith("selection_")},
                                  "stages": state.data["completed_steps"], "sampling": payload["sampling"]})
    if status == "completed":
        write_json(directory / "final" / "summary.json", result.to_dict())
    return result
