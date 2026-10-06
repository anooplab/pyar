"""Orchestration for solute-centred finite microsolvation shells."""

from __future__ import annotations

import math
import os
from pathlib import Path

from pyar.backend_errors import BackendExecutionError
from pyar.microsolvation.placement import assess_shell, generate_microsolvation_candidates
from pyar.microsolvation.request import MicrosolvationRequest
from pyar.state.microsolvation import MicrosolvationRunState
from pyar.state.grow import write_json
from pyar.workflow_results import MicrosolvationResult


def _quality_key(molecule):
    coverage = molecule.microsolvation_coverage["coverage_fraction"]
    energy = molecule.energy
    energy_key = float(energy) if energy is not None and math.isfinite(float(energy)) else math.inf
    return (-coverage, energy_key, molecule.name)


def _deduplicate_and_select(candidates, maximum_seeds):
    from pyar.selection.deduplication import deduplicate_structures
    candidates = sorted(candidates, key=_quality_key)
    result = deduplicate_structures(candidates, ordering="input")
    unique = result["kept"]
    preselected = unique[: max(maximum_seeds, maximum_seeds * 4)]
    selected = []
    remaining = list(preselected)
    while remaining and len(selected) < maximum_seeds:
        if not selected:
            selected.append(remaining.pop(0))
            continue
        from pyar.structure_comparison import GraphRMSDComparator
        comparator = GraphRMSDComparator(atom_mode="heavy")
        scored = []
        for candidate in remaining:
            distances = []
            for chosen in selected:
                comparison = comparator.compare(candidate, chosen)
                if comparison.distance is not None:
                    distances.append(comparison.distance)
                elif comparison.equivalent is False:
                    distances.append(float("inf"))
                else:
                    distances.append(0.0)
            scored.append((min(distances), _quality_key(candidate), candidate))
        # Quality ranks the eligible pool first; max-min structural separation
        # selects among the best-coverage/energy candidates.
        best_quality = min(item[1] for item in scored)
        quality_frontier = [item for item in scored if item[1][0] <= best_quality[0] + 0.03]
        chosen = max(quality_frontier, key=lambda item: (item[0], tuple(-ord(c) for c in item[2].name)))[2]
        selected.append(chosen)
        remaining.remove(chosen)
    return selected, {
        "selection_mode": "solute-surface-coverage-then-energy-and-structural-diversity",
        "energy_ranked": any(m.energy is not None for m in selected),
        "input_count": len(candidates),
        "unique_count": len(unique),
        "selected_count": len(selected),
        "maximum_number_of_seeds": maximum_seeds,
        "adaptive_rmsd_threshold_angstrom": result["rmsd_threshold"],
        "comparison_diagnostics": result["diagnostics"],
        "coverage_frontier_fraction_tolerance": 0.03,
    }


def _working_directory(path):
    class WorkingDirectory:
        def __enter__(self):
            self.old = os.getcwd()
            os.chdir(path)

        def __exit__(self, *_):
            os.chdir(self.old)
    return WorkingDirectory()


def _result(state):
    directory = state.directory
    paths = tuple(str((directory / reference["path"]).resolve())
                  for reference in state.data.get("final_seeds", []))
    request = state.data["request"]
    return MicrosolvationResult(
        workflow="microsolvation", status=state.data["status"],
        run_directory=str(directory), state_path=str(state.path), selected_paths=paths,
        metadata={
            "requested_solvents": request["count"],
            "completed_solvents": len(state.data.get("completed_steps", [])),
            "solute_atom_indices": state.data["solute_atom_indices"],
            "solvent_fragments": state.data.get("solvent_fragments", []),
            "coverage_history": state.data.get("coverage_history", []),
            "selection_policy": "coverage, then backend energy when present, then structural diversity",
            "sampling": request["sampling"],
            "backend": request["backend_parameters"].get("software"),
            "confinement": request["confinement"],
        },
    )


def microsolvate(request: MicrosolvationRequest, *, output="microsolvation"):
    """Generate a bounded ensemble by targeting original-solute surface sites.

    The solute atom prefix and its surface define placement targets at every
    step. The full cluster is used for steric rejection and backend energy.
    Optimized candidates are screened by a documented shell-validity
    post-filter; no optimization restraint is claimed or applied.
    """
    directory = Path(output).resolve()
    payload = request.to_state_dict()
    state = MicrosolvationRunState.load(directory, payload)
    if state is None:
        initial = request.solute.copy()
        initial.solute_atom_indices = tuple(range(len(initial)))
        initial.solvent_fragments = ()
        initial.microsolvation_coverage = {
            "surface_points": 0, "covered_points": 0, "open_points": 0,
            "coverage_fraction": 0.0, "whole_solute_coverage_fraction": 0.0,
        }
        state = MicrosolvationRunState.create(directory, payload, initial)
    if state.data["status"] == "completed":
        return _result(state)
    if state.data["status"] not in {"running", "stopped"}:
        raise RuntimeError(f"Cannot resume microsolvation state in status {state.data['status']!r}")
    from pyar.microsolvation.surface import build_solute_surface
    initial_surface = build_solute_surface(
        request.solute, points_per_atom=request.surface_points_per_atom,
        probe_radius=request.probe_radius,
    )
    write_json(directory / "surface" / "initial_surface.json", {
        "model": "vdw-probe-fibonacci-v1",
        "points_per_atom": request.surface_points_per_atom,
        "probe_radius_angstrom": request.probe_radius,
        "points_angstrom": initial_surface.points.tolist(),
        "parent_solute_atoms": initial_surface.parent_atoms.tolist(),
        "normals": initial_surface.normals.tolist(),
        "weights_angstrom2": initial_surface.weights.tolist(),
    })
    seeds = state.restore_pool()
    state.data["status"] = "running"
    for step in range(state.data["next_step"], request.count + 1):
        stage = directory / f"step_{step:03d}"
        candidate_directory = stage / "candidates"
        candidate_directory.mkdir(parents=True, exist_ok=True)
        candidates = []
        generated_count = 0
        shell_rejected = 0
        target_rejected = 0
        for seed_index, seed in enumerate(seeds):
            try:
                generated = generate_microsolvation_candidates(
                    seed, request, step=step * 10000 + seed_index,
                    solute_atom_count=len(request.solute),
                )
            except ValueError:
                target_rejected += 1
                continue
            for candidate in generated:
                generated_count += 1
                coverage = assess_shell(candidate, request)
                candidate.microsolvation_coverage = coverage
                (candidate_directory / f"{candidate.name}.xyz").parent.mkdir(parents=True, exist_ok=True)
                candidate.mol_to_xyz(str(candidate_directory / f"{candidate.name}.xyz"))
                if request.backend_parameters.get("software"):
                    qc = dict(request.backend_parameters)
                    qc.update(charge=candidate.charge, multiplicity=candidate.multiplicity,
                              scftype=candidate.scftype, opt_target="minimum", gamma=None)
                    try:
                        with _working_directory(candidate_directory):
                            from pyar.optimiser import is_usable, optimise
                            status = optimise(candidate, qc)
                        if not is_usable(status) or candidate.coordinates is None:
                            continue
                    except BackendExecutionError:
                        raise
                    except Exception as exc:
                        raise BackendExecutionError(f"Microsolvation optimization failed for {candidate.name}: {exc}") from exc
                    if candidate.coordinates is None or not candidate.solvent_fragments:
                        continue
                    # Reject remote solvent-droplet collapse after relaxation.
                    coverage = assess_shell(candidate, request)
                    candidate.microsolvation_coverage = coverage
                    if not coverage["shell_valid"]:
                        shell_rejected += 1
                        continue
                    candidate.mol_to_xyz(str(candidate_directory / f"{candidate.name}_optimized.xyz"))
                else:
                    if not coverage["shell_valid"]:
                        shell_rejected += 1
                        continue
                candidates.append(candidate)
        selected, diagnostics = _deduplicate_and_select(candidates, request.maximum_number_of_seeds)
        diagnostics.update({
            "selection_basis": "geometry-only" if not request.backend_parameters.get("software") else "surface-quality-then-energy",
            "surface_saturated_at_placement": any(
                getattr(candidate, "surface_saturated_at_placement", False) for candidate in candidates
            ),
            "generated_count": generated_count,
            "rejected_outside_shell": shell_rejected,
            "seeds_without_accessible_target_surface": target_rejected,
        })
        write_json(stage / "selection_diagnostics.json", diagnostics)
        if not selected:
            state.data["status"] = "no_candidates"
            state.save()
            break
        state.complete_step(step, selected, diagnostics)
        state.data["solute_atom_indices"] = list(selected[0].solute_atom_indices)
        state.data["solvent_fragments"] = [list(part) for part in selected[0].solvent_fragments]
        state.save()
        seeds = selected
    else:
        state.data["status"] = "completed"
    if state.data["status"] == "completed":
        references = state.snapshot_pool(seeds, directory / "final")
        state.data["final_seeds"] = references
        write_json(directory / "final" / "summary.json", {
            "workflow": "microsolvation",
            "status": "completed",
            "requested_solvents": request.count,
            "selected_count": len(references),
            "coverage": [reference["coverage"] for reference in references],
            "selection_policy": "surface validity and coverage precede energy and structural diversity",
        })
    state.save()
    return _result(state)
