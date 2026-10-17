"""Orchestration for solute-centred finite microsolvation shells."""

from __future__ import annotations

import math
import json
import os
from pathlib import Path

from pyar.microsolvation.placement import assess_shell, generate_microsolvation_candidates
from pyar.microsolvation.request import MicrosolvationRequest
from pyar.state.microsolvation import MicrosolvationRunState
from pyar.state.grow import write_json
from pyar.workflow_results import MicrosolvationResult


def _quality_key(molecule, request):
    coverage = molecule.microsolvation_coverage["coverage_fraction"]
    energy = molecule.energy
    energy_key = float(energy) if energy is not None and math.isfinite(float(energy)) else math.inf
    if request.objective == "cluster" or request.wall != "none":
        return (energy_key, -coverage, molecule.name)
    if request.placement_score == "energy":
        return (-coverage, energy_key, molecule.name)
    return (-coverage, molecule.name)


def _deduplicate_and_select(candidates, maximum_seeds, request):
    from pyar.selection.deduplication import deduplicate_structures
    candidates = sorted(candidates, key=lambda molecule: _quality_key(molecule, request))
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
            scored.append((min(distances), _quality_key(candidate, request), candidate))
        # Keep the best objective stratum and maximize structural separation
        # within it; coverage must not erase a low-energy clustered candidate.
        best_quality = min(item[1] for item in scored)
        if request.objective == "cluster" or request.wall != "none":
            quality_frontier = [item for item in scored if item[1][0] <= best_quality[0] + 0.005]
        else:
            quality_frontier = [item for item in scored if item[1][0] <= best_quality[0] + 0.03]
        chosen = max(quality_frontier, key=lambda item: (item[0], tuple(-ord(c) for c in item[2].name)))[2]
        selected.append(chosen)
        remaining.remove(chosen)
    return selected, {
        "selection_mode": ("fixed-composition-energy-then-diversity" if request.objective == "cluster"
                           else "unbiased-post-relaxation-energy-then-diversity" if request.wall != "none"
                           else "solute-surface-coverage-then-energy-and-structural-diversity"),
        "objective": request.objective,
        "placement_score": request.placement_score,
        "energy_unit": "hartree" if request.placement_score == "energy" or request.wall != "none" else None,
        "energy_ranked": any(m.energy is not None for m in selected),
        "input_count": len(candidates),
        "unique_count": len(unique),
        "selected_count": len(selected),
        "maximum_number_of_seeds": maximum_seeds,
        "adaptive_rmsd_threshold_angstrom": result["rmsd_threshold"],
        "comparison_diagnostics": result["diagnostics"],
        "coverage_frontier_fraction_tolerance": 0.03 if request.objective == "shell" and request.wall == "none" else None,
        "cluster_energy_frontier_hartree": 0.005 if request.objective == "cluster" else None,
        "unbiased_energy_frontier_hartree": 0.005 if request.wall != "none" else None,
    }


def _working_directory(path):
    class WorkingDirectory:
        def __enter__(self):
            Path(path).mkdir(parents=True, exist_ok=True)
            self.old = os.getcwd()
            os.chdir(path)

        def __exit__(self, *_):
            os.chdir(self.old)
    return WorkingDirectory()


def _energy_prescreen(candidates, request):
    """Score complete, fixed-composition candidates before expensive relaxations."""
    if request.placement_score != "energy":
        return list(candidates), []
    if not candidates:
        return [], []
    from pyar.data.units import angstrom2bohr
    from pyar.energy_gradient_providers import get_energy_gradient_provider
    from ase import Atoms

    first = candidates[0]
    qc = dict(request.backend_parameters)
    qc.update(charge=first.charge, multiplicity=first.multiplicity,
              scftype=first.scftype, gamma=None)
    provider = get_energy_gradient_provider(qc["software"], qc)
    records = []
    for candidate in candidates:
        candidate_qc = dict(qc, charge=candidate.charge, multiplicity=candidate.multiplicity,
                            scftype=candidate.scftype)
        try:
            candidate_provider = provider if candidate_qc == qc else get_energy_gradient_provider(
                candidate_qc["software"], candidate_qc)
            backend_molecule = Atoms(symbols=candidate.atoms_list, positions=candidate.coordinates)
            result = candidate_provider.evaluate(backend_molecule, angstrom2bohr(candidate.coordinates))
            energy = float(result.energy_hartree)
            if not math.isfinite(energy):
                raise ValueError("backend returned a non-finite prescreen energy")
            candidate.microsolvation_prescreen_energy = energy
            records.append({"name": candidate.name, "energy_hartree": energy, "status": "ok"})
        except Exception as exc:
            records.append({"name": candidate.name, "energy_hartree": None, "status": "failed",
                            "error": str(exc)[:500]})
    policy = request.to_state_dict()["selection_policy"]
    budget = request.maximum_number_of_seeds * policy["prescreen_pool_multiplier"]
    scored = [item for item in candidates if hasattr(item, "microsolvation_prescreen_energy")]
    unscored = [item for item in candidates if not hasattr(item, "microsolvation_prescreen_energy")]
    if request.objective == "cluster":
        ranked = sorted(scored, key=lambda item: (item.microsolvation_prescreen_energy, item.name))
        chosen = ranked[:budget] + unscored
    else:
        by_energy = sorted(scored, key=lambda item: (item.microsolvation_prescreen_energy, item.name))[:budget]
        by_coverage = sorted(candidates, key=lambda item: (
            -item.microsolvation_coverage["coverage_fraction"], item.name))[:budget]
        chosen = list({item.name: item for item in by_energy + by_coverage + unscored}.values())
    keep = {item.name for item in chosen}
    for record in records:
        record["selected_for_optimization"] = record["name"] in keep
    return chosen, records


def _optimize_candidates(candidates, request, candidate_directory):
    """Optimize the selected prescreen pool and apply objective-specific validity."""
    if not request.backend_parameters.get("software"):
        return list(candidates), 0, 0, []
    selected = []
    shell_rejected = 0
    failures = []
    for candidate in candidates:
        qc = dict(request.backend_parameters)
        qc.update(charge=candidate.charge, multiplicity=candidate.multiplicity,
                  scftype=candidate.scftype, opt_target="minimum", gamma=None)
        try:
            from pyar.optimiser import is_usable, optimise
            if request.wall != "none":
                qc["microsolvation_wall"] = {
                    "solute_atom_count": len(request.solute),
                    "solvent_fragments": [list(fragment) for fragment in candidate.solvent_fragments],
                    "radius_angstrom": float(request.wall_radius),
                    "force_constant_hartree_per_angstrom2": float(request.wall_force_constant),
                }
                with _working_directory(candidate_directory / "wall_search"):
                    status = optimise(candidate, qc)
                if not is_usable(status) or candidate.coordinates is None:
                    failures.append({"name": candidate.name, "stage": "wall_search",
                                     "status": str(status or "no coordinates")})
                    continue
                try:
                    wall_state = json.loads((candidate_directory / "wall_search" / f"job_{candidate.name}" /
                                             "pyar_geometric_state.json").read_text())
                    candidate.microsolvation_wall_objective_energy = float(wall_state["total_energy_hartree"])
                    candidate.microsolvation_wall_energy = float(wall_state["wall_energy_hartree"])
                    candidate.microsolvation_wall_energy_max = float(wall_state["wall_energy_max_hartree"])
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    raise RuntimeError(f"missing valid wall-search energy provenance for {candidate.name}: {exc}") from exc
                qc.pop("microsolvation_wall")
                with _working_directory(candidate_directory / "unbiased_relaxation"):
                    status = optimise(candidate, qc)
            else:
                with _working_directory(candidate_directory):
                    status = optimise(candidate, qc)
            if not is_usable(status) or candidate.coordinates is None:
                failures.append({"name": candidate.name, "status": str(status or "no coordinates")})
                continue
        except Exception as exc:
            failures.append({"name": candidate.name, "status": "error", "error": str(exc)[:500]})
            continue
        if candidate.coordinates is None or not candidate.solvent_fragments:
            continue
        coverage = assess_shell(candidate, request)
        candidate.microsolvation_coverage = coverage
        if request.wall != "none":
            candidate.microsolvation_wall_survived_unbiased_relaxation = bool(coverage["shell_valid"])
            candidate.microsolvation_physical_energy = float(candidate.energy)
        if request.objective == "shell" and not coverage["shell_valid"]:
            shell_rejected += 1
            continue
        candidate.mol_to_xyz(str(candidate_directory / f"{candidate.name}_optimized.xyz"))
        selected.append(candidate)
    return selected, shell_rejected, len(candidates), failures


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
            "selection_policy": request.get("selection_policy", {"name": "legacy-shell-v1"}),
            "sampling": request["sampling"],
            "backend": request["backend_parameters"].get("software"),
            "objective": request.get("selection_policy", {}).get("objective", "shell"),
            "confinement": request["confinement"],
        },
    )


def _microsolvate_v1(request: MicrosolvationRequest, *, output="microsolvation"):
    """Generate a bounded ensemble under the request's explicit objective."""
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
                candidates.append(candidate)
        prescreen_records = []
        if request.placement_score == "energy":
            candidates, prescreen_records = _energy_prescreen(candidates, request)
            write_json(stage / "energy_prescreen.json", {
                "method": request.backend_parameters,
                "energy_unit": "hartree",
                "fixed_composition_per_step": True,
                "candidates": prescreen_records,
            })
        if request.backend_parameters.get("software"):
            candidates, shell_rejected, optimized_count, optimization_failures = _optimize_candidates(
                candidates, request, candidate_directory)
        elif request.objective == "shell":
            retained = [candidate for candidate in candidates
                        if candidate.microsolvation_coverage["shell_valid"]]
            shell_rejected = len(candidates) - len(retained)
            candidates = retained
            optimized_count = 0
            optimization_failures = []
        else:
            optimized_count = 0
            optimization_failures = []
        selected, diagnostics = _deduplicate_and_select(candidates, request.maximum_number_of_seeds, request)
        diagnostics.update({
            "selection_basis": ("fixed-composition-energy" if request.objective == "cluster" else
                                "unbiased-post-relaxation-energy" if request.wall != "none" else
                                "geometry-only" if not request.backend_parameters.get("software") else
                                "surface-quality-then-energy"),
            "surface_saturated_at_placement": any(
                getattr(candidate, "surface_saturated_at_placement", False) for candidate in candidates
            ),
            "generated_count": generated_count,
            "prescreened_count": len(prescreen_records),
            "optimized_count": optimized_count,
            "optimization_failures": optimization_failures,
            "rejected_outside_shell": shell_rejected,
            "seeds_without_accessible_target_surface": target_rejected,
        })
        if request.wall != "none":
            diagnostics["wall_candidate_outcomes"] = [
                {
                    "name": candidate.name,
                    "wall_objective_energy_hartree": getattr(candidate, "microsolvation_wall_objective_energy", None),
                    "wall_energy_hartree": getattr(candidate, "microsolvation_wall_energy", None),
                    "wall_energy_max_hartree": getattr(candidate, "microsolvation_wall_energy_max", None),
                    "unbiased_physical_energy_hartree": (
                        None if candidate.energy is None else float(candidate.energy)),
                    "shell_survived_unbiased_relaxation": bool(
                        candidate.microsolvation_coverage.get("shell_valid", False)),
                }
                for candidate in candidates
            ]
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
            "selection_policy": request.to_state_dict().get("selection_policy", {"name": "legacy-shell-v1"}),
        })
    state.save()
    return _result(state)


def microsolvate(request: MicrosolvationRequest, *, output="microsolvation", recover_stale=False):
    """Create fresh v2 runs and keep established v1 continuation unchanged."""
    from pyar.persistence import detect_version, RunLock
    directory = Path(output).resolve()
    payload = request.to_state_dict()
    def legacy(root):
        return MicrosolvationRunState.load(root, payload) is not None
    version = detect_version(directory, 'microsolvation', validate_v1=legacy)
    if version == 1:
        return _microsolvate_v1(request, output=output)
    directory.parent.mkdir(parents=True, exist_ok=True)
    with RunLock(directory, recover_stale=recover_stale):
        version = detect_version(directory, 'microsolvation', validate_v1=legacy)
        if version == 1:
            return _microsolvate_v1(request, output=output)
        return _microsolvate_v2(request, output=output)


def _microsolvate_v2(request: MicrosolvationRequest, *, output="microsolvation"):
    """Generate a bounded ensemble under the request's explicit objective."""
    from pyar.microsolvation.layout_v2 import MicrosolvationRunStateV2
    from pyar.persistence import RunLayout
    directory = Path(output).resolve()
    layout = RunLayout(directory, "microsolvation")
    payload = request.to_state_dict()
    state = MicrosolvationRunStateV2.load(directory, payload) if layout.marker.exists() else None
    if state is None:
        initial = request.solute.copy()
        initial.solute_atom_indices = tuple(range(len(initial)))
        initial.solvent_fragments = ()
        initial.microsolvation_coverage = {
            "surface_points": 0, "covered_points": 0, "open_points": 0,
            "coverage_fraction": 0.0, "whole_solute_coverage_fraction": 0.0,
        }
        state = MicrosolvationRunStateV2.create(directory, payload, initial, request.solvent)
    if state.data["status"] == "completed":
        state.report_results()
        return _result(state)
    if state.data["status"] not in {"running", "stopped"}:
        raise RuntimeError(f"Cannot resume microsolvation state in status {state.data['status']!r}")
    from pyar.microsolvation.surface import build_solute_surface
    initial_surface = build_solute_surface(
        request.solute, points_per_atom=request.surface_points_per_atom,
        probe_radius=request.probe_radius,
    )
    state.publish_surface({
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
        if state.selection_path(step).exists():
            state.accept_step(step)
            seeds = state.restore_pool()
            continue
        stage = layout.path(f"steps/step_{step:03d}")
        candidate_directory = layout.path(f"steps/step_{step:03d}/attempts/attempt_000/execution/candidates")
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
                candidates.append(candidate)
        prescreen_records = []
        if request.placement_score == "energy":
            candidates, prescreen_records = _energy_prescreen(candidates, request)
            write_json(layout.path(f"steps/step_{step:03d}/energy_prescreen.json"), {
                "method": request.backend_parameters,
                "energy_unit": "hartree",
                "fixed_composition_per_step": True,
                "candidates": prescreen_records,
            })
        if request.backend_parameters.get("software"):
            candidates, shell_rejected, optimized_count, optimization_failures = _optimize_candidates(
                candidates, request, candidate_directory)
        elif request.objective == "shell":
            retained = [candidate for candidate in candidates
                        if candidate.microsolvation_coverage["shell_valid"]]
            shell_rejected = len(candidates) - len(retained)
            candidates = retained
            optimized_count = 0
            optimization_failures = []
        else:
            optimized_count = 0
            optimization_failures = []
        selected, diagnostics = _deduplicate_and_select(candidates, request.maximum_number_of_seeds, request)
        diagnostics.update({
            "selection_basis": ("fixed-composition-energy" if request.objective == "cluster" else
                                "unbiased-post-relaxation-energy" if request.wall != "none" else
                                "geometry-only" if not request.backend_parameters.get("software") else
                                "surface-quality-then-energy"),
            "surface_saturated_at_placement": any(
                getattr(candidate, "surface_saturated_at_placement", False) for candidate in candidates
            ),
            "generated_count": generated_count,
            "prescreened_count": len(prescreen_records),
            "optimized_count": optimized_count,
            "optimization_failures": optimization_failures,
            "rejected_outside_shell": shell_rejected,
            "seeds_without_accessible_target_surface": target_rejected,
        })
        if request.wall != "none":
            diagnostics["wall_candidate_outcomes"] = [
                {"name": candidate.name,
                 "wall_objective_energy_hartree": getattr(candidate, "microsolvation_wall_objective_energy", None),
                 "wall_energy_hartree": getattr(candidate, "microsolvation_wall_energy", None),
                 "wall_energy_max_hartree": getattr(candidate, "microsolvation_wall_energy_max", None),
                 "unbiased_physical_energy_hartree": None if candidate.energy is None else float(candidate.energy),
                 "shell_survived_unbiased_relaxation": bool(candidate.microsolvation_coverage.get("shell_valid", False))}
                for candidate in candidates
            ]
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
        state.finish_results()
        state.report_results()
    else:
        state.save()
    return _result(state)
