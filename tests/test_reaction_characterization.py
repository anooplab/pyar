import json
from types import SimpleNamespace
from unittest.mock import patch

from pyar.reaction_request import resolve_reaction_request
from pyar.workflows.reaction_characterization import characterize_reaction, validate_characterization_restart
from pyar.workflows.scan_path import run_path_stages
from pyar.neb import validate_stage_input_paths


def _xyz(path, atoms, coords):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{len(atoms)}\nfixture\n" + "".join(
        f"{atom} {x} {y} {z}\n" for atom, (x, y, z) in zip(atoms, coords)))
    return path


def test_run_path_stages_reuses_cumulative_engine_and_unbiased_qc(tmp_path):
    calls = []

    def run_neb(**kwargs):
        calls.append(kwargs)
        return {"stage": kwargs["stage"], "converged": True}

    with patch("pyar.neb.run_neb", side_effect=run_neb), patch(
            "pyar.neb._stage_gate_passed", return_value=True):
        result = run_path_stages(
            start="r.xyz", end="p.xyz", ts_guess="guess.xyz", output=tmp_path,
            qc_params={"software": "xtb", "gamma": 100, "bias_controller": "adaptive",
                       "bias_potential": "afir", "charge": 0, "multiplicity": 1,
                       "xtb_model": "gfn2", "scftype": "rhf", "scf_cycles": 99},
            through="neb", options={"images": 5},
        )
    assert result["status"] == "complete"
    assert [call["stage"] for call in calls] == ["relax", "neb"]
    assert all("gamma" not in call for call in calls)  # run_neb constructs an unbiased gamma=0 calculator.
    assert all("bias_controller" not in call and call["software"] == "xtb" for call in calls)
    assert calls[0]["backend_options"] == {"scftype": "rhf", "scf_cycles": 99}


def test_run_path_stages_all_uses_existing_cumulative_stage_sequence(tmp_path):
    calls = []
    with patch("pyar.neb.run_neb", side_effect=lambda **kw: calls.append(kw) or {"stage": kw["stage"]}), \
            patch("pyar.neb._stage_gate_passed", return_value=True):
        result = run_path_stages(start="r.xyz", end="p.xyz", ts_guess="guess.xyz", output=tmp_path,
                                 qc_params={"software": "xtb", "charge": 0, "multiplicity": 1},
                                 through="all", options={})
    assert [call["stage"] for call in calls] == [
        "relax", "neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency"]
    assert result["status"] == "complete"


def test_characterization_uses_latest_reactant_trace_and_records_route(tmp_path):
    reaction = tmp_path / "reaction"
    job = reaction / "gamma_0100" / "orientation_x" / "job_gamma_0100_x"
    candidates = job / "candidate_ts"
    candidates.mkdir(parents=True)
    product_path = _xyz(reaction / "products" / "accepted.xyz", ["H", "H"],
                        [(0, 0, 0), (0.74, 0, 0)])
    initial_path = _xyz(job.parent / "trial_gamma_0100_x.xyz", ["H", "H"],
                         [(0, 0, 0), (3.0, 0, 0)])
    _xyz(candidates / "highest_backend_energy.xyz", ["H", "H"],
         [(0, 0, 0), (1.0, 0, 0)])
    (candidates / "metadata.json").write_text("{}\n")
    trace = job / "reaction_trace"
    trace.mkdir()
    (trace / "trace.jsonl").write_text("trace metadata fixture\n")
    steps = trace / "steps"
    _xyz(steps / "step_000020.xyz", ["H", "H"], [(0, 0, 0), (2.0, 0, 0)])
    records = [
        {"step_index": 10, "symbols": ["H", "H"], "coordinates_angstrom": [[0, 0, 0], [2.4, 0, 0]],
         "current_bonds": []},
        {"step_index": 20, "symbols": ["H", "H"], "coordinates_angstrom": [[0, 0, 0], [2.0, 0, 0]],
         "current_bonds": []},
        {"step_index": 40, "symbols": ["H", "H"], "coordinates_angstrom": [[0, 0, 0], [0.74, 0, 0]],
         "current_bonds": [[0, 1]]},
    ]
    product = {"job_name": "gamma_0100_x", "gamma": 100.0,
               "path": "products/accepted.xyz", "inchi": "InChI=1/H2", "smiles": "[HH]",
               "trace_summary": {"candidate_ts_directory": str(candidates)}}
    request_data = {"backend_parameters": {"software": "xtb"}, "gamma_schedule": [100.0]}
    (reaction / "state.json").write_text(json.dumps({
        "version": 2, "workflow": "reaction", "status": "completed_products_found",
        "request": request_data, "products": [product], "completed_jobs": [],
    }))
    request = SimpleNamespace(restart_request=request_data,
                               qc_params={"software": "xtb", "xtb_model": "gfn2", "charge": 0,
                                          "multiplicity": 1, "scftype": "rhf", "nprocs": 1})
    options = {"pathway_ts_source": "highest-backend-energy", "images": 5}
    mock_path = {"status": "complete", "completed_stages": ["relax", "neb"]}
    with patch("pyar.workflows.reaction_characterization.load_trace_records", return_value=records), \
            patch("pyar.workflows.reaction_characterization._persistent_transition_index", return_value=2), \
            patch("pyar.workflows.reaction_characterization.run_path_stages", return_value=mock_path) as run:
        result = characterize_reaction(tmp_path, request, "neb", options)
        route_file = reaction / "pathways" / "product_001" / "route_001" / "pathway.json"
        route = json.loads(route_file.read_text())
        assert route["reactant_geometry_source"] == "latest_valid_pre_transition_trace_frame"
        assert route["reactant_trace_step"] == 20
        assert route["candidate_ts_source"] == "highest_backend_energy"
        assert route["candidate_ts_is_validated_transition_state"] is False
        assert result["status"] == "complete"
        assert run.call_args.kwargs["through"] == "neb"
        assert run.call_args.kwargs["ts_guess"].relative_to(run.call_args.kwargs["output"]).as_posix() == \
            "inputs/ts_guess.xyz"
        validate_stage_input_paths(run.call_args.kwargs["output"], "neb",
                                   ts_guess=run.call_args.kwargs["ts_guess"])
        from pyar.workflows.reaction_characterization import _xyz_data
        import numpy as np
        assert np.allclose(_xyz_data(run.call_args.kwargs["start"])[1],
                           _xyz_data(steps / "step_000020.xyz")[1])
        # A failed run created by the earlier bridge version stored the
        # waypoint at route_dir/ts_guess.xyz. Accept that read-only state and
        # migrate the immutable input before resuming NEB.
        old_waypoint = route_file.parent / "ts_guess.xyz"
        new_waypoint = route_file.parent / "inputs" / "ts_guess.xyz"
        old_waypoint.write_bytes(new_waypoint.read_bytes())
        new_waypoint.unlink()
        route.pop("pathway_input_files")
        route_file.write_text(json.dumps(route, indent=2, sort_keys=True) + "\n")
        # Identical reruns retain the same route files and invoke the engine in
        # the same output directory, where NEB's own stage hashes govern reuse.
        before = {path.relative_to(reaction): path.read_bytes() for path in (reaction / "pathways").rglob("*")
                  if path.is_file()}
        validate_characterization_restart(tmp_path, request, "neb",
                                          dict(options, through="neb"))
        after = {path.relative_to(reaction): path.read_bytes() for path in (reaction / "pathways").rglob("*")
                 if path.is_file()}
        assert before == after
        characterize_reaction(tmp_path, request, "neb", options)
        assert run.call_count == 2
        resumed_route = json.loads(route_file.read_text())
        assert resumed_route["pathway_input_files"]["ts_guess"] == "inputs/ts_guess.xyz"
        assert run.call_args.kwargs["ts_guess"] == new_waypoint
        assert run.call_args.kwargs["output"] == run.call_args_list[0].kwargs["output"]
