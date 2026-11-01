"""Real process-isolation contracts, independent of the inline mock fixture."""

import ast
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import random
import subprocess
import sys

import numpy as np
import pytest

from pyar.backends.subprocess_utils import run_command, run_output
from pyar.directory_execution import run_in_directory
from pyar.recovery.observations import capture_processes


def forbid_chdir(*args, **kwargs):
    pytest.fail("production attempted to change process cwd")


def test_production_contains_no_working_directory_mutation():
    root = Path(__file__).parents[1] / "pyar"
    offending = []
    for path in root.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
                if node.func.attr in ("chdir", "fchdir"):
                    offending.append(f"{path.relative_to(root)}:{node.lineno}")
    assert offending == []


def test_commands_own_stream_paths_and_leave_caller_untouched(tmp_path, monkeypatch):
    import os

    monkeypatch.setattr(os, "chdir", forbid_chdir)
    caller = Path.cwd()
    code = "from pathlib import Path; import sys; print(Path.cwd()); print('err', file=sys.stderr)"
    assert (
        run_command([sys.executable, "-c", code], "output", "error", cwd=tmp_path) == 0
    )
    assert (tmp_path / "output").read_text().strip() == str(tmp_path)
    assert (tmp_path / "error").read_text().strip() == "err"
    assert run_output([sys.executable, "-c", "print(42)"], cwd=tmp_path) == b"42\n"
    assert Path.cwd() == caller


def test_workers_are_isolated_in_parallel_without_parent_chdir(tmp_path, monkeypatch):
    import os

    monkeypatch.setattr(os, "chdir", forbid_chdir)
    caller = Path.cwd()
    directories = [tmp_path / name for name in ("first", "second")]

    def execute(directory):
        return run_in_directory(
            run_output,
            [
                sys.executable,
                "-c",
                "from pathlib import Path; Path('owned').write_text(str(Path.cwd())); print(Path.cwd())",
            ],
            directory=directory,
        )

    with ThreadPoolExecutor(2) as pool:
        results = list(pool.map(execute, directories))
    assert [result.decode().strip() for result in results] == list(
        map(str, directories)
    )
    assert [(directory / "owned").read_text() for directory in directories] == list(
        map(str, directories)
    )
    assert Path.cwd() == caller


def test_worker_relays_failure_and_retains_artifacts(tmp_path):
    caller = Path.cwd()
    with pytest.raises(subprocess.CalledProcessError) as caught:
        run_in_directory(
            run_output,
            [
                sys.executable,
                "-c",
                "from pathlib import Path; Path('failure-evidence').write_text('kept'); raise SystemExit(7)",
            ],
            directory=tmp_path,
        )
    assert caught.value.returncode == 7
    assert (tmp_path / "failure-evidence").read_text() == "kept"
    assert caught.value.__notes__ and "CalledProcessError" in caught.value.__notes__[0]
    assert Path.cwd() == caller


def test_worker_returns_actual_native_process_observation(tmp_path):
    with capture_processes() as processes:
        assert (
            run_in_directory(
                run_command, [sys.executable, "-c", "pass"], directory=tmp_path
            )
            == 0
        )
    assert len(processes) == 1
    assert processes[0].returncode == 0
    assert processes[0].program == Path(sys.executable).name


def test_worker_preserves_trial_geometry_and_random_progression(tmp_path):
    from pyar.core.molecule import Molecule
    from pyar.sampling.trial_generator import create_trial_geometries

    first = Molecule(["He"], [[0, 0, 0]], name="a")
    second = Molecule(["He"], [[0, 0, 0]], name="b")
    python_state, numpy_state = random.getstate(), np.random.get_state()
    try:
        random.seed(123)
        np.random.seed(123)
        expected = create_trial_geometries(
            "geom", first, second, 2, None, directory=tmp_path / "direct"
        )
        expected_python, expected_numpy = random.random(), np.random.random()
        random.seed(123)
        np.random.seed(123)
        actual = run_in_directory(
            create_trial_geometries,
            "geom",
            first,
            second,
            2,
            None,
            directory=tmp_path / "worker",
        )
        assert [m.name for m in actual] == [m.name for m in expected]
        for left, right in zip(actual, expected):
            np.testing.assert_array_equal(left.coordinates, right.coordinates)
        assert random.random() == expected_python
        assert np.random.random() == expected_numpy
        assert (tmp_path / "worker/trial_vectors.dat").read_bytes() == (
            tmp_path / "direct/trial_vectors.dat"
        ).read_bytes()
    finally:
        random.setstate(python_state)
        np.random.set_state(numpy_state)


def test_single_stage_growth_owns_snapshots_outside_caller(tmp_path, monkeypatch):
    from pyar.core.molecule import Molecule
    from pyar.growth import service

    caller = Path.cwd()
    directory = tmp_path / "owned"
    molecule = Molecule(["He", "He"], [[0, 0, 0], [0, 0, 3]], name="candidate")
    molecule.energy = -1.0
    monkeypatch.setattr(
        service.trial_generation,
        "create_trial_geometries",
        lambda *args, **kwargs: [molecule],
    )
    monkeypatch.setattr(
        service, "_optimize_growth_batch", lambda *args, **kwargs: [True]
    )
    monkeypatch.setattr(
        service.clustering, "choose_geometries", lambda *args, **kwargs: [molecule]
    )
    before = set(caller.iterdir())
    selected = service.add_one(
        "aggregate",
        [molecule],
        molecule,
        1,
        {"software": "aimnet_2", "_two_layer_optimization": False},
        1,
        None,
        connectivity_policy="off",
        directory=directory,
    )
    assert selected == [molecule]
    assert (directory / "selected/result_candidate.xyz").is_file()
    assert (directory / "selected/README.txt").is_file()
    assert set(caller.iterdir()) == before
    assert Path.cwd() == caller


def test_abrupt_worker_death_is_not_success_or_restart_authority(tmp_path):
    from pyar.backend_errors import BackendExecutionError
    caller = Path.cwd()
    code = (
        "from pathlib import Path; import os, signal; "
        "Path('partial').write_text('uncommitted'); "
        "os.kill(os.getppid(), signal.SIGKILL)"
    )
    with pytest.raises(BackendExecutionError, match='Directory operation failed'):
        run_in_directory(run_output, [sys.executable, '-c', code], directory=tmp_path)
    assert (tmp_path / 'partial').read_text() == 'uncommitted'
    assert Path.cwd() == caller
    assert not list(tmp_path.glob('*commit*'))
