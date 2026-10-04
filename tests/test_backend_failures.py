"""Backend failures must stop a job instead of looking like empty chemistry."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


INPUTS = Path(__file__).parent / "data" / "neb"


def _run_reaction(tmp_path, fake_xtb, model):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    executable = bin_dir / "xtb"
    executable.write_text("#!/bin/sh\n" + fake_xtb)
    executable.chmod(0o755)
    env = dict(os.environ, PATH=(
        f"{bin_dir}{os.pathsep}{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}"
    ))
    return subprocess.run(
        [sys.executable, "-m", "pyar.cli", "react", str(INPUTS / "hcn.xyz"),
         str(INPUTS / "hnc.xyz"), "-N", "1", "--software", "xtb",
         "--xtb-model", model, "--bias-min", "1000", "--bias-max", "1000",
         "--nprocs", "1"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )


def test_unsupported_xtb_model_exits_before_creating_reaction_state(tmp_path):
    if not Path(sys.executable).with_name("geometric-optimize").is_file():
        pytest.skip("geomeTRIC executable is unavailable")
    result = _run_reaction(tmp_path, "echo 'Usage: xtb --gfn 2'\n", "gxtb")
    assert result.returncode != 0
    assert "does not support --gxtb" in result.stderr
    assert not (tmp_path / "reaction" / "state.json").exists()


def test_backend_input_error_exits_and_persists_failed_state(tmp_path):
    if not Path(sys.executable).with_name("geometric-optimize").is_file() or shutil.which("obabel") is None:
        pytest.skip("reaction workflow dependencies are unavailable")
    result = _run_reaction(tmp_path, "echo 'invalid xTB input' >&2\nexit 2\n", "gfn2")
    assert result.returncode != 0
    state = json.loads((tmp_path / "reaction" / "state.json").read_text())
    assert state["status"] == "failed_backend"
    assert state["products"] == []
    assert state["completed_cycles"] == []
    assert "geometric.out" in state["failure"]
    assert "completed_no_candidates" not in state["status"]


def test_missing_optimizer_input_exits_nonzero(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "pyar.scripts.optimiser", str(tmp_path / "missing.xyz"),
         "-c", "0", "-m", "1", "--software", "xtb"],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "missing.xyz" in result.stderr
