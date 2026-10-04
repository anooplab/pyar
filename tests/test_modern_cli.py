"""Compatibility checks for the task-oriented clustering pilot."""

import importlib.metadata
import subprocess
import sys
import sysconfig
from pathlib import Path
from unittest import mock

import pytest

from pyar import __version__, modern_cli
from pyar.scripts import clustering


@pytest.mark.parametrize("argv, expected", [
    (["--help"], "clustering"),
    (["--version"], __version__),
    (["clustering", "--help"], "--maximum-number-of-seeds"),
])
def test_help_and_version(argv, expected, capsys):
    with pytest.raises(SystemExit) as exc:
        modern_cli.main(argv)
    assert exc.value.code == 0
    output = capsys.readouterr().out
    assert expected in output
    if argv == ["clustering", "--help"]:
        assert "pyar clustering" in output
        for action in clustering.build_parser()._actions:
            for option in action.option_strings:
                assert option in output


@pytest.fixture
def xyz_files(tmp_path):
    paths = []
    for name, distance in [("a", 0.7), ("b", 0.9), ("c", 1.1)]:
        path = tmp_path / f"{name}.xyz"
        path.write_text(f"2\n{name}: energy={distance}\nH 0 0 0\nH {distance} 0 0\n")
        paths.append(str(path))
    return paths


@pytest.mark.parametrize("options, expected", [
    ([], {"algorithm": "auto", "feature": "auto", "maximum_number_of_seeds": 12}),
    (["--algorithm", "agglomerative", "--feature", "distance-histogram",
      "--maximum-number-of-seeds", "5"],
     {"algorithm": "agglomerative", "feature": "distance-histogram",
      "maximum_number_of_seeds": 5}),
])
def test_shared_execution_and_defaults(xyz_files, options, expected, capsys):
    calls = []
    outputs = []
    for entry, prefix in [(clustering.main, []), (modern_cli.main, ["clustering"])]:
        with mock.patch.object(clustering.clustering, "choose_geometries",
                               side_effect=lambda molecules, **kwargs: molecules[:1]) as chooser:
            entry([*prefix, *xyz_files, *options])
        calls.append(chooser.call_args)
        outputs.append(capsys.readouterr().out)
    assert outputs[0] == outputs[1]
    assert [mol.relative_path for mol in calls[1].args[0]] == xyz_files
    assert calls[0].kwargs == calls[1].kwargs
    for key, value in expected.items():
        assert calls[1].kwargs[key] == value
    assert calls[1].kwargs["distance_metric"] == "euclidean"
    assert calls[1].kwargs["system_type"] == "auto"
    assert calls[1].kwargs["algorithm_options"] == {
        "min_samples": 2, "min_cluster_size": 2, "eps": None, "xi": 0.05,
    }
    assert calls[1].kwargs["distance_options"] == {}


@pytest.mark.parametrize("options, message", [
    (["--algorithm", "invalid"], "invalid choice"),
    (["--distance", "invalid"], "invalid choice"),
    (["--soap-cutoff", "-1"], "soap_cutoff"),
    (["--labels-output", "labels.csv"], "labels mode"),
    (["--mode", "filter", "--report-output", "report.json"], "cluster and labels"),
    (["--structure-report", "structure.json"], "analyze mode"),
    (["--mode", "analyze", "--coordinate-model", "distance-cutoff"], "requires --bond-cutoff"),
    (["--mode", "analyze", "--bond-cutoff", "2"], "requires --coordinate-model"),
])
def test_invalid_options(xyz_files, options, message, capsys):
    with pytest.raises(SystemExit) as exc:
        modern_cli.main(["clustering", *xyz_files, *options])
    assert exc.value.code != 0
    assert message in capsys.readouterr().err


def test_too_few_inputs(xyz_files, capsys):
    with pytest.raises(SystemExit) as exc:
        modern_cli.main(["clustering", xyz_files[0]])
    assert exc.value.code == 2
    assert "at least two" in capsys.readouterr().err


@pytest.mark.parametrize("contents", [None, "not an XYZ file\n"])
def test_missing_or_malformed_xyz(tmp_path, contents):
    path = tmp_path / "bad.xyz"
    if contents is not None:
        path.write_text(contents)
    with pytest.raises(SystemExit) as exc:
        modern_cli.main(["clustering", str(path), str(path)])
    assert exc.value.code != 0


@pytest.mark.parametrize("argv", [[], ["react"], ["clustering", "--unknown"]])
def test_invalid_command(argv):
    with pytest.raises(SystemExit) as exc:
        modern_cli.main(argv)
    assert exc.value.code == 2


def test_top_level_does_not_load_clustering():
    result = subprocess.run(
        [sys.executable, "-c",
         "import sys; from pyar.modern_cli import build_parser; build_parser(); "
         "assert 'pyar.scripts.clustering' not in sys.modules; "
         "assert 'pyar.selection.clusterers' not in sys.modules"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr


def test_installed_entry_points():
    # Source-tree egg-info can mask the editable installation's fresh metadata.
    distribution = next(dist for dist in importlib.metadata.distributions(
        path=[sysconfig.get_path("purelib")]) if dist.metadata["Name"] == "pyar-chem")
    entries = {ep.name: ep.value for ep in distribution.entry_points
               if ep.group == "console_scripts"}
    assert entries["pyar"] == "pyar.modern_cli:main"
    assert entries["pyar-cli"] == "pyar.cli:main"
    assert entries["pyar-clustering"] == "pyar.scripts.clustering:main"
    assert not Path("pyar.py").exists()
