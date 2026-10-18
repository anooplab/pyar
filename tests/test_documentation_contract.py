from types import SimpleNamespace

import pytest

from tools.check_documentation import check


def make_document_root(
    tmp_path,
    *,
    command="pyar react input.xyz --valid",
    python="3.10",
    extra="demo",
    scripts=("pyar", "pyar-cli"),
):
    (tmp_path / "docs").mkdir()
    (tmp_path / ".github/workflows").mkdir(parents=True)
    (tmp_path / "pyproject.toml").write_text(
        '[project]\nname="example"\nversion="1.0"\nrequires-python=">=' + python + '"\n'
        'classifiers=["Programming Language :: Python :: 3.10", "Programming Language :: Python :: 3.11"]\n'
        "[project.optional-dependencies]\n" + extra + "=[]\n"
        "[project.scripts]\n" + "".join(f'{name}="example:main"\n' for name in scripts)
    )
    (tmp_path / ".github/workflows/ci.yml").write_text(
        'python-version: ["3.10", "3.11"]\n'
    )
    (tmp_path / "docs/installation.rst").write_text(
        f"Python {python} or newer; tested Python 3.10–3.11.\nOptional: ``{extra}``.\n"
    )
    (tmp_path / "docs/cli_reference.rst").write_text(
        "\n".join(f"``{name}``" for name in scripts)
    )
    (tmp_path / "docs/architecture.rst").write_text(
        ".. code-block:: python\n\n   from pyar import Molecule\n"
    )
    (tmp_path / "README.md").write_text(
        f"Python 3.10–3.11.\n\n```bash\n{command}\n```\n"
    )
    return tmp_path


def fake_run(arguments, **_kwargs):
    if arguments[-1] == "--help":
        return SimpleNamespace(
            returncode=0, stdout="usage: command [--valid]\noptions: --valid", stderr=""
        )
    return SimpleNamespace(
        returncode=0, stdout="no calculations were performed", stderr=""
    )


def test_cli_documentation_rejects_removed_options(tmp_path):
    root = make_document_root(tmp_path, command="pyar react input.xyz --removed")
    with pytest.raises(ValueError, match="no documented option.*--removed"):
        check(root, run=fake_run)


def test_cli_documentation_rejects_unknown_extras(tmp_path):
    root = make_document_root(
        tmp_path, command='python -m pip install "pyar-chem[missing]"'
    )
    with pytest.raises(ValueError, match="undeclared extra \\[missing\\]"):
        check(root, run=fake_run)


def test_cli_documentation_rejects_stale_console_inventory(tmp_path):
    root = make_document_root(tmp_path, scripts=("pyar",))
    (root / "docs/cli_reference.rst").write_text("``pyar``\n``pyar-cli``\n")
    with pytest.raises(ValueError, match="CLI reference differs"):
        check(root, run=fake_run)


def test_cli_documentation_rejects_python_floor_drift(tmp_path):
    root = make_document_root(tmp_path, python="3.10")
    (root / "docs/installation.rst").write_text(
        "Python 3.11 or newer.\nOptional: ``demo``."
    )
    with pytest.raises(ValueError, match="Python support does not match"):
        check(root, run=fake_run)
