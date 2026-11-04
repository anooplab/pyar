"""Legacy graph-to-MOL2 contracts; element inference is not chemical validation."""

import json
import builtins

import pytest

pytest.importorskip("rdkit")
from pyar import checkvalidity as validity


def convert(tmp_path, text, n):
    source = tmp_path / "graph.txt"
    output = tmp_path / "graph.mol2"
    source.write_text(text)
    before = source.read_bytes()
    result = validity.guess_correct_molecules(source, output, n, 1)
    assert source.read_bytes() == before
    return result, output


def mol2_graph(path):
    atoms, bonds = path.read_text().split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>BOND\n")
    atom_types = {
        int(row.split()[0]): row.split()[-1]
        for row in atoms.splitlines()
        if row.strip()
    }
    edges = {
        tuple(sorted(map(int, row.split()[1:3])))
        for row in bonds.splitlines()
        if row.strip()
    }
    return atom_types, edges


@pytest.mark.parametrize(
    "text",
    [
        "0 1 {'weight': 1}\n1 2 {'weight': 1}\n",
        "1 0 {'weight': 1}\n1 2 {'weight': 1}\n",
        "1 2 {'weight': 1}\n1 0 {'weight': 1}\n",
        "2 1 {'weight': 1}\n1 0 {'weight': 1}\n",
        "1 0 {1}\n1 2 {1}\n",
    ],
)
def test_atom_labels_and_bonds_do_not_depend_on_edge_insertion_order(tmp_path, text):
    result, output = convert(tmp_path, text, 3)
    assert result is True
    # Labels 0 and2 have degree1; label1 has degree2, independently of serialization.
    assert mol2_graph(output) == ({1: "H", 2: "O", 3: "H"}, {(1, 2), (2, 3)})
    assert output.read_text().splitlines()[2] == "3 2"


@pytest.mark.parametrize(
    "degree,element", [(1, "H"), (2, "O"), (3, "N"), (4, "C"), (5, "UN"), (6, "S")]
)
def test_existing_degree_to_element_heuristic_is_preserved(tmp_path, degree, element):
    text = "".join(f"0 {index} {{'weight': 1}}\n" for index in range(1, degree + 1))
    result, output = convert(tmp_path, text, degree + 1)
    assert result is True
    atom_types, edges = mol2_graph(output)
    assert atom_types[1] == element
    assert all(atom_types[index] == "H" for index in range(2, degree + 2))
    assert edges == {(1, index) for index in range(2, degree + 2)}


def test_multiple_bond_weights_and_padding_preserve_output_schema(tmp_path):
    result, output = convert(tmp_path, "0 1 {'weight': 3}\n0 2 {'weight': 1}\n", 5)
    assert result is True
    assert mol2_graph(output) == ({1: "C", 2: "N", 3: "H"}, {(1, 2), (1, 3)})
    assert output.read_text() == (
        "@<TRIPOS>MOLECULE\nDummy Atom\n3 2\n\n@<TRIPOS>ATOM\n"
        "1\tC1\t0\t0\t0\tC\n2\tN1\t0\t0\t0\tN\n3\tH1\t0\t0\t0\tH\n"
        "\n@<TRIPOS>BOND\n1\t1\t2\t3\n2\t1\t3\t1\n"
    )


@pytest.mark.parametrize(
    "text,n",
    [
        ("0 1 {'weight': 1}\n", 0),
        ("0 1 {'weight': 1}\n", 1),
        ("-1 0 {'weight': 1}\n", 2),
        ("0 7 {'weight': 1}\n", 3),
    ],
)
def test_atom_labels_outside_requested_count_are_rejected_before_output(
    tmp_path, text, n
):
    output = tmp_path / "graph.mol2"
    output.write_bytes(b"previous output\n")
    result, _ = convert(tmp_path, text, n)
    assert result is False
    assert output.read_bytes() == b"previous output\n"


@pytest.mark.parametrize(
    "text,n,message",
    [
        ("", 3, "Null graph"),
        ("0 1 {'weight': 1}\n2 3 {'weight': 1}\n", 4, "Not connected"),
        ("0 1 {'weight': 7}\n", 2, None),
    ],
)
def test_rejected_graph_does_not_overwrite_existing_output(
    tmp_path, capsys, text, n, message
):
    output = tmp_path / "graph.mol2"
    output.write_bytes(b"previous output\n")
    result, path = convert(tmp_path, text, n)
    assert result is False
    assert path.read_bytes() == b"previous output\n"
    if message:
        assert message in capsys.readouterr().out


@pytest.mark.parametrize("operation", ["read_edgelist", "is_connected"])
def test_unexpected_graph_errors_are_not_converted_into_a_valid_or_null_graph(
    tmp_path, monkeypatch, operation
):
    error = RuntimeError("injected graph implementation failure")

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(validity.nx, operation, fail)
    with pytest.raises(RuntimeError) as caught:
        convert(tmp_path, "0 1 {'weight': 1}\n", 2)
    assert caught.value is error
    assert not (tmp_path / "graph.mol2").exists()


def test_valence_dictionary_reads_json_without_modifying_source(tmp_path):
    path = tmp_path / "valence.json"
    payload = {"C": [4], "O": [2]}
    path.write_text(json.dumps(payload))
    before = path.read_bytes()
    assert validity.getValenceDict(path) == payload
    assert path.read_bytes() == before


def test_valence_dictionary_preserves_json_parse_error(tmp_path):
    path = tmp_path / "invalid.json"
    path.write_text("{")
    with pytest.raises(json.JSONDecodeError):
        validity.getValenceDict(path)


def test_graph_conversion_closes_owned_streams_even_when_references_are_retained(
    tmp_path, monkeypatch
):
    handles = []

    def track(*args, **kwargs):
        stream = builtins.open(*args, **kwargs)
        handles.append(stream)
        return stream

    monkeypatch.setattr(validity, "open", track, raising=False)
    result, _ = convert(tmp_path, "0 1 {1}\n", 2)
    assert result is True
    assert len(handles) == 2
    assert all(stream.closed for stream in handles)
