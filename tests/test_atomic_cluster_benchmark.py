import json
import tarfile
import zipfile

import numpy as np

from pyar.scripts.atomic_cluster_benchmark import (
    _unwrap_nomad_record,
    load_lj13_archive,
    load_qcd_au13,
    pair_distance_fingerprint,
)


def test_load_lj13_reduced_units_and_structure_alignment(tmp_path):
    archive_path = tmp_path / "lj13.tar.bz2"
    with tarfile.open(archive_path, "w:bz2") as archive:
        for name, payload in {
            "min.data": " -44.0 0 0\n -40.0 0 0\n",
            "extractedmin": "".join(
                f"{i + frame * 20} 0 0\n" for frame in range(2) for i in range(13)
            ),
        }.items():
            path = tmp_path / name
            path.write_text(payload)
            archive.add(path, arcname=name)

    molecules, metadata = load_lj13_archive(archive_path)
    assert len(molecules) == len(metadata) == 2
    assert len(molecules[0].atoms_list) == 13
    assert molecules[1].coordinates[0, 0] == 20
    assert molecules[0].energy == -44.0
    assert metadata[0]["coordinate_unit"] == "sigma"
    assert metadata[0]["energy_unit"] == "epsilon"


def test_unwrap_nomad_record_converts_si_units_and_checks_au13():
    record = {
        "entry_id": "entry-1",
        "archive": {
            "metadata": {"entry_id": "entry-1", "license": "CC BY 4.0"},
            "run": [{
                "system": [{"atoms": {
                    "labels": ["Au"] * 13,
                    "positions": [[1e-10, 2e-10, 3e-10]] * 13,
                }}],
                "calculation": [{"energy": {"total": {"value": -2 * 1.602176634e-19}}}],
            }],
        },
    }
    metadata, symbols, positions, energy = _unwrap_nomad_record(record)
    assert metadata["license"] == "CC BY 4.0"
    assert symbols == ["Au"] * 13
    np.testing.assert_allclose(positions[0], [1.0, 2.0, 3.0])
    assert energy == -2.0


def test_load_qcd_au13_skips_manifest_and_preserves_entry_metadata(tmp_path):
    zip_path = tmp_path / "qcd.zip"
    entry = {
        "entry_id": "entry-1",
        "archive": {
            "metadata": {"entry_id": "entry-1", "mainfile": "Au/13/x.xml",
                         "upload_id": "upload", "license": "CC BY 4.0"},
            "run": [{
                "system": [{"atoms": {"labels": ["Au"] * 13,
                                        "positions": [[0.0, 0.0, 0.0]] * 13}}],
                "calculation": [{"energy": {"total": {"value": -1e-18}}}],
            }],
        },
    }
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("upload/entry.json", json.dumps(entry))
        archive.writestr("manifest.json", "[]")
    molecules, metadata = load_qcd_au13(zip_path)
    assert len(molecules) == len(metadata) == 1
    assert metadata[0]["source_id"] == "entry-1"
    assert metadata[0]["mainfile"] == "Au/13/x.xml"
    assert molecules[0].energy < 0


def test_pair_distance_fingerprint_ignores_rigid_motion_and_order():
    from pyar.core.molecule import Molecule

    xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    rotated_translated = xyz[[2, 0, 1]] @ np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    ) + [8.0, -4.0, 2.0]
    first = Molecule(["Au"] * 3, xyz)
    second = Molecule(["Au"] * 3, rotated_translated)
    np.testing.assert_allclose(pair_distance_fingerprint(first),
                               pair_distance_fingerprint(second))
