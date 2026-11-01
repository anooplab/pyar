"""Old-generation data must fail before it can authorize execution or mutation."""
import json
import pytest
from pyar.persistence import detect_version, LayoutError
from pyar.run_config import load_toml_document, RUN_CONFIG_KIND


def test_old_layout_marker_is_rejected_without_changing_bytes(tmp_path):
    marker = tmp_path / 'layout.json'
    marker.write_text(json.dumps(dict(format='pyar-run-layout', layout_version=2, workflow='grow')))
    before = marker.read_bytes()
    with pytest.raises(LayoutError):
        detect_version(tmp_path, 'grow')
    assert marker.read_bytes() == before
    assert list(tmp_path.iterdir()) == [marker]


def test_old_config_is_rejected_read_only(tmp_path):
    path = tmp_path / 'old.toml'
    path.write_text('schema_version = 1\nkind = "pyar-run-config"\nworkflow = "optimize"\n')
    before = path.read_bytes()
    with pytest.raises(ValueError):
        load_toml_document(path, RUN_CONFIG_KIND)
    assert path.read_bytes() == before
    assert list(tmp_path.iterdir()) == [path]
