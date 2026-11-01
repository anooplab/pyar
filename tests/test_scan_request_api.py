"""Invocation mappings do not replace signed scan/restart authority."""
from copy import deepcopy
from pathlib import Path

import pytest

from pyar.backends import write_xyz
from pyar.core.molecule import Molecule
from pyar.workflows import scan_bond


def invocation(tmp_path):
    first,second = tmp_path/'a.xyz',tmp_path/'b.xyz'
    for path in (first,second):
        write_xyz(['H'], [[0,0,0]], path)
    return dict(input_a=first,input_b=second,atoms=(0,0),orientations=1,
                qc_params={'software':'xtb'},output_dir=tmp_path/'scan',scan_end=1.,
                scan_points=3,through='scan')


def test_request_api_owns_a_copy_and_reuses_signed_layout(tmp_path, monkeypatch, deterministic_workflow_backend):
    request = invocation(tmp_path)
    before = deepcopy(request)
    molecule = Molecule(['H','H'], [[0,0,0],[2,0,0]], name='orientation_000', fragments=[[0],[1]])
    monkeypatch.setattr(scan_bond,'_make_orientations',lambda *args:[(0,molecule,None)])
    result = scan_bond.run_scan_bond_request(request)
    assert result['status'] == 'complete'
    assert request == before
    signed = Path(request['output_dir'])/'request.json'
    original = signed.read_bytes()
    repeated = scan_bond.run_scan_bond_request(request)
    assert repeated['status'] == 'complete'
    assert signed.read_bytes() == original
    assert len([c for c in deterministic_workflow_backend.calls if c['operation']=='scan']) == 1


def test_request_api_rejects_manifest_fields_before_allocating_run(tmp_path):
    request = invocation(tmp_path)
    request['request_signature'] = 'not-an-invocation-field'
    with pytest.raises(TypeError, match='request_signature'):
        scan_bond.run_scan_bond_request(request)
    assert not Path(request['output_dir']).exists()
