"""Planning contracts independent of backend execution and run directories."""
from copy import deepcopy
import pytest
from pyar.workflows.scan_path import _plan_characterisation
from pyar.workflows.scan_layout import _plan_registry_updates


@pytest.mark.parametrize('interpolation, refinement', [('linear', 'none'), ('idpp', 'idpp')])
def test_scan_seed_policy_preserves_representative_state(interpolation, refinement):
    definition = {'molecule': {'charge': -1, 'multiplicity': 2, 'scftype': 'uhf'},
                  'total_multiplicity': 2, 'total_multiplicity_source': 'explicit'}
    qc = {'software': 'xtb', 'charge': 0, 'multiplicity': 1}
    options = {'interpolation': interpolation, 'images': 7}
    original = deepcopy((definition, qc, options))
    plan = _plan_characterisation(definition, qc, options)
    assert plan.parameters == dict(qc, charge=-1, multiplicity=2, scftype='uhf',
        total_multiplicity=2, total_multiplicity_source='explicit')
    assert plan.options == dict(options, interpolation='linear', path_refinement=refinement)
    assert (definition, qc, options) == original


def test_scan_seed_policy_rejects_unsupported_route():
    with pytest.raises(ValueError, match='geodesic'):
        _plan_characterisation({'molecule': {}}, {}, {'interpolation': 'geodesic'})
    with pytest.raises(ValueError, match='absent'):
        _plan_characterisation(None, {}, {})


def test_registry_plan_keeps_unknown_unchanged_and_failed_separate():
    rows = [
        dict(orientation=4, status='scan_failed'),
        dict(orientation=3, status='scan_success_relax_success', identity_status='unknown'),
        dict(orientation=2, status='scan_success_relax_success', identity_status='success',
             relaxed_identity={'inchi':'InChI=1S/H2O/h1H2'}, product_identity_changed=False),
        dict(orientation=1, status='scan_success_relax_success', identity_status='success',
             relaxed_identity={'inchi':'InChI=1S/H2O/h1H2'}, product_identity_changed=True),
    ]
    sources = {row['orientation']:(deepcopy(row), {}) for row in rows if row['orientation'] != 4}
    before = deepcopy((rows, sources))
    plan = _plan_registry_updates(rows, sources)
    assert [item[0] for item in plan.groups[0]['members']] == [1]
    assert [item[0] for item in plan.unknown] == [3]
    assert plan.non_products == [{'orientation':2, 'reason':'identity_unchanged'}]
    assert plan.ineligible == [{'orientation':4, 'status':'scan_failed'}]
    assert (rows, sources) == before
