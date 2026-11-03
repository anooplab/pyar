import pytest

from pyar.backends import turbomole


def test_unbiased_optimization_uses_instance_gamma_and_jobex(monkeypatch):
    subject = turbomole.Turbomole.__new__(turbomole.Turbomole)
    subject.basis = "def2-SVP"
    subject.method = "bp86"
    subject.atoms_list = ["H"]
    subject.start_coords = [[0.0, 0.0, 0.0]]
    subject.charge = 0
    subject.multiplicity = 2
    subject.gamma = None
    calls = []
    subject.run_turbomole_jobex = lambda **kwargs: calls.append(kwargs) or True
    monkeypatch.setattr(turbomole, "make_coord", lambda *args, **kwargs: None)
    monkeypatch.setattr(turbomole, "prepare_control", lambda **kwargs: True)

    result = subject.optimize()

    assert result is True
    assert calls == [{"max_cycles": 200}]
