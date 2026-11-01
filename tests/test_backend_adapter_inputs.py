"""Input generation and status handling for the legacy QC adapters, without real executables."""
import numpy as np
import pytest

from pyar.core.molecule import Molecule


def _water(**state):
    return Molecule(['O', 'H', 'H'], [[0, 0, 0], [.758, 0, .504], [-.758, 0, .504]], name='w', **state)


def test_native_gaussian_requests_optimization_and_feeds_stdin(tmp_path, monkeypatch):
    from pyar.backends import gaussian
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gaussian, 'require_executable', lambda *args: 'g16')
    calls = []
    monkeypatch.setattr(gaussian, 'run_command', lambda command, **kwargs: calls.append((command, kwargs)) or 1)
    job = gaussian.Gaussian(_water(), dict(basis='def2-SVP', method='B3LYP', nprocs=1, scf_cycles=64, opt_cycles=30))
    assert job.optimize() is False
    route = next(line for line in (tmp_path / 'trial_w.com').read_text().splitlines() if line.startswith('#'))
    assert ' Opt=(MaxCycles=30) ' in route
    command, kwargs = calls[0]
    assert command == ['g16']
    assert kwargs['stdin_path'] == 'trial_w.com' and kwargs['stdout_path'] == 'trial_w.log'


def test_turbomole_refuses_unsupported_methods():
    from pyar.backends.turbomole import turbomole_functional
    assert turbomole_functional('BP86') == 'bp'
    assert turbomole_functional('b3lyp') == 'b3-lyp'
    with pytest.raises(ValueError, match='supports only'):
        turbomole_functional('pbe0')


def test_turbomole_define_requests_unpaired_electrons(tmp_path, monkeypatch):
    from pyar.backends import turbomole
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'coord').write_text('$coord\n$end\n')

    def fake_define(command, stdin, stdout, **kwargs):
        stdout.write('define : all done\n')
    monkeypatch.setattr(turbomole.subp, 'check_call', fake_define)
    assert turbomole.prepare_control(multiplicity=3) is True
    assert 'u 2\n' in (tmp_path / 'define.inp').read_text()


def test_turbomole_afir_gradient_header_matches_written_rows(tmp_path, monkeypatch):
    from pyar.backends import turbomole
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'gradient').write_text(
        '$grad\n'
        '  cycle =      1    SCF energy =     -1.0000000000   |dE/dxyz| =  0.100000\n'
        '    0.0 0.0 0.0 h\n'
        '    0.0 0.0 1.4 h\n'
        '   0.10D+00   0.0D+00   0.0D+00\n'
        '  -0.10D+00   0.0D+00   0.0D+00\n'
        '$end\n')
    force = np.array([[0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]])
    turbomole.rewrite_turbomole_gradient(2, 0.0, force)
    lines = (tmp_path / 'gradient').read_text().splitlines()
    header = next(line for line in lines if 'cycle' in line)
    rows = np.array([[float(value.replace('D', 'E')) for value in line.split()] for line in lines[4:6]])
    np.testing.assert_allclose(rows, [[0.05, 0, 0], [-0.05, 0, 0]])
    assert float(header.split()[-1]) == pytest.approx(np.linalg.norm(rows), abs=1e-6)


def test_psi4_uses_settings_and_valid_open_shell_reference(tmp_path, monkeypatch):
    from pyar.backends import psi4
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(psi4, 'require_executable', lambda *args: 'psi4')
    radical = Molecule(['O', 'H'], [[0, 0, 0], [0, 0, .97]], name='oh', multiplicity=2)
    psi4.Psi4(radical, dict(method='PBE0', basis='def2-TZVP'))
    text = (tmp_path / 'trial_oh.in').read_text()
    assert 'set basis def2-TZVP\nset reference uhf\n' in text
    assert 'optimize("PBE0")' in text
    assert 'ECP' not in text and 'UKS' not in text
    psi4.Psi4(_water(), {})
    text = (tmp_path / 'trial_w.in').read_text()
    assert 'set basis def2-SVP\n' in text and 'reference' not in text and 'optimize("B97-D")' in text


@pytest.mark.parametrize('method', ['HF', 'PBE0', 'MP2', 'MP3', 'CCSD'])
def test_psi4_open_shell_reference_supports_dft_and_correlated_methods(tmp_path, monkeypatch, method):
    from pyar.backends import psi4
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(psi4, 'require_executable', lambda *args: 'psi4')
    radical = Molecule(['O', 'H'], [[0, 0, 0], [0, 0, .97]], name='oh', multiplicity=2)
    psi4.Psi4(radical, dict(method=method, basis='def2-SVP'))
    text = (tmp_path / 'trial_oh.in').read_text()
    assert 'set reference uhf\n' in text
    assert f'optimize("{method}")' in text


def test_mopac_spin_keywords():
    from pyar.backends.mopac import mopac_spin_keywords
    assert mopac_spin_keywords(1) == ''
    assert mopac_spin_keywords(2) == ' UHF DOUBLET'
    assert mopac_spin_keywords(3) == ' UHF TRIPLET'
    with pytest.raises(ValueError):
        mopac_spin_keywords(9)


def test_orca_aiqm1_requires_a_user_supplied_wrapper(tmp_path, monkeypatch):
    from pyar.backends import orca_aiqm1
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(orca_aiqm1, 'require_executable', lambda *args: 'orca')
    monkeypatch.delenv(orca_aiqm1.EXTERNAL_PROGRAM_ENV, raising=False)
    with pytest.raises(ValueError, match=orca_aiqm1.EXTERNAL_PROGRAM_ENV):
        orca_aiqm1.OrcaAIQM1(_water(), {})
    wrapper = tmp_path / 'combine.py'
    wrapper.write_text('')
    monkeypatch.setenv(orca_aiqm1.EXTERNAL_PROGRAM_ENV, str(wrapper))
    job = orca_aiqm1.OrcaAIQM1(_water(), {})
    assert f'ProgExt "{wrapper}"' in job.keyword
    assert '/scratch/' not in job.keyword
