"""Reject accidental runtime assets, caches and development files in archives."""

from pathlib import Path
import sys
import tarfile
import zipfile


def check(directory):
    wheels = list(Path(directory).glob('*.whl'))
    sdists = list(Path(directory).glob('*.tar.gz'))
    assert wheels and sdists, 'Both wheel and sdist are required'
    for path in [*wheels, *sdists]:
        wheel = path.suffix == '.whl'
        if wheel:
            with zipfile.ZipFile(path) as archive:
                names = archive.namelist()
        else:
            with tarfile.open(path) as archive:
                names = [name.split('/', 1)[1] for name in archive.getnames() if '/' in name]
        for name in names:
            assert '__pycache__' not in name and not name.endswith(('.pyc', '.pyo', '.jpt', '.pt')), name
            assert name not in {'setup.py', 'build_backend.py', 'requirements.txt'}, name
            assert not name.startswith(('pyar/mlatom/MLatomF', 'pyar/mlatom/cs.so')), name
            if wheel:
                assert not name.startswith(('tests/', 'tools/', 'benchmarks/', '.github/', 'docs/')), name
        assert any(name.endswith('pyar/dependencies.py') for name in names), path
        print(f'{path.name}: {len(names)} entries inspected; contents passed.')


if __name__ == '__main__':
    check(sys.argv[1])
