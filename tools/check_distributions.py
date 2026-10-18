"""Reject accidental runtime assets, caches and development files in archives."""

from pathlib import Path
import sys
from email.parser import BytesParser

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
import tarfile
import zipfile


def check(directory):
    wheels = list(Path(directory).glob('*.whl'))
    sdists = list(Path(directory).glob('*.tar.gz'))
    assert len(wheels) == len(sdists) == 1, 'Exactly one wheel and sdist required'
    project = tomllib.loads(Path('pyproject.toml').read_text())['project']
    assert project['version'] in Path('CHANGELOG.md').read_text(), 'Version missing from changelog'
    for path in [*wheels, *sdists]:
        wheel = path.suffix == '.whl'
        if wheel:
            with zipfile.ZipFile(path) as archive:
                names = archive.namelist()
                metadata_name = next(n for n in names if n.endswith('.dist-info/METADATA'))
                metadata = BytesParser().parsebytes(archive.read(metadata_name))
        else:
            with tarfile.open(path) as archive:
                names = [name.split('/', 1)[1] for name in archive.getnames() if '/' in name]
                metadata_name = next(n for n in archive.getnames() if n.count('/') == 1 and n.endswith('/PKG-INFO'))
                metadata = BytesParser().parsebytes(archive.extractfile(metadata_name).read())
        assert metadata['Name'] == project['name']
        assert metadata['Version'] == project['version']
        assert metadata['Requires-Python'] == project['requires-python']
        assert set(metadata.get_all('Provides-Extra', [])) == set(project['optional-dependencies'])
        assert metadata.get_all('Requires-Dist'), 'Missing dependency metadata'
        assert any('LICENSE' in name for name in names), 'Missing license'
        assert any('THIRD_PARTY' in name for name in names), 'Missing third-party notices'
        for name in names:
            assert '__pycache__' not in name and not name.endswith(('.pyc', '.pyo', '.jpt', '.pt')), name
            assert name not in {'setup.py', 'build_backend.py', 'requirements.txt'}, name
            assert not name.startswith(('pyar/mlatom/MLatomF', 'pyar/mlatom/cs.so')), name
            if wheel:
                assert not name.startswith(('tests/', 'tools/', 'benchmarks/', '.github/', 'docs/')), name
        required = {
            'pyar/dependencies.py', 'pyar/persistence/__init__.py',
            'pyar/growth/layout_v2.py', 'pyar/microsolvation/layout_v2.py',
            'pyar/aggregation/layout_v2.py', 'pyar/state/reaction_v2.py',
            'pyar/workflows/scan_layout.py', 'pyar/workflows/path_layout.py',
            'pyar/conformer/layout_v2.py',
        }
        assert required <= set(names), (path, sorted(required - set(names)))
        print(f'{path.name}: {len(names)} entries inspected; contents passed.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Archive checks require assertions enabled')
    check(sys.argv[1])
