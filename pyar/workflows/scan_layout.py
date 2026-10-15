"""Scan-owned immutable records; native scan and path executors retain filenames."""
from pathlib import Path
from uuid import uuid4

from pyar.persistence import ArtifactRef, RunLayout, read_json, write_json, sha256_file
from pyar.persistence.artifacts import validate_uuid


class ScanLayout:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.layout = RunLayout(self.root, 'scan-bond')
        marker = self.root / 'layout.json'
        self.run_id = read_json(marker)['run_id'] if marker.exists() else str(uuid4())

    def initialize(self, inputs):
        if (self.root / 'layout.json').exists():
            for ref in read_json(self.root / 'inputs/metadata.json')['inputs']:
                ArtifactRef.from_dict(ref).verify(self.root)
            return
        from pyar.growth.layout_v2 import copy_result
        refs = []
        for name, source in zip(('fragment_a', 'fragment_b'), inputs):
            relative = f'inputs/{name}.xyz'
            copy_result(Path(source), self.layout.path(relative))
            refs.append(ArtifactRef.capture(self.root, relative, 'input').to_dict())
        write_json(self.layout.path('inputs/metadata.json'), dict(inputs=refs), overwrite=False)
        write_json(self.root / 'layout.json', dict(format='pyar-run-layout', layout_version=2,
                   workflow='scan-bond', run_id=self.run_id), overwrite=False)

    def commit_path(self, index):
        return self.layout.path(f'orientation_{index:03d}/scan_commit.json')

    def _paths(self, value, restore=False):
        if isinstance(value, dict):
            if restore and set(value) == {'artifact_path'}:
                return str(self.layout.path(value['artifact_path']))
            return {key: self._paths(item, restore) for key, item in value.items()}
        if isinstance(value, list):
            return [self._paths(item, restore) for item in value]
        if not restore and isinstance(value, str) and value.startswith(str(self.root) + '/'):
            return dict(artifact_path=Path(value).relative_to(self.root).as_posix())
        return value

    def restore(self, index, request):
        path = self.commit_path(index)
        if not path.exists():
            return None
        record = read_json(path)
        validate_uuid(record['record_id'])
        if (record['orientation'] != index or record['run_id'] != self.run_id or record['request_signature'] != request['request_signature']
                or record['status'] != 'completed' or record['stage_result_schema_version'] != 1):
            raise ValueError('Scan commit identity does not match the request')
        for ref in record['artifacts']:
            ArtifactRef.from_dict(ref).verify(self.root)
        return self._paths(record['result'], restore=True)

    def commit(self, index, request, result, directory, execution):
        if self.commit_path(index).exists():
            self.restore(index, request)
            return
        paths = [directory / 'start.xyz']
        paths.extend(path for subtree in ('scan', 'ts_candidates')
                     for path in sorted((execution / subtree).rglob('*')) if path.is_file())
        refs = [ArtifactRef.capture(self.root, path.relative_to(self.root).as_posix(),
                'scan-maximum-candidate' if 'ts_candidates' in path.parts else 'scan-evidence').to_dict()
                for path in paths]
        write_json(self.commit_path(index), dict(stage_result_schema_version=1, record_id=str(uuid4()),
            run_id=self.run_id, orientation=index, status='completed', request_signature=request['request_signature'],
            physical_method=request['qc_params'], artifacts=refs, result=self._paths(dict(result, status='scan_success'))), overwrite=False)

    def relaxed(self, index, result=None, paths=()):
        path = self.layout.path(f'orientation_{index:03d}/relax_commit.json')
        if path.exists():
            record = read_json(path)
            validate_uuid(record['record_id'])
            if (record['run_id'] != self.run_id or record['orientation'] != index
                    or record['status'] != 'completed' or record['stage_result_schema_version'] != 1):
                raise ValueError('Invalid free-relaxation commit identity')
            if record['scan_commit_sha256'] != sha256_file(self.commit_path(index)):
                raise ValueError('Relaxation depends on a different scan commit')
            for ref in record['artifacts']:
                ArtifactRef.from_dict(ref).verify(self.root)
            return self._paths(record['result'], restore=True)
        if result is not None:
            refs = [ArtifactRef.capture(self.root, file.relative_to(self.root).as_posix(),
                                       'free-relaxation-evidence').to_dict() for file in paths]
            write_json(path, dict(stage_result_schema_version=1, record_id=str(uuid4()),
                run_id=self.run_id, orientation=index, status='completed', scan_commit_sha256=sha256_file(self.commit_path(index)),
                artifacts=refs, result=self._paths(result)), overwrite=False)
        return None

    def manifest(self, results):
        artifacts = []
        for result in results:
            path = self.commit_path(result['orientation'])
            if path.exists():
                record = read_json(path)
                for ref in record['artifacts']:
                    if ref['role'] == 'scan-maximum-candidate' and ref['path'].endswith('.xyz'):
                        artifacts.append(dict(ref, source_record_id=record['record_id'],
                                              validated_transition_state=False))
        write_json(self.layout.path('results/manifest.json'), dict(format='pyar-results-manifest',
            manifest_schema_version=1, workflow='scan-bond', run_id=self.run_id, artifacts=artifacts))
