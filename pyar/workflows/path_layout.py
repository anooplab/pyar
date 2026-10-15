"""Immutable stage records around the cumulative native reaction-path executor."""
from pathlib import Path
from uuid import uuid4

from pyar.persistence import ArtifactRef, RunLayout, read_json, write_json, sha256_file
from pyar.growth.layout_v2 import copy_result


class PathLayout:
    def __init__(self, root, start, end, guess, qc_params, options, dependency):
        self.root = Path(root).resolve()
        self.layout = RunLayout(self.root, 'scan-bond')
        inputs = []
        for name, source in zip(('reactant', 'product', 'ts_guess'), (start, end, guess)):
            relative = f'inputs/{name}.xyz'
            target = self.layout.path(relative)
            digest = sha256_file(source)
            if target.exists():
                if sha256_file(target) != digest:
                    raise ValueError('Reaction path input differs from its immutable snapshot')
            else:
                copy_result(Path(source), target)
            inputs.append(dict(path=relative, sha256=digest))
        metadata = self.layout.path('inputs/metadata.json')
        if not metadata.exists():
            write_json(metadata, dict(inputs=[ArtifactRef.capture(self.root, ref['path'], 'input').to_dict()
                                             for ref in inputs]), overwrite=False)
        for ref in read_json(metadata)['inputs']:
            ArtifactRef.from_dict(ref).verify(self.root)
        self.inputs = [self.layout.path(ref['path']) for ref in inputs]
        self.request = dict(inputs=inputs, physical_method=qc_params, scientific_options=options,
                            upstream_sha256=dependency)
        self.index = 0
        while True:
            attempt = self.layout.path(f'attempts/attempt_{self.index:03d}')
            path = attempt / 'request.json'
            if not path.exists():
                write_json(path, dict(request_schema_version=1, request=self.request), overwrite=False)
                break
            if read_json(path)['request'] == self.request:
                break
            self.index += 1
        self.execution = attempt / 'execution'

    def commit_path(self, stage):
        return self.layout.path(f'stages/{stage}/attempts/attempt_{self.index:03d}/commit.json')

    def verify(self, stage):
        path = self.commit_path(stage)
        if not path.exists():
            return
        record = read_json(path)
        if record['request'] != self.request or record['stage'] != stage or record['status'] != 'completed':
            raise ValueError('Reaction path stage identity mismatch')
        for ref in record['artifacts']:
            ArtifactRef.from_dict(ref).verify(self.root)

    def commit(self, stage, result):
        path = self.commit_path(stage)
        if path.exists():
            self.verify(stage)
            return
        summary = self.execution / f'{stage}_summary.json'
        if not summary.is_file():
            raise ValueError('Successful native stage did not publish its summary')
        paths = {summary}
        for key in ('inputs', 'artifacts'):
            paths.update(Path(value) for value in result.get(key, {}))
        for name in result.get('dependencies', []):
            paths.add(self.execution / f'{name}_summary.json')
        refs = [ArtifactRef.capture(self.root, file.relative_to(self.root).as_posix(), 'path-stage-evidence').to_dict()
                for file in sorted(paths)]
        from pyar.neb import _STAGE_GATES
        gates = {key: result.get(key) for key in _STAGE_GATES.get(stage, ())}
        write_json(path, dict(stage_result_schema_version=1, record_id=str(uuid4()), stage=stage,
            scientific_validation=gates,
            request=self.request, status='completed', artifacts=refs), overwrite=False)
