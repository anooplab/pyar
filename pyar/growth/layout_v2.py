"""Grow-owned v2 records; legacy state and execution files remain independent."""
import os
from pathlib import Path
import shutil
import tempfile
from uuid import uuid4

from pyar.growth.reporting import snapshot_pool, restore_pool
from pyar.persistence import ArtifactRef, RunLayout, detect_version, read_json, write_json, sha256_file
from pyar.persistence.json_store import fsync_directory
from pyar.persistence.artifacts import validate_uuid
from pyar.state.grow import GrowStateError


def durable_pool(molecules, destination, root, role="selected"):
    refs = snapshot_pool(molecules, destination, root)
    for ref in refs:
        path = root / ref['path']
        with path.open('rb') as stream:
            os.fsync(stream.fileno())
        ref.update(artifact_id=str(uuid4()), role=role,
                   energy_record=None if ref['energy'] is None else
                   dict(value=ref['energy'], unit='hartree', convention='physical_potential'))
    parent = Path(destination)
    while parent.is_relative_to(root):
        fsync_directory(parent)
        if parent == root:
            break
        parent = parent.parent
    return refs


def verify_refs(refs, root):
    for ref in refs:
        try:
            ArtifactRef.from_dict({key: ref[key] for key in
                                  ('artifact_id', 'role', 'path', 'sha256')}).verify(root)
        except ValueError as exc:
            raise GrowStateError(f'Growth snapshot was modified or invalid: {exc}') from exc
    return restore_pool(refs, root)


def copy_result(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            with source.open('rb') as original:
                shutil.copyfileobj(original, stream)
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, destination)
        fsync_directory(destination.parent)
        fsync_directory(destination.parent.parent)
    finally:
        temporary.unlink(missing_ok=True)


class GrowRunStateV2:
    def __init__(self, directory, data):
        self.directory = Path(directory).resolve()
        self.path = self.directory / 'state.json'
        self.data = data

    @classmethod
    def load(cls, directory, request):
        directory = Path(directory).resolve()
        try:
            detect_version(directory, 'grow')
            data = read_json(directory / 'state.json')
            identity = read_json(directory / 'request.json')
            if data.get('request') != request or identity.get('scientific_request') != request:
                raise GrowStateError('Growth request differs from saved state; use a fresh output directory')
            for ref in identity['inputs']:
                ArtifactRef.from_dict(ref).verify(directory)
            if (len(data['initial_seeds']) != 1
                    or data['initial_seeds'][0]['sha256'] != identity['inputs'][0]['sha256']):
                raise GrowStateError('Initial pool differs from immutable input seed')
            state = cls(directory, data)
            steps = data['completed_steps']
            if (type(data['next_step']) is not int or data['next_step'] != len(steps) + 1
                    or [item['step'] for item in steps] != list(range(1, len(steps) + 1))
                    or len(steps) > request['count']
                    or data['status'] not in {'running', 'stopped', 'completed', 'no_candidates'}
                    or (data['status'] == 'completed' and len(steps) != request['count'])):
                raise GrowStateError('Invalid growth progress in saved state')
            previous = None
            for item in steps:
                if type(item['step']) is not int or type(item['selected_count']) is not int:
                    raise GrowStateError('Invalid growth step/count types')
                record = state.selection(item['step'])
                if (sha256_file(state.selection_path(item['step'])) != item['selection_sha256']
                        or record['upstream_sha256'] != previous
                        or item['selected_paths'] != [ref['path'] for ref in record['selected']]
                        or item['selected_count'] != len(record['selected'])):
                    raise GrowStateError('Invalid committed growth selection chain')
                previous = item['selection_sha256']
            expected = state.selection(len(steps))['selected'] if steps else data['initial_seeds']
            if data['current_seeds'] != expected:
                raise GrowStateError('Growth current pool does not match committed selection')
            verify_refs(data['initial_seeds'], directory)
            verify_refs(data['current_seeds'], directory)
            if data.get('final_refs') and (len(steps) != request['count']
                    or len(data['final_refs']) != len(data['current_seeds'])):
                raise GrowStateError('Incomplete final result pool')
            for index, ref in enumerate(data.get('final_refs', [])):
                source = data['current_seeds'][index]
                if (ref['sha256'] != source['sha256']
                        or ref['source_artifact_id'] != source['artifact_id']
                        or ref['rank'] != index
                        or ref['source_record_id'] != state.selection(len(steps))['record_id']):
                    raise GrowStateError('Final result provenance differs from committed selection')
                ArtifactRef.from_dict({key: ref[key] for key in
                                      ('artifact_id', 'role', 'path', 'sha256')}).verify(directory)
            if data['status'] == 'completed' and len(data.get('final_refs', [])) != len(data['current_seeds']):
                raise GrowStateError('Completed growth results are missing')
            return state
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise GrowStateError(f'Invalid v2 growth state: {exc}') from exc

    @classmethod
    def create(cls, directory, request, seed, monomer):
        directory.mkdir(parents=True, exist_ok=True)
        run_id = str(uuid4())
        inputs = durable_pool([seed, monomer], directory / 'inputs', directory, 'input')
        for index, name in enumerate(('seed.xyz', 'monomer.xyz')):
            source = directory / inputs[index]['path']
            target = directory / 'inputs' / name
            os.replace(source, target)
            inputs[index]['path'] = target.relative_to(directory).as_posix()
        fsync_directory(directory / 'inputs')
        initial = durable_pool([seed], directory / 'steps/step_000/selected', directory)
        data = dict(workflow_state_schema_version=2, workflow='grow', run_id=run_id,
                    request=request, status='running', next_step=1, completed_steps=[],
                    initial_seeds=initial, current_seeds=initial, final_refs=[])
        state = cls(directory, data)
        write_json(directory / 'inputs/metadata.json', dict(input_schema_version=1, run_id=run_id,
                    inputs=inputs), overwrite=False)
        write_json(directory / 'request.json', dict(request_schema_version=1, run_id=run_id,
                    workflow='grow', scientific_request=request,
                    inputs=[{key: ref[key] for key in ('artifact_id', 'role', 'path', 'sha256')} for ref in inputs]), overwrite=False)
        state.save()
        # Marker published last: partial unmarked initialization is rejected,
        # never mistaken for legacy state because it has no v1 version field.
        write_json(directory / 'layout.json', dict(format='pyar-run-layout', layout_version=2,
                    workflow='grow', run_id=run_id), overwrite=False)
        return state

    def selection_path(self, step):
        return RunLayout(self.directory, 'grow').path(f'steps/step_{step:03d}/selection.json')

    def selection(self, step):
        try:
            record = read_json(self.selection_path(step))
            validate_uuid(record['record_id'])
            if (record['run_id'] != self.data['run_id'] or type(record['step']) is not int or record['step'] != step
                    or record['scientific_request'] != self.data['request']
                    or record['physical_method'] != self.data['request']['backend_parameters']
                    or record['scientific_settings'] != self.data['request']
                    or record['status'] != 'completed' or record['stage_result_schema_version'] != 1
                    or len(record['selected']) > self.data['request']['maximum_number_of_seeds']):
                raise GrowStateError('Invalid growth selection record')
            verify_refs(record['selected'], self.directory)
            return record
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise GrowStateError(f'Invalid growth selection: {exc}') from exc

    def commit_step(self, step, refs, selection_mode):
        path = self.selection_path(step)
        record = dict(stage_result_schema_version=1, record_id=str(uuid4()),
                      run_id=self.data['run_id'], step=step, stage='grow-selection', status='completed',
                      scientific_request=self.data['request'], selected=refs, selection_mode=selection_mode,
                      physical_method=self.data['request']['backend_parameters'],
                      scientific_settings=self.data['request'], validation=dict(survivor_budget=True),
                      upstream_sha256=self.data['completed_steps'][-1]['selection_sha256']
                      if self.data['completed_steps'] else None)
        write_json(path, record, overwrite=False)
        self.accept_step(step)

    def accept_step(self, step):
        record = self.selection(step)
        expected = self.data['completed_steps'][-1]['selection_sha256'] if self.data['completed_steps'] else None
        if record['upstream_sha256'] != expected:
            raise GrowStateError('Stale growth selection dependency')
        self.data['completed_steps'].append(dict(step=step, selected_count=len(record['selected']),
                    selected_paths=[ref['path'] for ref in record['selected']],
                    selection_mode=record['selection_mode'], selection_sha256=sha256_file(self.selection_path(step))))
        self.data.update(current_seeds=record['selected'], next_step=step + 1)
        if not record['selected']:
            self.data['status'] = 'no_candidates'
        self.save()

    def finish(self):
        if not self.data.get('final_refs'):
            refs = []
            for index, source in enumerate(self.data['current_seeds']):
                target = RunLayout(self.directory, 'grow').path(f'results/structure_{index:03d}.xyz')
                copy_result(self.directory / source['path'], target)
                ref = ArtifactRef.capture(self.directory, target.relative_to(self.directory).as_posix(),
                                          'final-structure').to_dict()
                ref.update(source_artifact_id=source['artifact_id'], rank=index,
                           source_record_id=self.selection(len(self.data['completed_steps']))['record_id'],
                           energy_record=source['energy_record'])
                refs.append(ref)
            self.data['final_refs'] = refs
        self.data['status'] = 'completed'
        self.save()

    def save(self):
        write_json(self.path, self.data)
