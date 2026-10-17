"""Microsolvation-owned v2 records; shell science remains in the workflow."""
import os
from pathlib import Path
from uuid import uuid4

from pyar.growth.layout_v2 import durable_pool, copy_result
from pyar.persistence import ArtifactRef, RunLayout, detect_version, read_json, write_json, sha256_file
from pyar.persistence.artifacts import validate_uuid
from pyar.persistence.json_store import fsync_directory
from pyar.state.aggregate import _json_value
from pyar.state.microsolvation import MicrosolvationRunState, MicrosolvationStateError


class MicrosolvationRunStateV2:
    def __init__(self, directory, data):
        self.directory = Path(directory).resolve()
        self.path = self.directory / 'state.json'
        self.data = data
        self.layout = RunLayout(self.directory, 'microsolvation')

    # Reuse only this workflow's established physical snapshot interpretation.
    restore_pool = MicrosolvationRunState.restore_pool
    validate_progress = MicrosolvationRunState.validate_progress

    def save(self):
        write_json(self.path, _json_value(self.data))

    def snapshot_pool(self, molecules, destination):
        refs = _json_value(MicrosolvationRunState.snapshot_pool(self, molecules, destination))
        for ref in refs:
            with (self.directory / ref['path']).open('rb') as stream:
                os.fsync(stream.fileno())
            ref.update(artifact_id=str(uuid4()), role='selected', energy_record=None if ref['energy'] is None
                       else dict(value=ref['energy'], unit='hartree', convention='physical_potential'))
        parent = Path(destination)
        while parent.is_relative_to(self.directory):
            fsync_directory(parent)
            if parent == self.directory:
                break
            parent = parent.parent
        return refs

    def verify_refs(self, refs):
        for ref in refs:
            ArtifactRef.from_dict({key: ref[key] for key in
                                  ('artifact_id', 'role', 'path', 'sha256')}).verify(self.directory)

    @classmethod
    def create(cls, directory, request, initial, solvent):
        directory.mkdir(parents=True, exist_ok=True)
        run_id = str(uuid4())
        inputs = durable_pool([initial, solvent], directory / 'inputs', directory, 'input')
        for index, name in enumerate(('solute.xyz', 'solvent.xyz')):
            source = directory / inputs[index]['path']
            target = directory / 'inputs' / name
            os.replace(source, target)
            inputs[index]['path'] = target.relative_to(directory).as_posix()
        fsync_directory(directory / 'inputs')
        state = cls(directory, dict(workflow_state_schema_version=2, workflow='microsolvation',
            run_id=run_id, request=request, status='running', next_step=1, completed_steps=[],
            solute_atom_indices=list(range(len(initial))), solvent_fragments=[], coverage_history=[],
            current_seeds=[], final_seeds=[], surface_sha256=None, selection_hashes=[]))
        refs = state.snapshot_pool([initial], state.layout.path('steps/step_000/selected'))
        state.data.update(initial_seeds=refs, current_seeds=refs)
        write_json(directory / 'inputs/metadata.json', dict(input_schema_version=1, run_id=run_id, inputs=inputs), overwrite=False)
        write_json(directory / 'request.json', dict(request_schema_version=1, run_id=run_id,
            workflow='microsolvation', scientific_request=request,
            inputs=[{key: ref[key] for key in ('artifact_id', 'role', 'path', 'sha256')} for ref in inputs]), overwrite=False)
        state.save()
        write_json(directory / 'layout.json', dict(format='pyar-run-layout', layout_version=2,
            run_id=run_id, workflow='microsolvation'), overwrite=False)
        return state

    @classmethod
    def load(cls, directory, request):
        directory = Path(directory).resolve()
        try:
            detect_version(directory, 'microsolvation')
            data = read_json(directory / 'state.json')
            identity = read_json(directory / 'request.json')
            if data['request'] != request or identity['scientific_request'] != request:
                raise MicrosolvationStateError('Microsolvation request differs from saved state')
            state = cls(directory, data)
            resampled_current = data.get('current_seeds')
            if data.get('resampling') is not None and data.get('completed_steps'):
                growth_current = state.selection(len(data['completed_steps']))['selected']
                data['current_seeds'] = growth_current
            state.validate_progress()
            data['current_seeds'] = resampled_current
            state.verify_refs(identity['inputs'])
            state.verify_refs(data['initial_seeds'])
            if len(data['initial_seeds']) != 1 or data['initial_seeds'][0]['sha256'] != identity['inputs'][0]['sha256']:
                raise MicrosolvationStateError('Initial solute differs from immutable input')
            steps = data['completed_steps']
            if len(steps) > request['count'] or type(data['next_step']) is not int:
                raise MicrosolvationStateError('Invalid microsolvation progress')
            if len(data['selection_hashes']) != len(steps):
                raise MicrosolvationStateError('Selection hash inventory differs from progress')
            if steps and data['surface_sha256'] is None:
                raise MicrosolvationStateError('Committed microsolvation steps lack original surface provenance')
            previous = None
            for index, item in enumerate(steps):
                record = state.selection(item['step'])
                if (sha256_file(state.selection_path(item['step'])) != data['selection_hashes'][index]
                        or type(item['step']) is not int
                        or type(item['selected_count']) is not int or item['selected_count'] != len(record['selected'])
                        or record['step_record'] != item
                        or record['upstream_sha256'] != previous):
                    raise MicrosolvationStateError('Invalid microsolvation selection chain')
                previous = sha256_file(state.selection_path(item['step']))
            expected = state.selection(len(steps))['selected'] if steps else data['initial_seeds']
            if (data['solute_atom_indices'] != list(range(len(request['solute']['atoms'])))
                    or (expected and data['solvent_fragments'] != expected[0]['solvent_fragments'])):
                raise MicrosolvationStateError('Original solute/solvent identity differs from committed pool')
            resampling = data.get('resampling')
            if resampling is None and data['current_seeds'] != expected:
                raise MicrosolvationStateError('Current pool differs from immutable selection')
            if resampling is not None:
                if 'resampling_policy' not in request:
                    raise MicrosolvationStateError('Unexpected MS4 state for a request without resampling')
                if resampling.get('schema_version') != 1 or resampling.get('current_references') != data['current_seeds']:
                    raise MicrosolvationStateError('Invalid resampling current-pool references')
                state.verify_refs(data['current_seeds'])
                attempts = resampling.get('attempts')
                if not isinstance(attempts, list) or len(attempts) > request.get('resampling_policy', {}).get('moves', 0):
                    raise MicrosolvationStateError('Invalid resampling attempt inventory')
                for attempt_index, attempt in enumerate(attempts):
                    if (attempt.get('attempt') != attempt_index
                            or attempt.get('status') not in {'running', 'interrupted', 'rejected_clash',
                                                              'rejected_identity', 'completed',
                                                              'backend_failure', 'budget_exhausted'}):
                        raise MicrosolvationStateError('Invalid resampling attempt sequence/status')
                    if attempt.get('proposal') is not None:
                        state.verify_refs([attempt['proposal']])
                    result_path = state.layout.path(
                        f"resampling/attempts/attempt_{attempt['attempt']:05d}/result.json")
                    if result_path.exists() and read_json(result_path).get('status') != attempt.get('status'):
                        raise MicrosolvationStateError('Resampling attempt result differs from saved state')
            if data['coverage_history'] != [item['coverage'] for item in steps]:
                raise MicrosolvationStateError('Coverage history differs from committed steps')
            if data.get('surface_sha256') is not None:
                if sha256_file(state.layout.path('surface/initial_surface.json')) != data['surface_sha256']:
                    raise MicrosolvationStateError('Original solute surface was modified')
            if data['status'] == 'completed':
                if len(steps) != request['count']:
                    raise MicrosolvationStateError('Completed microsolvation results are incomplete')
                if not data['final_seeds'] and 'resampling_policy' not in request:
                    raise MicrosolvationStateError('Completed microsolvation results are incomplete')
                if data['final_seeds']:
                    if len(data['final_seeds']) != len(data['current_seeds']):
                        raise MicrosolvationStateError('Completed microsolvation results are incomplete')
                    state.verify_refs(data['final_seeds'])
                    for index, (ref, source) in enumerate(zip(data['final_seeds'], data['current_seeds'])):
                        if (ref['sha256'] != source['sha256'] or ref['source_artifact_id'] != source['artifact_id']
                                or ref['path'] != f'results/structure_{index:03d}.xyz'
                                or ref['role'] != 'final-cluster' or type(ref['rank']) is not int or ref['rank'] != index
                                or any(ref[key] != source[key] for key in
                                       ('energy', 'energy_record', 'coverage', 'solute_atom_indices', 'solvent_fragments'))
                                or ref['source_record_id'] != (data.get('resampling', {}).get('record_id')
                                                               or state.selection(len(steps))['record_id'])):
                            raise MicrosolvationStateError('Final result differs from committed survivor')
            return state
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise MicrosolvationStateError(f'Invalid v2 microsolvation state: {exc}') from exc

    def publish_surface(self, payload):
        path = self.layout.path('surface/initial_surface.json')
        if path.exists():
            if read_json(path) != payload:
                raise MicrosolvationStateError('Original solute surface differs from request')
        else:
            write_json(path, payload, overwrite=False)
        digest = sha256_file(path)
        if self.data['surface_sha256'] not in (None, digest):
            raise MicrosolvationStateError('Original solute surface was modified')
        self.data['surface_sha256'] = digest
        self.save()

    def selection_path(self, step):
        return self.layout.path(f'steps/step_{step:03d}/selection.json')

    def selection(self, step):
        try:
            return self._read_selection(step)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise MicrosolvationStateError(f"Invalid microsolvation selection: {exc}") from exc

    def _read_selection(self, step):
        record = read_json(self.selection_path(step))
        validate_uuid(record['record_id'])
        if (type(record['stage_result_schema_version']) is not int or record['stage_result_schema_version'] != 1 or record['run_id'] != self.data['run_id']
                or type(record['step']) is not int or record['step'] != step
                or record['scientific_request'] != self.data['request'] or record['status'] != 'completed'
                or record['physical_method'] != self.data['request']['backend_parameters']
                or record['scientific_settings'] != self.data['request']
                or record['surface_sha256'] != self.data['surface_sha256']
                or record['selected'] != record['step_record']['selected_references']
                or len(record['selected']) > self.data['request']['maximum_number_of_seeds']):
            raise MicrosolvationStateError('Invalid microsolvation selection')
        for ref in record['selected']:
            if (ref['solute_atom_indices'] != list(range(len(self.data['request']['solute']['atoms'])))
                    or len(ref['solvent_fragments']) != step
                    or (self.data['request'].get('selection_policy', {}).get('objective', 'shell') == 'shell'
                        and not ref['coverage']['shell_valid'])):
                raise MicrosolvationStateError('Invalid original-solute shell survivor metadata')
        self.verify_refs(record['selected'])
        return record

    def complete_step(self, step, molecules, diagnostics):
        if step != self.data['next_step']:
            raise MicrosolvationStateError('Microsolvation step completion is out of sequence')
        refs = self.snapshot_pool(molecules, self.layout.path(f'steps/step_{step:03d}/selected'))
        diagnostics = _json_value(diagnostics)
        coverage = [dict(ref['coverage'], name=ref['name']) for ref in refs]
        step_record = dict(step=step, solvent_count=step, selected_count=len(refs),
                          selected_paths=[ref['path'] for ref in refs], selected_references=refs,
                          coverage=coverage, selection=diagnostics)
        write_json(self.layout.path(f'steps/step_{step:03d}/coverage.json'), step_record)
        write_json(self.selection_path(step), dict(stage_result_schema_version=1, record_id=str(uuid4()),
            run_id=self.data['run_id'], step=step, status='completed', stage='microsolvation-selection',
            scientific_request=self.data['request'], physical_method=self.data['request']['backend_parameters'],
            scientific_settings=self.data['request'], surface_sha256=self.data['surface_sha256'],
            validation=dict(shell_filtered=self.data['request'].get('selection_policy', {}).get(
                'objective', 'shell') == 'shell', survivor_budget=True), selected=refs, step_record=step_record,
            upstream_sha256=sha256_file(self.selection_path(step - 1)) if step > 1 else None), overwrite=False)
        self.accept_step(step)

    def accept_step(self, step):
        record = self.selection(step)
        expected = sha256_file(self.selection_path(step - 1)) if step > 1 else None
        if record['upstream_sha256'] != expected:
            raise MicrosolvationStateError('Stale microsolvation selection dependency')
        self.data['selection_hashes'].append(sha256_file(self.selection_path(step)))
        self.data['completed_steps'].append(record['step_record'])
        self.data['coverage_history'].append(record['step_record']['coverage'])
        self.data.update(current_seeds=record['selected'], next_step=step + 1)
        if record['selected']:
            self.data.update(solute_atom_indices=record['selected'][0]['solute_atom_indices'],
                             solvent_fragments=record['selected'][0]['solvent_fragments'])
        self.save()

    def finish_results(self):
        if not self.data['final_seeds']:
            refs = []
            for index, source in enumerate(self.data['current_seeds']):
                target = self.layout.path(f'results/structure_{index:03d}.xyz')
                copy_result(self.directory / source['path'], target)
                ref = dict(source, **ArtifactRef.capture(self.directory,
                           target.relative_to(self.directory).as_posix(), 'final-cluster').to_dict(),
                           source_artifact_id=source['artifact_id'],
                           source_record_id=(self.data.get('resampling', {}).get('record_id')
                                             or self.selection(len(self.data['completed_steps']))['record_id']), rank=index)
                refs.append(ref)
            self.data['final_seeds'] = refs
        self.data['status'] = 'completed'
        self.save()

    def report_results(self):
        write_json(self.layout.path('results/manifest.json'), dict(format='pyar-results-manifest',
            manifest_schema_version=1, workflow='microsolvation', run_id=self.data['run_id'],
            artifacts=self.data['final_seeds']))
        write_json(self.layout.path('summary.json'), dict(workflow='microsolvation', status=self.data['status'],
            requested_solvents=self.data['request']['count'], selected_count=len(self.data['final_seeds']),
            coverage=[ref['coverage'] for ref in self.data['final_seeds']]))
