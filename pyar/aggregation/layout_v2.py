"""Aggregation-owned records; logical legacy names preserve final selection order."""
import os
import string
from collections import Counter
from pathlib import Path
from uuid import uuid4

from pyar.growth.layout_v2 import durable_pool, verify_refs as _verify_growth_refs, copy_result
from pyar.state.grow import GrowStateError
from pyar.persistence import RunLayout, ArtifactRef, read_json, write_json, sha256_file, detect_version
from pyar.persistence.artifacts import validate_uuid
from pyar.persistence.json_store import fsync_directory
from pyar.state.aggregate import AggregateStateError, _json_value


def verify_refs(refs, root):
    """Keep shared snapshot validation within the aggregation error contract."""
    try:
        return _verify_growth_refs(refs, root)
    except GrowStateError as exc:
        raise AggregateStateError(f'Invalid aggregation snapshot: {exc}') from exc


class AggregateRunStateV2:
    def __init__(self, root, data):
        self.root_directory = Path(root).resolve()
        self.aggregate_directory = self.root_directory / 'aggregates'
        self.state_file = self.aggregate_directory / 'state.json'
        self.layout = RunLayout(self.aggregate_directory, 'aggregate')
        self.data = data

    @classmethod
    def create(cls, root, request, molecules, labels, sampling):
        directory = Path(root).resolve() / 'aggregates'
        directory.mkdir(parents=True, exist_ok=True)
        run_id = str(uuid4())
        inputs = durable_pool(molecules, directory / 'inputs', directory, 'input')
        for index, ref in enumerate(inputs):
            target = directory / f'inputs/component_{index:03d}.xyz'
            os.replace(directory / ref['path'], target)
            ref['path'] = target.relative_to(directory).as_posix()
        fsync_directory(directory / 'inputs')
        state = cls(root, dict(workflow_state_schema_version=2, workflow='aggregate', run_id=run_id,
            request=request, status='running', pathways=list(labels), completed_pathways=[],
            sampling=sampling, steps={}, pathway_records={}, final_selected_results=[]))
        write_json(directory / 'request.json', dict(request_schema_version=1, workflow='aggregate',
            run_id=run_id, scientific_request=request, inputs=[{key: ref[key] for key in
            ('artifact_id','role','path','sha256')} for ref in inputs]), overwrite=False)
        write_json(directory / 'inputs/metadata.json', dict(input_schema_version=1, run_id=run_id, inputs=inputs), overwrite=False)
        state.initialize_pathways()
        write_json(directory / 'layout.json', dict(format='pyar-run-layout', layout_version=2,
            workflow='aggregate', run_id=run_id), overwrite=False)
        return state

    @classmethod
    def load(cls, root, request):
        directory = Path(root).resolve() / 'aggregates'
        try:
            detect_version(directory, 'aggregate')
            data = read_json(directory / 'state.json')
            identity = read_json(directory / 'request.json')
            if data['request'] != request or identity['scientific_request'] != request:
                raise AggregateStateError('Existing aggregation state does not match this invocation')
            state = cls(root, data)
            if data['status'] != 'running':
                raise AggregateStateError(f"Aggregation state is already {data['status']!r}; start a new calculation")
            for ref in identity['inputs']:
                ArtifactRef.from_dict(ref).verify(directory)
            for index, label in enumerate(data['pathways']):
                meta = state.pathway(index, label)
                recorded = data['pathway_records'][str(index)]
                if sha256_file(meta) != recorded:
                    raise AggregateStateError('Aggregation pathway identity was modified')
                for step, digest in data['steps'].get(str(index), {}).items():
                    record = state.selection(index, int(step))
                    if sha256_file(state.selection_path(index, int(step))) != digest:
                        raise AggregateStateError('Committed aggregation selection was modified')
                    expected = data['steps'].get(str(index), {}).get(str(int(step) - 1)) if int(step) else None
                    if record['upstream_sha256'] != expected:
                        raise AggregateStateError('Stale aggregation step dependency')
                    state.verify_selected(record)
            completed = data['completed_pathways']
            if [entry['index'] for entry in completed] != list(range(len(completed))):
                raise AggregateStateError('Aggregation completed pathways are out of order')
            for entry in completed:
                if type(entry['index']) is not int or entry['label'] != data['pathways'][entry['index']]:
                    raise AggregateStateError('Aggregation pathway label mismatch')
            if state.layout.path('final_selection.json').exists():
                state.final_record()
            return state
        except (OSError, ValueError, TypeError, KeyError) as exc:
            raise AggregateStateError(f'Invalid v2 aggregation state: {exc}') from exc

    def save(self):
        write_json(self.state_file, _json_value(self.data))

    def pathway(self, index, label):
        path = self.layout.path(f'pathways/pathway_{index:03d}/pathway.json')
        record = read_json(path)
        validate_uuid(record['pathway_id'])
        inputs = read_json(self.layout.path('request.json'))['inputs']
        sequence = [inputs[string.ascii_lowercase.index(char)]['artifact_id'] for char in label]
        if (record['run_id'] != self.data['run_id'] or record['label'] != label
                or type(record['index']) is not int or record['index'] != index
                or record['scientific_request'] != self.data['request']
                or record['addition_sequence'] != sequence):
            raise AggregateStateError('Aggregation pathway identity mismatch')
        return path

    def initialize_pathways(self):
        inputs = read_json(self.layout.path('request.json'))['inputs']
        for index, label in enumerate(self.data['pathways']):
            path = self.layout.path(f'pathways/pathway_{index:03d}/pathway.json')
            if path.exists():
                self.pathway(index, label)
            else:
                write_json(path, dict(pathway_schema_version=1, pathway_id=str(uuid4()),
                    run_id=self.data['run_id'], index=index, label=label,
                    addition_sequence=[inputs[string.ascii_lowercase.index(char)]['artifact_id'] for char in label],
                    composition=[dict(component_artifact_id=ref['artifact_id'], count=count)
                                 for ref, count in zip(inputs, self.data['request']['aggregate_sizes'])],
                    scientific_request=self.data['request']), overwrite=False)
            self.data['pathway_records'][str(index)] = sha256_file(path)
        self.save()

    def selection_path(self, index, step):
        return self.layout.path(f'pathways/pathway_{index:03d}/steps/step_{step:03d}/selection.json')

    def selection(self, index, step):
        record = read_json(self.selection_path(index, step))
        validate_uuid(record['record_id'])
        if (record['run_id'] != self.data['run_id'] or record['scientific_request'] != self.data['request']
                or type(record['stage_result_schema_version']) is not int or record['stage_result_schema_version'] != 1
                or type(record['pathway']) is not int or record['pathway'] != index
                or type(record['step']) is not int or record['step'] != step or record['status'] != 'completed'
                or record['physical_method'] != self.data['request']['backend_parameters']
                or record['scientific_settings'] != self.data['request']
                or len(record['selected']) > self.data['request']['maximum_number_of_seeds']):
            raise AggregateStateError('Invalid aggregation selection')
        return record

    def verify_selected(self, record):
        verify_refs(record['selected'], self.aggregate_directory)
        for item in record['native_sources']:
            ArtifactRef.from_dict(item['artifact']).verify(self.aggregate_directory)

    def commit_step(self, index, step, molecules, native, logical_home):
        directory = self.layout.path(f'pathways/pathway_{index:03d}/steps/step_{step:03d}/selected')
        refs = durable_pool(molecules, directory, self.aggregate_directory)
        sources = []
        if native is not None:
            for source_index, path in enumerate(sorted(native.rglob('selected/result_*.xyz'))):
                relative = path.relative_to(self.aggregate_directory).as_posix()
                self.layout.path(relative)
                snapshot = directory / f'native_source_{source_index:03d}.xyz'
                copy_result(path, snapshot)
                sources.append(dict(legacy_path=str(Path(logical_home) / path.relative_to(native)),
                    artifact=ArtifactRef.capture(self.aggregate_directory, snapshot.relative_to(self.aggregate_directory).as_posix(),
                                                'cross-path-selection-input').to_dict()))
        previous = self.data['steps'].get(str(index), {}).get(str(step - 1)) if step else None
        write_json(self.selection_path(index, step), dict(stage_result_schema_version=1,
            record_id=str(uuid4()), run_id=self.data['run_id'], pathway=index, step=step,
            status='completed', scientific_request=self.data['request'], selected=refs,
            physical_method=self.data['request']['backend_parameters'], scientific_settings=self.data['request'],
            native_sources=sources, upstream_sha256=previous), overwrite=False)
        self.adopt_step(index, step)

    def adopt_step(self, index, step):
        record = self.selection(index, step)
        previous = self.data['steps'].get(str(index), {}).get(str(step - 1)) if step else None
        if record['upstream_sha256'] != previous:
            raise AggregateStateError('Stale aggregation step dependency')
        self.verify_selected(record)
        self.data['steps'].setdefault(str(index), {})[str(step)] = sha256_file(self.selection_path(index, step))
        self.save()
        return verify_refs(record['selected'], self.aggregate_directory)

    def complete_pathway(self, index, label):
        self.data['completed_pathways'].append(dict(index=index, label=label, status='completed'))
        self.save()

    def native_sources(self):
        sources = {}
        for index, steps in self.data['steps'].items():
            for step in steps:
                for item in self.selection(int(index), int(step))['native_sources']:
                    path = ArtifactRef.from_dict(item['artifact']).verify(self.aggregate_directory)
                    if item['legacy_path'] in sources.values():
                        raise AggregateStateError('Duplicate logical aggregation source')
                    sources[path] = item['legacy_path']
        return sources

    def finish(self, selected):
        record_id = str(uuid4())
        target_composition = Counter()
        for fragment, count in zip(self.data['request']['fragments'], self.data['request']['aggregate_sizes']):
            target_composition.update({atom: fragment['atoms'].count(atom) * count for atom in set(fragment['atoms'])})
        source_ids = {}
        for pathway, steps in self.data['steps'].items():
            for step in steps:
                for item in self.selection(int(pathway), int(step))['native_sources']:
                    key = str(Path(item['legacy_path']).with_suffix('')).replace(os.sep, '_')
                    source_ids[key] = item['artifact']['artifact_id']
        refs = []
        for index, molecule in enumerate(selected):
            source = self.layout.path('selection/attempts/attempt_000/execution/selected') / molecule.relative_path
            target = self.layout.path(f'results/aggregate_{index:03d}.xyz')
            copy_result(source, target)
            ref = ArtifactRef.capture(self.aggregate_directory, target.relative_to(self.aggregate_directory).as_posix(),
                                      'final-cross-path-selection').to_dict()
            ref.update(rank=index, selection_name=molecule.name, source_record_id=record_id,
                       source_artifact_id=source_ids[molecule.name],
                       matches_requested_composition=Counter(molecule.atoms_list) == target_composition,
                       composition={atom: molecule.atoms_list.count(atom) for atom in sorted(set(molecule.atoms_list))},
                       energy_record=None if molecule.energy is None else
                       dict(value=molecule.energy, unit='hartree', convention='physical_potential'))
            refs.append(ref)
        write_json(self.layout.path('final_selection.json'), dict(stage_result_schema_version=1,
            record_id=record_id, run_id=self.data['run_id'], scientific_request=self.data['request'],
            physical_method=self.data['request']['backend_parameters'], scientific_settings=self.data['request'],
            requested_composition=dict(target_composition), dependencies=self.data['steps'], selected=refs, status='completed'), overwrite=False)
        self.adopt_final()

    def final_record(self):
        record = read_json(self.layout.path('final_selection.json'))
        validate_uuid(record['record_id'])
        if (record['run_id'] != self.data['run_id'] or record['scientific_request'] != self.data['request']
                or record['dependencies'] != self.data['steps'] or record['status'] != 'completed'):
            raise AggregateStateError('Final aggregation selection provenance mismatch')
        for index, ref in enumerate(record['selected']):
            if (ref['source_record_id'] != record['record_id'] or ref['rank'] != index
                    or ref['path'] != f'results/aggregate_{index:03d}.xyz'):
                raise AggregateStateError('Final aggregation result identity mismatch')
            ArtifactRef.from_dict({key: ref[key] for key in ('artifact_id','role','path','sha256')}).verify(self.aggregate_directory)
        return record

    def adopt_final(self):
        record = self.final_record()
        self.data.update(status='completed', final_selected_results=[ref['path'] for ref in record['selected']],
                         final_selection_sha256=sha256_file(self.layout.path('final_selection.json')))
        self.save()
        write_json(self.layout.path('results/manifest.json'), dict(format='pyar-results-manifest',
            manifest_schema_version=1, workflow='aggregate', run_id=self.data['run_id'], artifacts=record['selected']))
