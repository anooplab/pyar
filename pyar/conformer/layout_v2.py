"""Conformer completion records; inspection metadata does not provide resume."""
from pathlib import Path
from uuid import uuid4

from pyar.persistence import ArtifactRef, RunLayout, read_json, write_json, detect_version
from pyar.persistence.json_store import fsync_directory
from pyar.persistence.artifacts import validate_uuid
from pyar.growth.layout_v2 import copy_result


class ConformerLayout:
    def __init__(self, root, request):
        self.root = Path(root).resolve()
        self.layout = RunLayout(self.root, 'conformer')
        self.run_id = str(uuid4())
        from pyar.workflows.conformer import _json_value
        request = _json_value(request)
        self.request = request
        self.input_refs = []
        from pyar.workflows.conformer import _resolve_input_format
        source = request['input']
        if _resolve_input_format(source, request['input_format']) != 'smiles':
            target = self.layout.path('inputs/source' + Path(source).suffix.lower())
            copy_result(Path(source), target)
        else:
            target = self.layout.path('inputs/input.smiles')
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(str(source) + '\n')
            with target.open('rb') as stream:
                import os
                os.fsync(stream.fileno())
        self.input_refs.append(self.capture(target.relative_to(self.root).as_posix(), 'input').to_dict())
        write_json(self.layout.path('inputs/metadata.json'), dict(inputs=self.input_refs, request=request), overwrite=False)
        write_json(self.root / 'state.json', dict(workflow_state_schema_version=4, workflow='conformer',
            run_id=self.run_id, status='running', request=request))
        write_json(self.layout.marker, dict(format='pyar-run-layout', layout_version=2,
            workflow='conformer', run_id=self.run_id), overwrite=False)

    def capture(self, relative, role):
        path = self.layout.path(relative)
        with path.open('rb') as stream:
            import os
            os.fsync(stream.fileno())
        parent = path.parent
        while parent.is_relative_to(self.root):
            fsync_directory(parent)
            if parent == self.root:
                break
            parent = parent.parent
        return ArtifactRef.capture(self.root, relative, role)

    def relative(self, value):
        if isinstance(value, dict):
            return {key: self.relative(item) for key, item in value.items()}
        if isinstance(value, list):
            return [self.relative(item) for item in value]
        if isinstance(value, str) and value.startswith(str(self.root) + '/'):
            return Path(value).relative_to(self.root).as_posix()
        return value

    def publish(self, state):
        from pyar.workflows.conformer import _json_value
        state = self.relative(_json_value(state))
        state.pop('version')
        state.update(workflow_state_schema_version=4, run_id=self.run_id)
        artifacts = list(self.input_refs)
        known = {ref['path'] for ref in artifacts}
        for path in sorted((self.root / 'inputs').iterdir()):
            relative = path.relative_to(self.root).as_posix()
            if path.is_file() and relative not in known and path.name != 'metadata.json':
                artifacts.append(self.capture(relative, 'converted-input').to_dict())
        generated = []
        selected = []
        ranks = {path: rank for rank, path in enumerate(state['selected_results'])}
        for record in state['records']:
            record = dict(record, conformer_id=str(uuid4()))
            record['rdkit_energy_record'] = dict(value=record['rdkit_energy'], unit='kcal/mol',
                                                  method=record['force_field'])
            record['backend_energy_record'] = (None if record['backend_energy'] is None else
                dict(value=record['backend_energy'], unit='hartree', convention='physical_potential'))
            source_id = None
            for field, role in (('rdkit_path', 'generated-conformer'), ('backend_path', 'backend-refined-conformer')):
                if record[field]:
                    ref = self.capture(record[field], role).to_dict()
                    artifacts.append(ref)
                    record[field.replace('_path', '_artifact')] = ref
                    source_id = ref['artifact_id']
            if record['selected_path']:
                ref = self.capture(record['selected_path'], 'selected-conformer').to_dict()
                ref.update(rank=ranks[record['selected_path']], conformer_id=record['conformer_id'],
                           source_artifact_id=source_id, seed=record['seed'], source_conf_id=record['source_conf_id'],
                           selection_reason=record['refinement_reason'],
                           rdkit_energy_record=record['rdkit_energy_record'],
                           backend_energy_record=record['backend_energy_record'])
                artifacts.append(ref)
                selected.append(ref)
            generated.append(record)
        selected.sort(key=lambda ref: ref['rank'])
        record_id = str(uuid4())
        selected = [dict(ref, source_record_id=record_id) for ref in selected]
        write_json(self.layout.path('generated/records.json'), dict(record_schema_version=1,
            run_id=self.run_id, conformers=generated), overwrite=False)
        artifacts.append(self.capture('generated/records.json', 'generation-metadata').to_dict())
        write_json(self.layout.path('completion.json'), dict(stage_result_schema_version=1, record_id=record_id,
            run_id=self.run_id, status=state['status'], request=self.request,
            comparison_policy=state['comparison_policy'], artifacts=artifacts, selected=selected), overwrite=False)
        write_json(self.layout.path('results/manifest.json'), dict(format='pyar-results-manifest',
            manifest_schema_version=1, workflow='conformer', run_id=self.run_id, artifacts=selected))
        write_json(self.root / 'state.json', state)
        write_json(self.root / 'summary.json', dict(status=state['status'], generated=len(generated),
                   selected=len(selected), run_id=self.run_id))


def selected_artifacts(root):
    """Read checked v2 selections without inferring results from filenames."""
    root = Path(root).resolve()
    detect_version(root, 'conformer')
    completion = read_json(root / 'completion.json')
    manifest = read_json(root / 'results/manifest.json')
    marker = read_json(root / 'layout.json')
    validate_uuid(completion['record_id'])
    if (completion.get('stage_result_schema_version') != 1
            or manifest.get('manifest_schema_version') != 1
            or manifest.get('format') != 'pyar-results-manifest'
            or manifest.get('workflow') != 'conformer'
            or completion.get('status') not in {'completed', 'completed_no_backend_success', 'failed_no_conformers'}
            or completion['run_id'] != marker['run_id'] or manifest['run_id'] != marker['run_id']
            or manifest['artifacts'] != completion['selected']):
        raise ValueError('Conformer completion and manifest identities differ')
    source_ids = {ref['artifact_id'] for ref in completion['artifacts']}
    for rank, ref in enumerate(manifest['artifacts']):
        validate_uuid(ref['conformer_id'])
        if (ref['rank'] != rank or ref['role'] != 'selected-conformer'
                or ref['source_record_id'] != completion['record_id']
                or ref['source_artifact_id'] not in source_ids
                or ref['path'] != f'results/conf_{rank:04d}.xyz'):
            raise ValueError('Invalid selected conformer provenance')
    for ref in completion['artifacts']:
        ArtifactRef.from_dict({key: ref[key] for key in ('artifact_id','role','path','sha256')}).verify(root)
    return manifest['artifacts']
