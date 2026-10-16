"""Reaction-owned discovery commits and checked immutable geometry snapshots."""
from copy import deepcopy
from pathlib import Path
from uuid import uuid4

from pyar import __version__
from pyar.persistence import (ArtifactRef, RunLayout, detect_version, read_json,
                              write_json, sha256_file, fsync_artifacts)
from pyar.growth.layout_v2 import copy_result
from pyar.persistence.artifacts import validate_uuid
from pyar.persistence.json_store import fsync_directory
from pyar.state.reaction import ReactionRunState, ReactionStateError, _json_value


class ReactionRunStateV2(ReactionRunState):
    def __init__(self, root_directory, data):
        super().__init__(root_directory, data)
        self.state_directory = self.reaction_directory / 'checkpoints'
        self.geometry_directory = self.state_directory / 'geometries'
        self.layout = RunLayout(self.reaction_directory, 'reaction')

    @classmethod
    def create(cls, root_directory, request, orientations, reactants, sampling=None):
        data = dict(workflow_state_schema_version=3, workflow='reaction', run_id=str(uuid4()),
            package_version=__version__, status='running', request=_json_value(request),
            gamma_schedule=[float(value) for value in request['gamma_schedule']], current_cycle=0,
            completed_cycles=[], completed_jobs=[], products=[], pending_orientations=[],
            current_survivors=[], input_reactants=[], next_snapshot=0, cycle_records={})
        if sampling is not None:
            data['sampling'] = _json_value(sampling)
        state = cls(root_directory, data)
        state.geometry_directory.mkdir(parents=True, exist_ok=True)
        for molecule, label in zip(reactants, ('a', 'b')):
            ref = state._snapshot_molecule(molecule, 'input')
            target = state.layout.path(f'inputs/reactant_{label}.xyz')
            copy_result(state.reaction_directory / ref['path'], target)
            (state.reaction_directory / ref['path']).unlink()
            ref['path'] = target.relative_to(state.reaction_directory).as_posix()
            state.data['input_reactants'].append(ref)
        state._replace_pending_orientations(orientations)
        write_json(state.layout.path('inputs/metadata.json'), dict(inputs=state.data['input_reactants']), overwrite=False)
        write_json(state.layout.path('request.json'), dict(request_schema_version=1, workflow='reaction',
            run_id=data['run_id'], scientific_request=data['request']), overwrite=False)
        state.save()
        write_json(state.layout.marker, dict(format='pyar-run-layout', layout_version=2,
            workflow='reaction', run_id=data['run_id']), overwrite=False)
        return state

    @classmethod
    def load(cls, root_directory, expected_request, *, completed=False):
        directory = Path(root_directory).resolve() / 'reaction'
        try:
            detect_version(directory, 'reaction')
            state = cls(root_directory, read_json(directory / 'state.json'))
            if completed:
                record = state.discovery_record()
                state = cls(root_directory, record['completed_state'])
            else:
                state.validate_progress()
                # Product publication can precede the per-job progress save.
                for path in sorted(state.layout.path('discovery/products').glob('product_*/product.json')):
                    product_record = read_json(path)
                    if product_record['run_id'] != state.data['run_id']:
                        raise ReactionStateError('Accepted product identity mismatch')
                    product = product_record['product']
                    validate_uuid(product['product_id'])
                    state.verify_ref(product['artifact'])
                    for ref in product.get('evidence', []):
                        state.verify_ref(ref)
                    if not any(item['product_id'] == product['product_id'] for item in state.data['products']):
                        state.data['products'].append(product)
                # A committed cycle may precede its mutable progress publication.
                while state.cycle_path(state.data['current_cycle']).exists():
                    index = state.data['current_cycle']
                    record = read_json(state.cycle_path(index))
                    state.verify_cycle(record, index)
                    state = cls(root_directory, record['checkpoint'])
                    state.data['cycle_records'][str(index)] = sha256_file(state.cycle_path(index))
                if state.layout.path('discovery/commit.json').exists():
                    raise ReactionStateError('Reaction discovery is completed; use characterization or a fresh directory')
                if state.data['status'] != 'running':
                    raise ReactionStateError(f"Reaction state is already {state.data['status']!r}")
            state.validate_progress()
            state.validate_request(expected_request)
            identity = read_json(state.layout.path('request.json'))
            if identity['scientific_request'] != state.data['request']:
                raise ReactionStateError('Reaction request differs from immutable input identity')
            for role in ('pending_orientations', 'current_survivors', 'input_reactants'):
                for ref in state.data[role]:
                    state.verify_ref(ref)
            for index, digest in state.data['cycle_records'].items():
                if sha256_file(state.cycle_path(int(index))) != digest:
                    raise ReactionStateError('Reaction cycle commit was modified')
                state.verify_cycle(read_json(state.cycle_path(int(index))), int(index))
            return state
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise ReactionStateError(f'Invalid v2 reaction state: {exc}') from exc

    @classmethod
    def read_completed(cls, root_directory, expected_request):
        return cls.load(root_directory, expected_request, completed=True)

    def validate_progress(self):
        index = self.data['current_cycle']
        schedule = self.data['gamma_schedule']
        if (type(index) is not int or not 0 <= index <= len(schedule)
                or schedule != self.data['request']['gamma_schedule']
                or self.data['workflow_state_schema_version'] != 3
                or len(self.data['completed_cycles']) != index
                or set(self.data['cycle_records']) != {str(i) for i in range(index)}
                or [cycle['gamma'] for cycle in self.data['completed_cycles']] != schedule[:index]):
            raise ReactionStateError('Invalid reaction cycle progress')

    def verify_ref(self, reference):
        return ArtifactRef.from_dict({key: reference[key] for key in
            ('artifact_id', 'role', 'path', 'sha256')}).verify(self.reaction_directory)

    def _snapshot_molecule(self, molecule, role):
        ref = super()._snapshot_molecule(molecule, role)
        artifact = ArtifactRef.capture(self.reaction_directory, ref['path'], role)
        with (self.reaction_directory / ref['path']).open('rb') as stream:
            import os
            os.fsync(stream.fileno())
        parent = self.geometry_directory
        while parent.is_relative_to(self.reaction_directory):
            fsync_directory(parent)
            if parent == self.reaction_directory:
                break
            parent = parent.parent
        ref.update(artifact.to_dict())
        ref['energy_record'] = (None if ref['energy'] is None else
                                dict(value=ref['energy'], unit='hartree', convention='physical_potential'))
        return ref

    def _load_molecule(self, reference):
        try:
            self.verify_ref(reference)
            return super()._load_molecule(reference)
        except (OSError, ValueError, KeyError) as exc:
            raise ReactionStateError(f'Invalid reaction snapshot: {exc}') from exc

    def save(self):
        write_json(self.state_file, _json_value(self.data))

    def cycle_path(self, index):
        return self.layout.path(f'discovery/cycles/cycle_{index:03d}/cycle.json')

    def cycle_execution(self, gamma_id, gamma):
        index = self.data['current_cycle']
        root = self.layout.path(f'discovery/cycles/cycle_{index:03d}')
        request = dict(run_id=self.data['run_id'], cycle=index, gamma=float(gamma),
                       physical_settings=self.data['request']['backend_parameters'])
        path = root / 'request.json'
        if path.exists():
            if read_json(path) != request:
                raise ReactionStateError('Reaction cycle settings mismatch')
        else:
            write_json(path, request, overwrite=False)
        return root / f'attempts/attempt_000/execution/gamma_{gamma_id}'

    def verify_cycle(self, record, index):
        validate_uuid(record['record_id'])
        if (record['stage_result_schema_version'] != 1
                or record['run_id'] != self.data['run_id'] or record['cycle'] != index
                or record['gamma'] != self.data['gamma_schedule'][index]
                or record['checkpoint']['current_cycle'] != index + 1
                or record['checkpoint']['request'] != self.data['request']):
            raise ReactionStateError('Reaction cycle identity mismatch')
        for ref in record['artifacts']:
            self.verify_ref(ref)

    def complete_cycle(self, gamma, next_orientations):
        index = self.data['current_cycle']
        self.data['completed_cycles'].append(dict(gamma=float(gamma), status='completed'))
        self.data['current_cycle'] += 1
        self._replace_pending_orientations(next_orientations)
        self.data['current_survivors'] = []
        cycle = self.cycle_path(index).parent
        evidence_paths = [path for path in sorted(cycle.rglob('*'))
                          if path.is_file() and path != self.cycle_path(index)]
        fsync_artifacts(evidence_paths, self.reaction_directory)
        refs = [ArtifactRef.capture(self.reaction_directory, path.relative_to(self.reaction_directory).as_posix(),
                'discovery-cycle-evidence').to_dict() for path in evidence_paths]
        refs.extend(self.data['pending_orientations'])
        write_json(self.cycle_path(index), dict(stage_result_schema_version=1, record_id=str(uuid4()),
            run_id=self.data['run_id'], cycle=index, gamma=float(gamma), artifacts=refs,
            checkpoint=deepcopy(self.data)), overwrite=False)
        self.data['cycle_records'][str(index)] = sha256_file(self.cycle_path(index))
        self.save()

    def _relative_values(self, value):
        if isinstance(value, dict):
            return {key: self._relative_values(item) for key, item in value.items()}
        if isinstance(value, list):
            return [self._relative_values(item) for item in value]
        if isinstance(value, str) and value.startswith(str(self.reaction_directory) + '/'):
            return Path(value).relative_to(self.reaction_directory).as_posix()
        return value

    def record_product(self, job_name, gamma, inchi, smiles, path, trace_summary=None):
        index = len(self.data['products']) + 1
        relative = f'discovery/products/product_{index:03d}/structure.xyz'
        copy_result(Path(path), self.layout.path(relative))
        ref = ArtifactRef.capture(self.reaction_directory, relative, 'accepted-unbiased-product').to_dict()
        product = dict(job_name=job_name, gamma=float(gamma), inchi=inchi, smiles=smiles, path=relative,
                       product_id=str(uuid4()), artifact=ref)
        if trace_summary:
            summary = _json_value(trace_summary)
            candidate = Path(summary.get('candidate_ts_directory', ''))
            if not candidate.is_absolute():
                candidate = self.reaction_directory / candidate
            job = candidate.parent
            destination = self.layout.path(f'discovery/products/product_{index:03d}/sources')
            evidence = []
            paths = [file for folder in (job / 'reaction_trace', candidate)
                     for file in folder.rglob('*') if file.is_file()]
            paths.extend(file for file in (job.parent / f'trial_{job_name}.xyz',
                         job.parent / 'trial_relax.xyz', job.parent / 'result_relax.xyz',
                         job / 'pyar_geometric_state.json') if file.is_file())
            paths.extend(file for pattern in ('release_*.json', 'release_*.jsonl', 'pre_release_*.xyz')
                         for file in job.parent.glob(pattern) if file.is_file())
            for source in paths:
                target = destination / source.relative_to(job.parent)
                copy_result(source, target)
                evidence.append(ArtifactRef.capture(self.reaction_directory,
                    target.relative_to(self.reaction_directory).as_posix(), 'accepted-product-evidence').to_dict())
            summary['candidate_ts_directory'] = str(destination / candidate.relative_to(job.parent))
            product['trace_summary'] = self._relative_values(summary)
            product['evidence'] = evidence
        record = dict(product_schema_version=1, record_id=str(uuid4()), run_id=self.data['run_id'], product=product)
        write_json(self.layout.path(f'discovery/products/product_{index:03d}/product.json'), record, overwrite=False)
        self.data['products'].append(product)
        self.save()

    def selected_paths(self):
        return tuple(f'results/product_{index:03d}.xyz' for index in range(1, len(self.data['products']) + 1))

    def discovery_record(self):
        record = read_json(self.layout.path('discovery/commit.json'))
        validate_uuid(record['record_id'])
        if (record['stage_result_schema_version'] != 1
                or record['completed_state']['status'] not in {'completed', 'completed_products_found',
                    'completed_no_products', 'completed_no_candidates'}
                or record['run_id'] != self.data['run_id'] or record['status'] != 'completed'
                or record['completed_state']['run_id'] != self.data['run_id']):
            raise ReactionStateError('Discovery commit identity mismatch')
        for ref in record['artifacts']:
            self.verify_ref(ref)
        return record

    def finish(self, status='completed'):
        path = self.layout.path('discovery/commit.json')
        if path.exists():
            self.data = self.discovery_record()['completed_state']
            self.save()
            return
        self.data.update(status=status, pending_orientations=[])
        refs = []
        discovery_files = [native for native in sorted(self.layout.path('discovery').rglob('*'))
                           if native.is_file() and native != path]
        fsync_artifacts(discovery_files, self.reaction_directory)
        for native in discovery_files:
            if native.is_file():
                refs.append(ArtifactRef.capture(self.reaction_directory,
                    native.relative_to(self.reaction_directory).as_posix(), 'discovery-evidence').to_dict())
        results = []
        record_id = str(uuid4())
        for index, product in enumerate(self.data['products'], 1):
            target = self.layout.path(f'results/product_{index:03d}.xyz')
            copy_result(self.reaction_directory / product['path'], target)
            results.append(dict(ArtifactRef.capture(self.reaction_directory,
                target.relative_to(self.reaction_directory).as_posix(), 'accepted-unbiased-product').to_dict(),
                source_artifact_id=product['artifact']['artifact_id'], product_id=product['product_id'],
                source_record_id=record_id))
        write_json(path, dict(stage_result_schema_version=1, record_id=record_id, run_id=self.data['run_id'],
            status='completed', physical_settings=self.data['request'], artifacts=refs,
            completed_state=deepcopy(self.data), products=self.data['products']), overwrite=False)
        write_json(self.layout.path('results/manifest.json'), dict(format='pyar-results-manifest',
            manifest_schema_version=1, workflow='reaction', run_id=self.data['run_id'], artifacts=results))
        self.save()
