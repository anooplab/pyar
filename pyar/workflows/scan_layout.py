"""Scan-owned immutable records; native scan and path executors retain filenames."""
import hashlib
import json
import math
from pathlib import Path
from uuid import uuid4

from pyar.persistence import ArtifactRef, RunLayout, read_json, write_json, sha256_file, fsync_artifacts
from pyar.persistence.artifacts import validate_uuid
from pyar.reaction_identity import same_molecular_identity


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
        fsync_artifacts(paths, self.root)
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
            fsync_artifacts(paths, self.root)
            refs = [ArtifactRef.capture(self.root, file.relative_to(self.root).as_posix(),
                                       'free-relaxation-evidence').to_dict() for file in paths]
            write_json(path, dict(stage_result_schema_version=1, record_id=str(uuid4()),
                run_id=self.run_id, orientation=index, status='completed', scan_commit_sha256=sha256_file(self.commit_path(index)),
                artifacts=refs, result=self._paths(result)), overwrite=False)
        return None

    def _product_sources(self, index, request):
        """Verify the immutable scan and free-relaxation evidence for an orientation."""
        scan_path = self.commit_path(index)
        scan_result = self.restore(index, request)
        if scan_result is None:
            raise ValueError(f"Eligible product orientation {index} has no scan commit")
        scan_record = read_json(scan_path)
        relax_path = self.layout.path(f'orientation_{index:03d}/relax_commit.json')
        relax_result = self.relaxed(index)
        if relax_result is None:
            raise ValueError(f"Eligible product orientation {index} has no relaxation commit")
        relax_record = read_json(relax_path)
        # Reuse the immutable commit payloads, not the mutable result.json files,
        # as the registry's scientific source.
        if scan_record['result'] != self._paths(scan_result):
            raise ValueError(f"Scan report does not match its commit for orientation {index}")
        if self._paths(relax_record['result'], restore=True) != relax_result:
            raise ValueError(f"Relaxation report does not match its commit for orientation {index}")
        return {
            "scan_result": scan_result,
            "relax_result": relax_result,
            "scan_commit": scan_record,
            "relax_commit": relax_record,
            "scan_commit_sha256": sha256_file(scan_path),
            "relax_commit_sha256": sha256_file(relax_path),
            "scan_artifacts": scan_record['artifacts'],
            "relaxation_artifacts": relax_record['artifacts'],
        }

    def _trace_selection(self, index, request, source):
        """Score a route only for representative selection, never as a barrier."""
        from pyar.backends.orca_scan import parse_xyz_trajectory

        symbols = request['orientation_definitions'][index]['molecule']['atoms']
        scan_root = self.root / f"orientation_{index:03d}/attempts/attempt_000/execution/scan"
        trajectory = scan_root / "scan_trajectory.xyz"
        profile_path = scan_root / "scan_profile.json"
        try:
            frames = parse_xyz_trajectory(trajectory, len(symbols), symbols)
            profile = json.loads(profile_path.read_text(encoding="utf-8"))['points']
            expected = int(source['relax_result']['scan_points'])
            energies = [float(item['energy_hartree']) for item in profile]
            trace_complete = (len(frames) == expected == len(profile)
                              and all(math.isfinite(value) for value in energies))
        except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
            profile = []
            energies = []
            trace_complete = False
        candidate_paths = {ref['path'] for ref in source['scan_artifacts']
                           if ref['role'] == 'scan-maximum-candidate'}
        metadata_path = self.root / f"orientation_{index:03d}/attempts/attempt_000/execution/ts_candidates/metadata.json"
        highest_path = self.root / f"orientation_{index:03d}/attempts/attempt_000/execution/ts_candidates/highest_scan_energy.xyz"
        internal_maximum = False
        if trace_complete and metadata_path.relative_to(self.root).as_posix() in candidate_paths \
                and highest_path.relative_to(self.root).as_posix() in candidate_paths:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            internal_maximum = bool(metadata.get('maximum_is_internal'))
        peak_delta = (max(energies) - energies[0]) * 627.509474 if trace_complete else math.inf
        return {
            "complete_validated_trace": trace_complete,
            "usable_internal_scan_maximum": internal_maximum,
            "scan_peak_delta_kcal_mol": peak_delta if math.isfinite(peak_delta) else None,
        }

    def _write_immutable_product(self, product_id, record):
        path = self.layout.path(f"products/{product_id}/product.json")
        if path.exists():
            if read_json(path) != record:
                raise ValueError(f"Immutable scan product record changed: {product_id}")
        else:
            write_json(path, record, overwrite=False)
        return path

    def product_registry(self, results, request):
        """Group successful relaxed scan products and rebuild the v2 registry.

        Only changed, successfully identified products are grouped. Unknown
        identities remain separate orientation records; unchanged reactants and
        failed orientations remain visible in their own index sections.
        """
        ordered = sorted(results, key=lambda row: int(row.get('orientation', -1)))
        groups = []
        unknown = []
        non_products = []
        ineligible = []

        for row in ordered:
            index = int(row.get('orientation', -1))
            if row.get('status') != 'scan_success_relax_success':
                ineligible.append({"orientation": index, "status": row.get('status')})
                continue
            source = self._product_sources(index, request)
            committed = source['relax_result']
            if committed != row:
                # Registry pointers are derived and may be present in a previous
                # report; compare scientific identity/status independently.
                for key in ('status', 'identity_status', 'product_identity_changed',
                            'relaxed_identity', 'target_bond_present_after_relaxation'):
                    if row.get(key) != committed.get(key):
                        raise ValueError(f"Orientation {index} result differs from its relaxation commit")
            identity = committed.get('relaxed_identity') or {}
            if committed.get('identity_status') != 'success' or not identity.get('inchi'):
                unknown.append((index, committed, source))
                continue
            if committed.get('product_identity_changed') is not True:
                non_products.append({"orientation": index, "reason": "identity_unchanged"})
                continue
            match = next((group for group in groups
                          if same_molecular_identity(group['identity'], identity)), None)
            if match is None:
                groups.append({"identity": identity, "members": [(index, committed, source)]})
            else:
                match['members'].append((index, committed, source))

        product_entries = []
        orientation_updates = {}
        seen_ids = {}
        for group in groups:
            identity = group['identity']
            inchi = identity['inchi']
            digest = hashlib.sha256(inchi.encode('utf-8')).hexdigest()
            product_id = f"product_{digest[:12]}"
            if product_id in seen_ids and seen_ids[product_id] != inchi:
                raise ValueError(f"InChI-derived product ID collision: {product_id}")
            seen_ids[product_id] = inchi
            members = sorted(group['members'], key=lambda item: item[0])
            scored = []
            for index, result, source in members:
                score = self._trace_selection(index, request, source)
                scored.append(((-int(score['complete_validated_trace']),
                                -int(score['usable_internal_scan_maximum']),
                                (score['scan_peak_delta_kcal_mol'] if score['scan_peak_delta_kcal_mol'] is not None else math.inf),
                                index), index, result, source, score))
            _sort_key, representative, _result, _source, representative_score = min(scored)
            orientation_ids = [item[0] for item in members]
            target_statuses = sorted({
                "target" if result.get('target_bond_present_after_relaxation') is True else
                "off_target" if result.get('target_bond_present_after_relaxation') is False else "unknown"
                for _index, result, _source in members
            })
            target_status = target_statuses[0] if len(target_statuses) == 1 else "mixed"
            smiles_variants = sorted({item[1].get('relaxed_identity', {}).get('smiles', '')
                                      for item in members
                                      if item[1].get('relaxed_identity', {}).get('smiles')})
            source_records = []
            for index, _result, source in members:
                source_records.append({
                    "orientation": index,
                    "scan_commit_record_id": source['scan_commit']['record_id'],
                    "scan_commit_sha256": source['scan_commit_sha256'],
                    "relax_commit_record_id": source['relax_commit']['record_id'],
                    "relax_commit_sha256": source['relax_commit_sha256'],
                    "scan_artifacts": source['scan_artifacts'],
                    "relaxation_artifacts": source['relaxation_artifacts'],
                })
            product_record = {
                "product_schema_version": 1,
                "run_id": self.run_id,
                "product_id": product_id,
                "identity": {"inchi": inchi,
                             "smiles": smiles_variants[0] if smiles_variants else None,
                             "smiles_variants": smiles_variants},
                "orientation_ids": orientation_ids,
                "representative_orientation": representative,
                "alternative_orientations": [value for value in orientation_ids if value != representative],
                "target_status": target_status,
                "target_status_by_orientation": {
                    str(index): ("target" if result.get('target_bond_present_after_relaxation') is True else
                                 "off_target" if result.get('target_bond_present_after_relaxation') is False else "unknown")
                    for index, result, _source in members
                },
                "representative_selection": representative_score,
                "sources": source_records,
            }
            product_path = self._write_immutable_product(product_id, product_record)
            product_entries.append({
                "product_id": product_id,
                "identity": product_record['identity'],
                "orientation_ids": orientation_ids,
                "representative_orientation": representative,
                "alternative_orientations": product_record['alternative_orientations'],
                "target_status": target_status,
                "characterization_status": "not_started",
                "characterization_path": None,
                "product_record_path": product_path.relative_to(self.root).as_posix(),
                "product_record_sha256": sha256_file(product_path),
            })
            for index, _result, _source in members:
                orientation_updates[index] = {
                    "product_id": product_id,
                    "product_group_size": len(members),
                    "product_role": "representative" if index == representative else "duplicate_route",
                    "representative_orientation": representative,
                    "characterization_reference": f"products/{product_id}/characterization",
                }

        unknown_entries = []
        for index, result, source in unknown:
            group_id = f"identity_unknown_orientation_{index:03d}"
            record = {
                "product_schema_version": 1,
                "run_id": self.run_id,
                "unknown_group_id": group_id,
                "product_id": None,
                "orientation": index,
                "identity_status": result.get('identity_status'),
                "identity_error": result.get('identity_error'),
                "scan_commit_record_id": source['scan_commit']['record_id'],
                "scan_commit_sha256": source['scan_commit_sha256'],
                "relax_commit_record_id": source['relax_commit']['record_id'],
                "relax_commit_sha256": source['relax_commit_sha256'],
                "scan_artifacts": source['scan_artifacts'],
                "relaxation_artifacts": source['relaxation_artifacts'],
            }
            path = self._write_immutable_product(group_id, record)
            unknown_entries.append({"unknown_group_id": group_id, "orientation": index,
                                    "product_record_path": path.relative_to(self.root).as_posix(),
                                    "product_record_sha256": sha256_file(path)})
            orientation_updates[index] = {
                "product_id": None, "product_group_id": group_id, "product_group_size": 1,
                "product_role": "identity_unknown", "representative_orientation": index,
                "characterization_reference": None,
            }

        for row in ordered:
            index = int(row.get('orientation', -1))
            update = orientation_updates.get(index)
            if update is not None:
                row.update(update)
            elif row.get('status') == 'scan_success_relax_success':
                row.update(product_id=None, product_group_id=None,
                           product_group_size=0, product_role="not_a_product",
                           representative_orientation=None, characterization_reference=None)
            write_json(self.layout.path(f"orientation_{index:03d}/result.json"), row)

        product_entries.sort(key=lambda item: item['product_id'])
        index_record = {
            "format": "pyar-scan-product-registry",
            "product_registry_schema_version": 1,
            "run_id": self.run_id,
            "request_signature": request['request_signature'],
            "products": product_entries,
            "identity_unknown": unknown_entries,
            "non_products": non_products,
            "ineligible_orientations": ineligible,
        }
        write_json(self.layout.path('products/index.json'), index_record)
        return index_record

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
