"""Audit fixture evidence, then run isolated, bidirectional comparator checks.

Run from the checkout with its environment:
python benchmarks/structure_comparison/systematic_test.py --output /tmp/pyar-comparison
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
METHODS = ('coulomb_eigenvalue_energy_rule', 'coulomb_eigenvalue_prefilter_rmsd',
           'graph_rmsd', 'irmsd', 'graph_irmsd_policy')
THRESHOLDS = (0.05, 0.1, 0.2, 0.5)


def load_dataset(root=ROOT):
    lines = (root / 'structures.xyz').read_text().splitlines()
    frames = {}
    i = 0
    while i < len(lines):
        n = int(lines[i])
        if n <= 0 or i + n + 2 > len(lines):
            raise ValueError(f'Invalid frame at line {i + 1}')
        meta = dict(item.split('=', 1) for item in lines[i + 1].split())
        key = meta['structure_id']
        if key in frames:
            raise ValueError(f'Duplicate structure ID: {key}')
        rows = [line.split() for line in lines[i + 2:i + n + 2]]
        if any(len(row) != 4 for row in rows):
            raise ValueError(f'Invalid XYZ rows: {key}')
        coords = np.array([[float(v) for v in row[1:]] for row in rows])
        if not np.isfinite(coords).all():
            raise ValueError(f'Nonfinite geometry: {key}')
        frames[key] = SimpleNamespace(name=key, atoms_list=[row[0] for row in rows],
                                      coordinates=coords, metadata=meta)
        i += n + 2
    with (root / 'pairs.csv').open(newline='') as stream:
        pairs = list(csv.DictReader(stream))
    seen, seen_edges = set(), set()
    for pair in pairs:
        key = pair['pair_id']
        edge = tuple(sorted((pair['structure_a'], pair['structure_b'])))
        if key in seen or edge in seen_edges or edge[0] == edge[1]:
            raise ValueError(f'Duplicate or self pair: {pair}')
        seen.add(key)
        seen_edges.add(edge)
        for sid in edge:
            if sid not in frames:
                raise ValueError(f'Unresolved frame {sid}')
    return frames, pairs


def audit_dataset(frames, pairs):
    # These checks do not call the comparators under evaluation.
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    errors, checks, geometry = [], [], []
    identities = {}
    for key, frame in frames.items():
        xyz = frame.coordinates
        distances = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=2)
        upper = distances[np.triu_indices(len(xyz), 1)]
        minimum = float(upper.min())
        if minimum < 1e-8:
            errors.append(f'Coincident atoms: {key}')
        row = {'structure_id': key, 'minimum_separation_angstrom': minimum,
               'severe_close_contact_below_0.5_A': minimum < 0.5,
               'category': frame.metadata.get('system_class', 'organic')}
        smiles = frame.metadata.get('smiles')
        if smiles:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            if [atom.GetSymbol() for atom in mol.GetAtoms()] != frame.atoms_list:
                errors.append(f'SMILES/XYZ atom order mismatch: {key}')
            identity = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=False)
            identities[key] = (rdMolDescriptors.CalcMolFormula(mol), identity)
            conf = Chem.Conformer(len(xyz))
            for atom, point in enumerate(xyz):
                conf.SetAtomPosition(atom, point.tolist())
            mol.AddConformer(conf)
            props = AllChem.MMFFGetMoleculeProperties(mol)
            ff = AllChem.MMFFGetMoleculeForceField(mol, props)
            residual = abs(ff.CalcEnergy() - float(frame.metadata['mmff_energy_kcal_mol']))
            row['energy_recomputed_error_kcal_mol'] = residual
            row['max_gradient_component_kcal_mol_A'] = float(np.abs(ff.CalcGrad()).max())
            if residual > 1e-4 or frame.metadata.get('mmff_converged') != 'true':
                errors.append(f'MMFF provenance invalid: {key}')
        # Fixed water O-H geometry and intermolecular contacts independently
        # establish that these are disconnected molecular-cluster fixtures.
        if key.startswith('water_'):
            for offset in range(0, len(xyz), 3):
                if not np.allclose(np.linalg.norm(xyz[offset + 1:offset + 3] - xyz[offset], axis=1),
                                   [0.9572, np.hypot(0.239, 0.927)], atol=2e-8):
                    errors.append(f'Water monomer changed: {key}')
            cross = [distances[a, b] for a in range(len(xyz)) for b in range(a + 1, len(xyz))
                     if a // 3 != b // 3]
            row['minimum_interfragment_separation_angstrom'] = float(min(cross))
            if min(cross) < 1.2:
                errors.append(f'Water intermolecular clash: {key}')
        geometry.append(row)
    for pair in pairs:
        a, b = (frames[pair[k]] for k in ('structure_a', 'structure_b'))
        label = pair['label']
        check = {'pair_id': pair['pair_id'], 'label': label, 'binary_expected': None}
        if label == 'same_structure':
            order = [int(i) for i in b.metadata['source_atom_order'].split(',')]
            if sorted(order) != list(range(len(a.atoms_list))):
                errors.append(f'Invalid permutation: {pair["pair_id"]}')
                continue
            ref = a.coordinates[order]
            if [a.atoms_list[i] for i in order] != b.atoms_list:
                errors.append(f'Element-changing permutation: {pair["pair_id"]}')
            # Fit only the recorded atom correspondence, with a proper rotation.
            x, y = ref - ref.mean(axis=0), b.coordinates - b.coordinates.mean(axis=0)
            u, _, vt = np.linalg.svd(x.T @ y)
            correction = np.eye(3)
            correction[-1, -1] = np.linalg.det(u @ vt)
            residual = float(np.max(np.linalg.norm(x @ u @ correction @ vt - y, axis=1)))
            check['recorded_mapping_max_residual_angstrom'] = residual
            check['binary_expected'] = True
            if residual > 1e-7:
                errors.append(f'Rigid transform not verified: {pair["pair_id"]}')
        elif label == 'different_connectivity':
            formula_a, identity_a = identities[a.name]
            formula_b, identity_b = identities[b.name]
            if formula_a != formula_b or identity_a == identity_b:
                errors.append(f'Invalid constitutional isomer label: {pair["pair_id"]}')
            check['binary_expected'] = False
            check['evidence'] = 'independent_SMILES_formula_and_connectivity'
        elif label == 'different_composition':
            if Counter(a.atoms_list) == Counter(b.atoms_list):
                errors.append(f'Composition not different: {pair["pair_id"]}')
            check['binary_expected'] = False
        elif label in {'different_cluster_motif', 'different_cluster_arrangement'}:
            # A mismatch proves noncongruence, but does not specify a basin or
            # impose a binary decision at an arbitrary RMSD threshold.
            def spectrum(coords):
                d = np.linalg.norm(coords[:, None] - coords[None, :], axis=2)
                return np.sort(d[np.triu_indices(len(coords), 1)])
            difference = float(np.max(np.abs(spectrum(a.coordinates) - spectrum(b.coordinates))))
            check['distance_spectrum_max_difference_angstrom'] = difference
            if difference < 1e-6:
                errors.append(f'Unproven different motif: {pair["pair_id"]}')
        elif label not in {'unreviewed_conformer_pair', 'coordinate_perturbation', 'ambiguous', 'distorted_cluster'}:
            errors.append(f'Unknown label: {label}')
        checks.append(check)
    return {'valid': not errors, 'errors': errors, 'structures': len(frames), 'pairs': len(pairs),
            'label_counts': dict(Counter(p['label'] for p in pairs)),
            'pair_evidence': checks, 'geometry_checks': geometry,
            'interpretation': 'Only constructed invariance, composition, and SMILES isomer labels are scored. Others are diagnostic.'}


def worker(method, a, b):
    sys.path.insert(0, str(REPO))
    from pyar.structure_comparison import (
        CoulombEigenvalueRMSDComparator, GraphFirstDeduplicationComparator,
        GraphRMSDComparator, IRMSDComparator,
    )
    frames, _ = load_dataset()
    before = [frames[s].coordinates.copy() for s in (a, b)]
    started = time.perf_counter()
    if method == 'coulomb_eigenvalue_energy_rule':
        from pyar.selection.deduplication import calc_fingerprint_distance

        first, second = frames[a], frames[b]

        def known_mmff_energy(frame):
            current = frame
            visited = set()
            while current.name not in visited:
                visited.add(current.name)
                if 'mmff_energy_kcal_mol' in current.metadata:
                    return float(current.metadata['mmff_energy_kcal_mol'])
                source = current.metadata.get('source_structure')
                if not source or source not in frames:
                    return None
                current = frames[source]
            return None

        same_atom_count = len(first.atoms_list) == len(second.atoms_list)
        first_energy, second_energy = known_mmff_energy(first), known_mmff_energy(second)
        if not same_atom_count:
            result = SimpleNamespace(compatible=False, distance=None, equivalent=False,
                                     metadata={'rule_status': 'incompatible_atom_count'})
        elif first_energy is None or second_energy is None:
            result = SimpleNamespace(compatible=True, distance=None, equivalent=None,
                                     metadata={'rule_status': 'energy_unavailable'})
        else:
            fingerprint_distance = calc_fingerprint_distance(first, second)
            equivalent = (Counter(first.atoms_list) == Counter(second.atoms_list)
                          and abs(first_energy - second_energy) < 1e-5
                          and abs(fingerprint_distance) < 1.0)
            result = SimpleNamespace(
                compatible=True, distance=None, equivalent=equivalent,
                metadata={'rule_status': 'energy_and_sorted_coulomb_eigenvalue_distance',
                          'energy_difference': abs(first_energy - second_energy),
                          'fingerprint_distance': float(fingerprint_distance)},
            )
    else:
        comparator = {'graph_rmsd': GraphRMSDComparator, 'irmsd': IRMSDComparator,
                      'graph_irmsd_policy': GraphFirstDeduplicationComparator,
                      'coulomb_eigenvalue_prefilter_rmsd': CoulombEigenvalueRMSDComparator}[method](threshold=0.1)
        result = comparator.compare(frames[a], frames[b])
    elapsed = time.perf_counter() - started
    payload = dict(compatible=bool(result.compatible), distance=result.distance,
                   equivalent=None if result.equivalent is None else bool(result.equivalent),
                   metadata=dict(result.metadata), seconds=elapsed,
                   inputs_unchanged=all(np.array_equal(old, frames[s].coordinates)
                                        for old, s in zip(before, (a, b))))
    print('RESULT_JSON:' + json.dumps(payload, allow_nan=False), flush=True)


def run_case(case):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    try:
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker',
                               case['method'], '--a', case['a'], '--b', case['b']],
                              capture_output=True, text=True, timeout=30, env=env)
        lines = proc.stdout.splitlines()
        payloads = [line.removeprefix('RESULT_JSON:') for line in lines if line.startswith('RESULT_JSON:')]
        diagnostics = '\n'.join(line for line in lines if not line.startswith('RESULT_JSON:')) + proc.stderr
        row = dict(case, diagnostics=diagnostics.strip(), returncode=proc.returncode)
        if proc.returncode or len(payloads) != 1:
            return dict(row, status='error')
        row.update(json.loads(payloads[0]))
        if 'falling back' in diagnostics.lower():
            status = 'backend_fallback'
        elif diagnostics.strip():
            status = 'backend_diagnostic'
        elif not row['inputs_unchanged']:
            status = 'input_mutation'
        elif (row['metadata'].get('comparison_complete') is False
              and row['metadata'].get('fallback_status') != 'ok'):
            status = 'incomplete'
        elif row['distance'] is not None and (not np.isfinite(row['distance']) or row['distance'] < 0):
            status = 'invalid_distance'
        else:
            status = 'ok'
        return dict(row, status=status)
    except subprocess.TimeoutExpired:
        return dict(case, status='timeout')


def summarize(rows, audit):
    expected = {p['pair_id']: p['binary_expected'] for p in audit['pair_evidence']}
    summary = {}
    for method in METHODS:
        selected = [r for r in rows if r['method'] == method]
        forward = [r for r in selected if r['direction'] == 'forward']
        self_rows = [r for r in selected if r['direction'] == 'self']
        sweeps = {}
        thresholds = ('energy_eigenvalue_rule',) if method == 'coulomb_eigenvalue_energy_rule' else THRESHOLDS
        for threshold in thresholds:
            counts = Counter()
            for r in forward:
                truth = expected[r['case_id']]
                if truth is None:
                    counts['unscored'] += 1
                    continue
                if r['status'] != 'ok':
                    counts['abstained_positive' if truth else 'abstained_negative'] += 1
                    continue
                if method == 'coulomb_eigenvalue_energy_rule':
                    decision = r['equivalent'] is True
                else:
                    decision = r['compatible'] and (r['distance'] is not None and r['distance'] < threshold)
                if r['distance'] is None and r['compatible'] and r['equivalent'] is None:
                    counts['abstained_positive' if truth else 'abstained_negative'] += 1
                else:
                    counts['true_positive' if truth and decision else 'false_split' if truth
                           else 'false_merge' if decision else 'true_negative'] += 1
            sweeps[str(threshold)] = dict(counts)
        by_key = {(r['case_id'], r['direction']): r for r in selected}
        symmetry = []
        for r in forward:
            reverse = by_key[(r['case_id'], 'reverse')]
            if r['status'] == reverse['status'] == 'ok':
                if method == 'coulomb_eigenvalue_energy_rule' and r['equivalent'] != reverse['equivalent']:
                    symmetry.append({'pair_id': r['case_id'], 'kind': 'decision'})
                elif r['compatible'] != reverse['compatible'] or (r['distance'] is None) != (reverse['distance'] is None):
                    symmetry.append({'pair_id': r['case_id'], 'kind': 'decision_or_availability'})
                elif r['distance'] is not None and abs(r['distance'] - reverse['distance']) > 1e-6:
                    symmetry.append({'pair_id': r['case_id'], 'kind': 'distance',
                                     'difference_angstrom': abs(r['distance'] - reverse['distance'])})
        summary[method] = {'forward_status': dict(Counter(r['status'] for r in forward)),
                           'all_call_status': dict(Counter(r['status'] for r in selected)),
                           'threshold_sweep': sweeps, 'directional_disagreements': symmetry,
                           'self_status': dict(Counter(r['status'] for r in self_rows)),
                           'nonzero_or_rejected_self': [r['case_id'] for r in self_rows if r['status'] == 'ok'
                               and ((method == 'coulomb_eigenvalue_energy_rule' and r['equivalent'] is False)
                                    or (method != 'coulomb_eigenvalue_energy_rule'
                                        and (not r['compatible'] or r['distance'] is None
                                             or r['distance'] > 1e-6)))],
                           'forward_comparison_seconds': sum(r.get('seconds', 0) for r in forward)}
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--worker', choices=METHODS)
    parser.add_argument('--a')
    parser.add_argument('--b')
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.a, args.b)
        return
    frames, pairs = load_dataset()
    audit = audit_dataset(frames, pairs)
    args.output.mkdir(parents=True, exist_ok=True)
    versions = {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'networkx', 'rdkit', 'irmsd')}
    audit['versions'] = dict(versions, python=sys.version)
    audit['sha256'] = {str(path.relative_to(REPO)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in [ROOT / 'structures.xyz', ROOT / 'pairs.csv', ROOT / 'generate_dataset.py',
                                    Path(__file__).resolve(), REPO / 'pyar/selection/deduplication.py',
                                    *sorted((REPO / 'pyar/structure_comparison').glob('*.py'))]}
    (args.output / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(f'Audit: {audit["structures"]} frames, {audit["pairs"]} pairs, {len(audit["errors"])} errors', flush=True)
    if not audit['valid']:
        raise RuntimeError(audit['errors'])
    if args.audit_only:
        return
    cases = []
    for method in METHODS:
        for p in pairs:
            for direction, a, b in [('forward', p['structure_a'], p['structure_b']),
                                    ('reverse', p['structure_b'], p['structure_a'])]:
                cases.append(dict(method=method, case_id=p['pair_id'], label=p['label'], direction=direction, a=a, b=b))
        for sid in frames:
            cases.append(dict(method=method, case_id=sid, label='self', direction='self', a=sid, b=sid))
    rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_case, case) for case in cases]
        for future in as_completed(futures):
            rows.append(future.result())
            if len(rows) % 100 == 0:
                print(f'Compared {len(rows)}/{len(cases)}', flush=True)
    rows.sort(key=lambda r: (r['method'], r['case_id'], r['direction']))
    (args.output / 'pairwise_results.jsonl').write_text(''.join(json.dumps(r, allow_nan=False) + '\n' for r in rows))
    summary = summarize(rows, audit)
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
