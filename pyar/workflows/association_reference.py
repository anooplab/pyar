"""Scientific references for bimolecular reaction endpoints."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np

from pyar.backends import write_xyz


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=str) + '\n')
    temporary.replace(path)


def prepare_separated_reactant_reference(reactants, qc_params, output, options):
    """Optimize and frequency-check each isolated reactant once.

    Returns a persistent E_A + E_B reference. Failed or non-minimum fragments
    remain explicit failures and never receive an energy-sum claim.
    """
    from pyar.backends.geometric import PyarGeometricCalculator
    from pyar.neb import _frequency, _optimize_geometry
    from pyar.reaction_identity import (molecule_identity_from_xyz, same_molecular_identity,
                                        separated_reactant_identity)

    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    request_rows = []
    sources = []
    for index, molecule in enumerate(reactants):
        source = output / f'fragment_{index + 1}_input.xyz'
        if source.exists():
            existing_hash = _sha256(source)
            write_xyz(molecule.atoms_list, molecule.coordinates, output / f'.fragment_{index + 1}_check.xyz', precision=12)
            check = output / f'.fragment_{index + 1}_check.xyz'
            expected_hash = _sha256(check)
            check.unlink()
            if existing_hash != expected_hash:
                raise ValueError('Separated-reactant input snapshot changed; use a new characterization directory')
        else:
            write_xyz(molecule.atoms_list, molecule.coordinates, source, precision=12)
        sources.append(source)
        request_rows.append({
            'atoms': list(molecule.atoms_list), 'source_sha256': _sha256(source),
            'charge': int(molecule.charge), 'multiplicity': int(molecule.multiplicity),
            'scftype': getattr(molecule, 'scftype', 'rhf'),
        })
    physical = {key: qc_params.get(key) for key in
                ('software', 'method', 'basis', 'xtb_model', 'nprocs', 'scf_cycles')}
    combined_identity = separated_reactant_identity(reactants[0], reactants[1])
    request = {'schema_version': 1, 'fragments': request_rows, 'physical_method': physical,
               'combined_reactant_identity': combined_identity,
               'optimization': {'max_steps': options.get('product_relaxation_max_steps', 200),
                                 'fmax_ev_per_angstrom': options.get('product_relaxation_fmax', .05)},
               'imaginary_frequency_threshold_cm-1': options.get('imaginary_frequency_threshold', 20.0)}
    request_hash = hashlib.sha256(json.dumps(request, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    manifest = output / 'reference.json'
    if manifest.is_file():
        saved = json.loads(manifest.read_text())
        if saved.get('request_sha256') != request_hash or saved.get('request') != request:
            raise ValueError('Separated-reactant reference settings or inputs changed; use a new characterization directory')
        if saved.get('status') == 'complete':
            for item in saved['fragments']:
                path = output / item['optimized_geometry']
                if not path.is_file() or _sha256(path) != item['optimized_geometry_sha256']:
                    raise ValueError(f"Separated-reactant artifact changed: {path}")
            return saved

    record = {'schema_version': 1, 'request': request, 'request_sha256': request_hash,
              'status': 'running', 'reference_energy_hartree': None, 'fragments': []}
    _write_json(manifest, record)
    try:
        for index, (molecule, source, fragment_request) in enumerate(zip(reactants, sources, request_rows), start=1):
            label = f'fragment_{index}'
            settings = dict(physical, charge=fragment_request['charge'],
                            multiplicity=fragment_request['multiplicity'],
                            scftype=fragment_request['scftype'], geometry_optimizer='geometric', gamma=0.0)
            calculator = PyarGeometricCalculator(qc_params=settings)
            frames, energies, converged = _optimize_geometry(
                molecule.atoms_list, np.asarray(molecule.coordinates, dtype=float), calculator,
                output, label, request['optimization']['max_steps'],
                fmax=request['optimization']['fmax_ev_per_angstrom'],
            )
            if not converged:
                raise RuntimeError(f'isolated reactant fragment {index} did not converge')
            optimized = output / f'{label}_optimized.xyz'
            frequency = _frequency(
                molecule.atoms_list, frames[-1], calculator, output, label,
                request['imaginary_frequency_threshold_cm-1'],
            )
            if not frequency['minimum_confirmed']:
                raise RuntimeError(f'isolated reactant fragment {index} is not a confirmed minimum')
            final_coordinates = frequency['evaluated_coordinates_angstrom']
            write_xyz(molecule.atoms_list, final_coordinates, optimized,
                      energy=frequency['energy_hartree'], precision=12)
            initial_identity = molecule_identity_from_xyz(source)
            final_identity = molecule_identity_from_xyz(optimized)
            if not same_molecular_identity(initial_identity, final_identity):
                raise RuntimeError(f'isolated reactant fragment {index} changed chemical identity during optimization')
            record['fragments'].append({
                'fragment': index, 'input_sha256': fragment_request['source_sha256'],
                'optimized_geometry': optimized.name,
                'optimized_geometry_sha256': _sha256(optimized),
                'charge': fragment_request['charge'], 'multiplicity': fragment_request['multiplicity'],
                'identity': final_identity, 'energy_hartree': frequency['energy_hartree'],
                'optimization_converged': converged, 'frequency': frequency,
            })
            record['status'] = 'running'
            _write_json(manifest, record)
        record['reference_energy_hartree'] = sum(item['energy_hartree'] for item in record['fragments'])
        record['status'] = 'complete'
    except Exception as exc:
        record.update(status='failed', error=f'{type(exc).__name__}: {exc}')
    _write_json(manifest, record)
    return record


def classify_association_endpoints(symbols, endpoint_geometries, product_path,
                                   reactant_reference, reactants, tolerance,
                                   reference_directory=None):
    """Assign IRC branches to product or isolated reactants independent of pose."""
    from pyar.neb import _aligned_rmsd, _bond_set, read_xyz
    from pyar.reaction_identity import molecule_identity_from_xyz, same_molecular_identity
    import tempfile

    if reactant_reference.get('status') != 'complete':
        return {'status': 'reactant_reference_failed', 'endpoint_assignment': None,
                'reactant_product_connection_confirmed': False}
    product_symbols, product_coordinates, _ = read_xyz(product_path)
    product_identity = molecule_identity_from_xyz(product_path)
    reactant_identity = reactant_reference.get('request', {}).get('combined_reactant_identity')
    if not reactant_identity:
        raise ValueError('Separated-reactant reference has no combined reactant identity')
    fragment_lengths = [len(molecule.atoms_list) for molecule in reactants]
    if (len(reactants) != 2 or len(reactant_reference.get('fragments', [])) != len(reactants)
            or sum(fragment_lengths) != len(symbols) or product_symbols != symbols):
        raise ValueError('Association endpoint atoms do not match the saved fragment mapping')
    boundaries = np.cumsum([0, *fragment_lengths])
    for index, molecule in enumerate(reactants):
        if list(symbols[boundaries[index]:boundaries[index + 1]]) != list(molecule.atoms_list):
            raise ValueError('Association endpoint atom order differs from the saved fragment mapping')
    reference_directory = Path(reference_directory or Path(product_path).parent)
    reference_geometries = [read_xyz(reference_directory / item['optimized_geometry'])[1]
                            for item in reactant_reference['fragments']]
    branch_classes = []
    branch_details = []
    for branch, coordinates in enumerate(endpoint_geometries):
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape != (len(symbols), 3) or not np.isfinite(coordinates).all():
            raise ValueError('Malformed IRC endpoint geometry')
        # Use an isolated temporary directory so identity checks leave no run artifacts.
        with tempfile.TemporaryDirectory(prefix='pyar_association_identity_') as directory:
            from pyar.backends import write_xyz as _write_xyz
            full_path = Path(directory) / 'endpoint.xyz'
            _write_xyz(symbols, coordinates, full_path)
            whole_identity = molecule_identity_from_xyz(full_path)
            product_identity_match = same_molecular_identity(product_identity, whole_identity)
            product_match = (product_identity_match
                             and _aligned_rmsd(coordinates, product_coordinates) <= tolerance)
            fragment_matches, fragment_rmsds = [], []
            for index, molecule in enumerate(reactants):
                lower, upper = int(boundaries[index]), int(boundaries[index + 1])
                piece = Path(directory) / f'fragment_{index + 1}.xyz'
                _write_xyz(molecule.atoms_list, coordinates[lower:upper], piece)
                identity = molecule_identity_from_xyz(piece)
                expected = reactant_reference['fragments'][index]['identity']
                fragment_matches.append(same_molecular_identity(expected, identity))
                fragment_rmsds.append(_aligned_rmsd(coordinates[lower:upper], reference_geometries[index]))
            # A stable encounter complex may acquire an interfragment bond
            # under XYZ bond perception. Reactant-side assignment therefore
            # rests on the separately mapped fragment identities and their
            # internal geometries. Every conformer of the product is excluded,
            # even when its geometry falls outside the product RMSD tolerance.
            reactant_match = (not product_identity_match and all(fragment_matches)
                              and all(value <= tolerance for value in fragment_rmsds))
        if product_match and not same_molecular_identity(product_identity, reactant_identity):
            label = 'product'
        elif reactant_match and not product_match:
            interfragment = any((i < boundaries[1] <= j) or (j < boundaries[1] <= i)
                                for i, j in _bond_set(symbols, coordinates)) if len(reactants) == 2 else False
            label = 'stable_pre_reactive_minimum' if interfragment else 'separated_fragments_minima'
        elif product_match and reactant_match:
            label = 'ambiguous'
        else:
            label = 'unmatched'
        branch_classes.append(label)
        branch_details.append({'branch_index': branch, 'classification': label,
                               'fragment_identity_matches': fragment_matches,
                               'fragment_aligned_rmsds_angstrom': fragment_rmsds,
                               'product_identity_and_geometry_match': product_match})
    assignments = [(i, j) for i, first in enumerate(branch_classes)
                   for j, second in enumerate(branch_classes) if i != j
                   and first in {'stable_pre_reactive_minimum', 'separated_fragments_minima'}
                   and second == 'product']
    mapping = assignments[0] if len(assignments) == 1 else None
    return {'status': 'matched' if mapping is not None else 'unmatched',
            'branch_results': branch_details,
            'reactant_product_assignment': list(mapping) if mapping is not None else None,
            'reactant_product_connection_confirmed': mapping is not None}


def classify_reactant_relaxation_failure(stage_result, execution_directory, product_path):
    """Diagnose a failed association reactant relaxation without guessing minima."""
    error = str(stage_result.get("error") or "")
    lowered = error.lower()
    if "did not converge" in lowered:
        return {"status": "reactant_relaxation_nonconverged",
                "classification": "nonconverged", "error": error,
                "reactant_product_connection_confirmed": False}
    if "connectivity changed" not in lowered:
        return {"status": "reactant_relaxation_failed", "classification": "unclassified_failure",
                "error": error, "reactant_product_connection_confirmed": False}
    from pyar.backends import write_xyz
    from pyar.neb import read_xyz
    from pyar.reaction_identity import molecule_identity_from_xyz, same_molecular_identity
    from pyar.xyz_trajectory import parse_xyz_trajectory
    trajectory = Path(execution_directory) / "reactant_path.xyz"
    try:
        frames = parse_xyz_trajectory(trajectory)
        symbols, coordinates = frames[-1][0], frames[-1][1]
        product_symbols, _product_coordinates, _ = read_xyz(product_path)
        if list(symbols) != list(product_symbols):
            return {"status": "reactant_topology_changed_unclassified",
                    "classification": "topology_changed", "trajectory": str(trajectory),
                    "reactant_product_connection_confirmed": False}
        with tempfile.TemporaryDirectory(prefix="pyar_relax_identity_") as directory:
            final_path = Path(directory) / "relaxed_reactant.xyz"
            write_xyz(symbols, coordinates, final_path)
            product_identity = molecule_identity_from_xyz(product_path)
            final_identity = molecule_identity_from_xyz(final_path)
        is_product = same_molecular_identity(product_identity, final_identity)
        return {"status": ("product_formed_during_reactant_relaxation" if is_product
                           else "reactant_topology_changed_unclassified"),
                "classification": ("product" if is_product else "topology_changed"),
                "trajectory": str(trajectory), "product_identity_match": bool(is_product),
                "error": error, "reactant_product_connection_confirmed": False}
    except Exception as exc:
        return {"status": "reactant_relaxation_outcome_unresolved",
                "classification": "unknown", "trajectory": str(trajectory),
                "error": f"{type(exc).__name__}: {exc}",
                "reactant_product_connection_confirmed": False}
