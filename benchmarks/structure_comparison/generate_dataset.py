"""Generate the compact structure-comparison fixture dataset.

Requires the optional RDKit dependency: ``pip install 'pyar-chem[conformer]'``.
The checked-in XYZ/CSV files are the benchmark inputs; this script records how
they were generated so rebuilding them is not required to run comparisons.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski, rdMolTransforms


ROOT = Path(__file__).resolve().parent
XYZ_PATH = ROOT / "structures.xyz"
PAIRS_PATH = ROOT / "pairs.csv"
SEED = 20260929

# Formula-matched constitutional-isomer pairs, with straightforward neutral
# valence. The labels describe graph topology, independent of a distance score.
ISOMER_FAMILIES = {
    "C4H10": ("CCCC", "CC(C)C"),
    "C2H6O": ("CCO", "COC"),
    "C3H8O": ("CCCO", "CC(O)C"),
    "C3H6O": ("CC(=O)C", "CCC=O"),
    "C4H8O_ring_chain": ("C1CCCO1", "CCCC=O"),
    "C4H10O": ("CCCCO", "CCC(O)C"),
    "C5H12": ("CCCCC", "CC(C)CC"),
    "C4H8": ("C1CCC1", "C=CCC"),
}


def _stable_seed(text):
    digest = hashlib.sha256(text.encode()).digest()
    return int.from_bytes(digest[:4], "little") & 0x7FFFFFFF


def _rotatable_torsions(molecule):
    torsions = []
    for left, right in molecule.GetSubstructMatches(Lipinski.RotatableBondSmarts):
        left_neighbors = sorted(
            atom.GetIdx() for atom in molecule.GetAtomWithIdx(left).GetNeighbors()
            if atom.GetIdx() != right
        )
        right_neighbors = sorted(
            atom.GetIdx() for atom in molecule.GetAtomWithIdx(right).GetNeighbors()
            if atom.GetIdx() != left
        )
        if left_neighbors and right_neighbors:
            torsions.append((left_neighbors[0], left, right, right_neighbors[0]))
    return torsions


def _torsion_signature(molecule, conformer_id):
    conformer = molecule.GetConformer(conformer_id)
    signature = []
    for torsion in _rotatable_torsions(molecule):
        angle = rdMolTransforms.GetDihedralDeg(conformer, *torsion) % 360.0
        signature.append(int(round(angle / 45.0)) % 8)
    return tuple(signature)


def _generate(smiles, name):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Could not parse SMILES for {name}: {smiles}")
    molecule = Chem.AddHs(molecule)
    params = AllChem.ETKDGv3()
    params.randomSeed = _stable_seed(f"{SEED}:{name}")
    params.numThreads = 1
    conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=24, params=params))
    if not conformer_ids:
        raise RuntimeError(f"RDKit failed to embed {name}")
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        results = AllChem.MMFFOptimizeMoleculeConfs(molecule, numThreads=1, maxIters=500)
        energies = {cid: float(result[1]) for cid, result in zip(conformer_ids, results) if result[0] == 0}
        conformer_ids = [cid for cid in conformer_ids if cid in energies]
        if not conformer_ids:
            raise RuntimeError(f"No converged MMFF geometries for {name}")
    else:
        raise RuntimeError(f"Missing MMFF parameters for {name}")

    torsions = _rotatable_torsions(molecule)
    selected = []
    signatures = set()
    for cid in sorted(conformer_ids, key=lambda item: (energies[item], item)):
        signature = _torsion_signature(molecule, cid)
        if torsions and signature in signatures:
            continue
        selected.append(cid)
        signatures.add(signature)
        if len(selected) == 3:
            break
    if len(selected) < 2 and len(conformer_ids) > 1:
        selected = sorted(conformer_ids, key=lambda item: (energies[item], item))[:2]
    return molecule, selected, energies


def _coordinates(molecule, conformer_id):
    conformer = molecule.GetConformer(conformer_id)
    return np.asarray(conformer.GetPositions(), dtype=float)


def _write_xyz_frame(stream, structure_id, symbols, coordinates, metadata):
    stream.write(f"{len(symbols)}\n")
    detail = " ".join(f"{key}={value}" for key, value in metadata.items())
    stream.write(f"structure_id={structure_id} {detail}\n")
    for symbol, (x, y, z) in zip(symbols, coordinates):
        stream.write(f"{symbol:3s} {x: .9f} {y: .9f} {z: .9f}\n")


def _random_rotation(rng):
    matrix = rng.normal(size=(3, 3))
    q, r = np.linalg.qr(matrix)
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def _au13_motifs():
    """Return idealized icosahedral and cuboctahedral 13-site geometries."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    ico = []
    for x, y in itertools.product((-1.0, 1.0), repeat=2):
        ico.extend(((0.0, x, y * phi), (x, y * phi, 0.0), (x * phi, 0.0, y)))
    ico = np.unique(np.asarray(ico, dtype=float), axis=0)
    ico /= np.linalg.norm(ico, axis=1).mean()
    ico *= 3.8

    # Permutations of (±1, ±1, 0) form the cuboctahedral shell.
    cuboct = set()
    for zero_axis in range(3):
        nonzero_axes = [axis for axis in range(3) if axis != zero_axis]
        for signs in itertools.product((-1.0, 1.0), repeat=2):
            point = [0.0, 0.0, 0.0]
            point[nonzero_axes[0]], point[nonzero_axes[1]] = signs
            cuboct.add(tuple(point))
    cuboct = np.asarray(sorted(cuboct), dtype=float)
    cuboct /= np.linalg.norm(cuboct, axis=1).mean()
    cuboct *= 3.8
    return np.vstack((np.zeros((1, 3)), ico)), np.vstack((np.zeros((1, 3)), cuboct))


def _water_monomer(angle=0.0, origin=(0.0, 0.0, 0.0)):
    """Return one rigid H2O monomer with a chosen in-plane orientation."""
    base = np.asarray([[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0], [-0.2390, 0.9270, 0.0]])
    cosine, sine = np.cos(angle), np.sin(angle)
    rotation = np.asarray([[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]])
    return base @ rotation.T + np.asarray(origin, dtype=float)


def _water_clusters():
    """Construct labelled water cluster arrangements without energy claims."""
    hbond = np.vstack((
        _water_monomer(0.0, (0.0, 0.0, 0.0)),
        _water_monomer(0.0, (2.8, 0.0, 0.0)),
    ))
    side = np.vstack((
        _water_monomer(0.0, (0.0, 0.0, 0.0)),
        _water_monomer(np.pi / 2.0, (0.0, 3.6, 0.0)),
    ))

    ring_positions = np.asarray([[0.0, 0.0, 0.0], [2.8, 0.0, 0.0], [1.4, 2.424871, 0.0]])
    ring = []
    for index, origin in enumerate(ring_positions):
        target = ring_positions[(index + 1) % len(ring_positions)] - origin
        ring.append(_water_monomer(np.arctan2(target[1], target[0]), origin))
    chain = np.vstack((
        _water_monomer(0.0, (0.0, 0.0, 0.0)),
        _water_monomer(0.0, (2.8, 0.0, 0.0)),
        _water_monomer(0.0, (5.6, 0.0, 0.0)),
    ))
    return {
        "water_dimer_hbond": hbond,
        "water_dimer_side": side,
        "water_trimer_ring": np.vstack(ring),
        "water_trimer_chain": chain,
    }


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    structures = {}
    molecule_records = {}
    with XYZ_PATH.open("w") as xyz:
        for family, smiles_pair in ISOMER_FAMILIES.items():
            for isomer_index, smiles in enumerate(smiles_pair):
                name = f"{family}_{isomer_index + 1}"
                molecule, conformers, energies = _generate(smiles, name)
                symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
                for rank, conformer_id in enumerate(conformers, start=1):
                    structure_id = f"{name}_conf{rank}"
                    coords = _coordinates(molecule, conformer_id)
                    structures[structure_id] = (symbols, coords)
                    molecule_records[name] = {
                        "smiles": smiles,
                        "molecule": molecule,
                        "symbols": symbols,
                        "conformer_ids": conformers,
                        "energies": energies,
                        "family": family,
                        "torsion_signatures": {
                            rank: _torsion_signature(molecule, cid)
                            for rank, cid in enumerate(conformers, start=1)
                        },
                        "has_rotatable_torsion": bool(_rotatable_torsions(molecule)),
                    }
                    _write_xyz_frame(xyz, structure_id, symbols, coords, {
                        "family": family,
                        "smiles": smiles,
                        "generator": "RDKit-ETKDGv3-MMFF",
                        "seed": _stable_seed(f"{SEED}:{name}"),
                        "conformer": rank,
                        "mmff_converged": "true",
                        "mmff_energy_kcal_mol": f"{energies[conformer_id]:.8f}",
                    })

        # Exact same-geometry pairs after arbitrary rigid transforms and atom
        # order changes. The relation label follows directly from construction.
        rng = np.random.default_rng(SEED)
        base_names = sorted(molecule_records)
        for index in range(24):
            base_id = f"{base_names[index % len(base_names)]}_conf1"
            symbols, coords = structures[base_id]
            permutation = rng.permutation(len(symbols))
            transformed = (coords - coords.mean(axis=0)) @ _random_rotation(rng)
            transformed += rng.normal(size=3) * 4.0
            transformed_id = f"transform_{index:03d}"
            structures[transformed_id] = ([symbols[i] for i in permutation], transformed[permutation])
            _write_xyz_frame(xyz, transformed_id, *structures[transformed_id], {
                "source_structure": base_id,
                "augmentation": "rigid_transform_and_atom_permutation",
                "source_atom_order": ",".join(map(str, permutation)),
            })

        # Small unoptimized Cartesian perturbations are intentionally labeled
        # ambiguous: they test threshold sensitivity, not a basin assertion.
        for index in range(12):
            base_id = f"{base_names[index % len(base_names)]}_conf1"
            symbols, coords = structures[base_id]
            noisy_id = f"nearby_ambiguous_{index:03d}"
            noisy = coords + rng.normal(scale=0.08, size=coords.shape)
            structures[noisy_id] = (symbols, noisy)
            _write_xyz_frame(xyz, noisy_id, symbols, noisy, {
                "source_structure": base_id,
                "augmentation": "cartesian_noise_sigma_0.08_angstrom",
            })

        # Idealized cluster motifs are deliberately not described as optimized
        # minima. They exercise permutation, composition, and motif comparisons.
        ico, cuboct = _au13_motifs()
        au_symbols = ["Au"] * 13
        ar_symbols = ["Ar"] * 13
        cluster_geometries = {
            "au13_ico": (au_symbols, ico),
            "au13_cubocta": (au_symbols, cuboct),
            "ar13_ico": (ar_symbols, ico),
        }
        for structure_id, (symbols, coords) in cluster_geometries.items():
            structures[structure_id] = (symbols, coords)
            _write_xyz_frame(xyz, structure_id, symbols, coords, {
                "system_class": "atomic_cluster",
                "construction": "idealized_13_site_motif",
                "geometry": structure_id.rsplit("_", 1)[-1],
            })
        extra_rng = np.random.default_rng(SEED + 1)
        permutation = extra_rng.permutation(len(au_symbols))
        au_rotated = (ico - ico.mean(axis=0)) @ _random_rotation(extra_rng)
        au_rotated += np.asarray([8.0, -3.0, 2.0])
        structures["au13_ico_rotperm"] = ( [au_symbols[i] for i in permutation], au_rotated[permutation])
        _write_xyz_frame(xyz, "au13_ico_rotperm", *structures["au13_ico_rotperm"], {
            "system_class": "atomic_cluster",
            "source_structure": "au13_ico",
            "augmentation": "rigid_transform_and_atom_permutation",
                "source_atom_order": ",".join(map(str, permutation)),
        })
        au_distorted = ico + extra_rng.normal(scale=0.18, size=ico.shape)
        structures["au13_ico_distorted"] = (au_symbols, au_distorted)
        _write_xyz_frame(xyz, "au13_ico_distorted", au_symbols, au_distorted, {
            "system_class": "atomic_cluster",
            "source_structure": "au13_ico",
            "augmentation": "cartesian_noise_sigma_0.18_angstrom",
        })

        # Rigid water clusters cover disconnected molecular fragments and
        # non-covalent arrangements. Their labels follow the construction.
        for structure_id, coords in _water_clusters().items():
            count = len(coords) // 3
            symbols = ["O", "H", "H"] * count
            structures[structure_id] = (symbols, coords)
            _write_xyz_frame(xyz, structure_id, symbols, coords, {
                "system_class": "molecular_non_covalent_cluster",
                "construction": "rigid_water_monomers",
                "arrangement": structure_id.removeprefix("water_"),
            })
        swapped_water = np.vstack((structures["water_dimer_hbond"][1][3:], structures["water_dimer_hbond"][1][:3]))
        structures["water_dimer_hbond_fragperm"] = (["O", "H", "H"] * 2, swapped_water)
        _write_xyz_frame(xyz, "water_dimer_hbond_fragperm", *structures["water_dimer_hbond_fragperm"], {
            "system_class": "molecular_non_covalent_cluster",
            "source_structure": "water_dimer_hbond",
            "augmentation": "whole_monomer_permutation",
            "source_atom_order": "3,4,5,0,1,2",
        })

        # Coordinate perturbations are stress probes. Source identity does
        # not establish the physical connectivity of the perturbed geometry.
        distortion_rng = np.random.default_rng(SEED + 2)
        distortion_sources = base_names[:8]
        for source_name in distortion_sources:
            source_id = f"{source_name}_conf1"
            symbols, coords = structures[source_id]
            for sigma in (0.05, 0.15, 0.35, 0.65):
                suffix = str(int(round(sigma * 100))).zfill(2)
                structure_id = f"distorted_{source_name}_noise{suffix}"
                distorted = coords + distortion_rng.normal(scale=sigma, size=coords.shape)
                structures[structure_id] = (symbols, distorted)
                _write_xyz_frame(xyz, structure_id, symbols, distorted, {
                    "system_class": "distorted_molecule",
                    "source_structure": source_id,
                    "augmentation": f"cartesian_noise_sigma_{sigma:.2f}_angstrom",
                })

    pairs = []
    def add_pair(left, right, label, provenance, confidence="high"):
        pairs.append({
            "pair_id": f"pair_{len(pairs) + 1:04d}",
            "structure_a": left,
            "structure_b": right,
            "label": label,
            "confidence": confidence,
            "provenance": provenance,
        })

    # Each exact transform has a known relation independent of any comparator.
    for index in range(24):
        base_id = f"{base_names[index % len(base_names)]}_conf1"
        add_pair(base_id, f"transform_{index:03d}", "same_structure", "constructed_rigid_transform_and_permutation")

    # Formula-matched connectivity isomers: cross all available conformers.
    for family in ISOMER_FAMILIES:
        left = f"{family}_1"
        right = f"{family}_2"
        left_ids = molecule_records[left]["conformer_ids"]
        right_ids = molecule_records[right]["conformer_ids"]
        for i, j in itertools.product(range(1, len(left_ids) + 1), range(1, len(right_ids) + 1)):
            add_pair(
                f"{left}_conf{i}", f"{right}_conf{j}",
                "different_connectivity", f"formula_matched_isomer_family:{family}",
            )

    # Indexed torsion bins can distinguish symmetry-equivalent rotamers.
    # These pairs have no validated geometric-equivalence or basin label.
    for name, record in molecule_records.items():
        conformer_count = len(record["conformer_ids"])
        if not record["has_rotatable_torsion"]:
            continue
        for left, right in itertools.combinations(range(1, conformer_count + 1), 2):
            if record["torsion_signatures"][left] == record["torsion_signatures"][right]:
                continue
            add_pair(
                f"{name}_conf{left}", f"{name}_conf{right}",
                "unreviewed_conformer_pair", "different_indexed_torsion_bins_not_validated_for_symmetry",
                confidence="medium",
            )

    for index in range(12):
        base_id = f"{base_names[index % len(base_names)]}_conf1"
        add_pair(base_id, f"nearby_ambiguous_{index:03d}", "ambiguous", "small_unoptimized_cartesian_perturbation", confidence="low")

    add_pair("au13_ico", "au13_ico_rotperm", "same_structure", "atomic_cluster_rigid_transform_and_permutation")
    add_pair("au13_ico", "au13_cubocta", "different_cluster_motif", "idealized_icosahedral_vs_cuboctahedral_13_site_geometry", confidence="medium")
    add_pair("au13_ico", "au13_ico_distorted", "distorted_cluster", "cartesian_noise_sigma_0.18_angstrom", confidence="medium")
    add_pair("au13_ico", "ar13_ico", "different_composition", "same_icosahedral_geometry_different_element")

    add_pair("water_dimer_hbond", "water_dimer_hbond_fragperm", "same_structure", "whole_water_monomer_permutation")
    add_pair("water_dimer_hbond", "water_dimer_side", "different_cluster_arrangement", "constructed_hbond_vs_side_by_side_water_dimer", confidence="medium")
    add_pair("water_trimer_ring", "water_trimer_chain", "different_cluster_arrangement", "constructed_cyclic_vs_chain_water_trimer", confidence="medium")

    for source_name in distortion_sources:
        source_id = f"{source_name}_conf1"
        for sigma in (0.05, 0.15, 0.35, 0.65):
            suffix = str(int(round(sigma * 100))).zfill(2)
            structure_id = f"distorted_{source_name}_noise{suffix}"
            add_pair(source_id, structure_id, "coordinate_perturbation", f"source_geometry_cartesian_noise_sigma_{sigma:.2f}_angstrom", confidence="medium")

    fieldnames = ["pair_id", "structure_a", "structure_b", "label", "confidence", "provenance"]
    with PAIRS_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(pairs)
    print(f"Wrote {len(structures)} structures and {len(pairs)} labeled pairs")
    print(f"Pair labels: {dict((label, sum(pair['label'] == label for pair in pairs)) for label in sorted({p['label'] for p in pairs}))}")


if __name__ == "__main__":
    main()
