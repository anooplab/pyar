"""Generate independently labelled, constructed clustering-distance fixtures.

Labels describe deliberately constructed geometry families, not PES basins.
No comparator output is used to assign a label or transformation witness.
"""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def build_dataset():
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    rng = np.random.default_rng(20261001)
    water = np.array([[0., 0., 0.], [.9572, 0., 0.], [-.239, .927, 0.]])
    angles = np.arange(6) * np.pi / 3
    benzene = np.vstack((np.column_stack((1.397 * np.cos(angles), 1.397 * np.sin(angles), np.zeros(6))),
                         np.column_stack((2.477 * np.cos(angles), 2.477 * np.sin(angles), np.zeros(6)))))
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    if AllChem.EmbedMolecule(molecule, randomSeed=31) < 0 or AllChem.MMFFOptimizeMolecule(molecule) != 0:
        raise RuntimeError("Could not construct the butane template")
    pools, frames = [], []

    def add_pool(name, system_type, variants):
        records, witnesses = [], []
        for index, (family, atoms, coordinates) in enumerate(variants):
            identifier = f"{name}_{index:02d}"
            frames.append((identifier, list(atoms), coordinates))
            records.append({"id": identifier, "family": family, "energy": float(index)})
            order = rng.permutation(len(atoms))
            rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            rotation[:, -1] *= np.linalg.det(rotation)
            translation = rng.uniform(-5, 5, size=3)
            transformed = coordinates[order] @ rotation + translation
            transform_id = identifier + "_permuted"
            frames.append((transform_id, np.asarray(atoms)[order].tolist(), transformed))
            records.append({"id": transform_id, "family": family, "energy": float(index) + .1})
            witnesses.append({"first": identifier, "second": transform_id,
                              "order": order.tolist(), "rotation": rotation.tolist(),
                              "translation": translation.tolist()})
        pools.append({"id": name, "system_type": system_type, "records": records,
                      "rigid_transform_witnesses": witnesses})

    variants = []
    for family, base in (("compact", 3.4), ("expanded", 6.4)):
        for delta in (-.1, 0, .1):
            variants.append((family, ["O", "H", "H"] * 2,
                             np.vstack((water, water + [base + delta, 0, 0]))))
    add_pool("water_dimer_packing", "molecular-aggregate", variants)

    variants = []
    for family, base in (("separated", 15.), ("more-separated", 20.)):
        for delta in (-.1, 0, .1):
            variants.append((family, ["O", "H", "H"] * 2,
                             np.vstack((water, water + [base + delta, 0, 0]))))
    add_pool("water_dimer_beyond_cutoff", "molecular-aggregate", variants)

    variants = []
    for family, base in (("co-oriented", 0), ("turned", 110)):
        for delta in (-10, 0, 10):
            angle = np.deg2rad(base + delta)
            rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                                 [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
            variants.append((family, ["O", "H", "H"] * 2,
                             np.vstack((water, water @ rotation + [3.5, 0, 0]))))
    add_pool("water_dimer_orientation", "molecular-aggregate", variants)

    ammonia = np.array([[0., 0., 0.], [.95, 0., -.35],
                        [-.475, .823, -.35], [-.475, -.823, -.35]])
    variants = []
    for family, base in (("compact", 3.4), ("expanded", 6.4)):
        for delta in (-.1, 0, .1):
            variants.append((family, ["O", "H", "H", "N", "H", "H", "H"],
                             np.vstack((water, ammonia + [base + delta, 0, 0]))))
    add_pool("water_ammonia_packing", "molecular-aggregate", variants)

    variants = []
    for family in ("triangle", "chain"):
        for spacing in (3.4, 3.5, 3.6):
            centers = ([[0, 0, 0], [spacing, 0, 0], [spacing / 2, spacing * np.sqrt(3) / 2, 0]]
                       if family == "triangle" else [[-spacing, 0, 0], [0, 0, 0], [spacing, 0, 0]])
            variants.append((family, ["O", "H", "H"] * 3,
                             np.vstack([water + center for center in centers])))
    add_pool("water_trimer_packing", "molecular-aggregate", variants)

    variants = []
    for family, shift in (("stacked", 0.), ("slipped", 3.)):
        for height in (3.4, 3.5, 3.6):
            variants.append((family, ["C"] * 6 + ["H"] * 6 + ["C"] * 6 + ["H"] * 6,
                             np.vstack((benzene, benzene + [shift, 0, height]))))
    add_pool("benzene_dimer_packing", "molecular-aggregate", variants)

    variants = []
    for family, torsions in (("anti", (170, 180, 190)), ("gauche", (50, 60, 70))):
        for torsion in torsions:
            conformer = Chem.Mol(molecule)
            rdMolTransforms.SetDihedralDeg(conformer.GetConformer(), 0, 1, 2, 3, torsion)
            variants.append((family, [atom.GetSymbol() for atom in conformer.GetAtoms()],
                             np.asarray(conformer.GetConformer().GetPositions())))
    add_pool("butane_torsions", "molecular", variants)
    return {"schema_version": 1, "seed": 20261001,
            "structures_file": "structures.xyz", "pools": pools,
            "label_interpretation": "constructed packing/torsion families; not optimized PES basins",
            "energy_interpretation": "synthetic ranking metadata only"}, frames


def write_dataset(root=ROOT):
    manifest, frames = build_dataset()
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    lines = []
    for identifier, atoms, coordinates in frames:
        lines.extend([str(len(atoms)), "structure_id=" + identifier])
        lines.extend(f"{atom} {x:.12f} {y:.12f} {z:.12f}"
                     for atom, (x, y, z) in zip(atoms, coordinates))
    (root / "structures.xyz").write_text("\n".join(lines) + "\n")
    (root / "dataset.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    write_dataset()
