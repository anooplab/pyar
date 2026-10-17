"""Opt-in atom-contact evidence and operational solvent-shell labels."""

from __future__ import annotations

import numpy as np


def assess_contacts(molecule, request):
    """Measure explicit-distance contacts without assigning chemical bonds.

    Site occupancy counts solvent molecules with at least one atom pair under
    the user-selected distance. A molecule may count at multiple named sites.
    """
    if not request.site_groups:
        return None
    solute_count = len(request.solute)
    atoms = list(molecule.atoms_list)
    coordinates = np.asarray(molecule.coordinates, dtype=float)
    solvent_fragments = tuple(tuple(int(i) for i in part)
                              for part in getattr(molecule, "solvent_fragments", ()))
    if (coordinates.shape != (len(atoms), 3) or not np.isfinite(coordinates).all()
            or tuple(getattr(molecule, "solute_atom_indices", ())) != tuple(range(solute_count))
            or atoms[:solute_count] != list(request.solute.atoms_list)):
        return {"status": "unknown", "reason": "atom coordinates or original-solute mapping unavailable"}
    solvent_indices = [i for fragment in solvent_fragments for i in fragment]
    expected = list(range(solute_count, len(atoms)))
    if sorted(solvent_indices) != expected or len(solvent_indices) != len(set(solvent_indices)):
        return {"status": "unknown", "reason": "solvent fragment identity is incomplete or overlapping"}

    cutoff = float(request.contact_distance)

    # Report changes in a covalent-radius proximity graph without treating
    # this heuristic as a bond-order or reaction assignment.
    from pyar.property import get_bond_matrix
    solute_radii = np.asarray(request.solute.covalent_radius, dtype=float)
    initial_bonds = get_bond_matrix(np.asarray(request.solute.coordinates), solute_radii)
    observed_bonds = get_bond_matrix(coordinates[:solute_count], solute_radii)
    original_pairs = {tuple(int(i) for i in pair) for pair in np.argwhere(np.triu(initial_bonds, 1))}
    observed_pairs = {tuple(int(i) for i in pair) for pair in np.argwhere(np.triu(observed_bonds, 1))}

    def pairs(left, right):
        if not left or not right:
            return []
        delta = coordinates[np.asarray(left), None, :] - coordinates[np.asarray(right)][None, :, :]
        distances = np.linalg.norm(delta, axis=2)
        found = []
        for i, j in np.argwhere(distances <= cutoff):
            found.append({"atom_a": int(left[i]), "atom_b": int(right[j]),
                          "distance_angstrom": round(float(distances[i, j]), 6)})
        return sorted(found, key=lambda item: (item["distance_angstrom"], item["atom_a"], item["atom_b"]))

    solute_indices = list(range(solute_count))
    by_fragment = []
    first = []
    for fragment_index, fragment in enumerate(solvent_fragments):
        contact_pairs = [{"solute_atom": pair["atom_a"], "solvent_atom": pair["atom_b"],
                          "distance_angstrom": pair["distance_angstrom"]}
                         for pair in pairs(solute_indices, list(fragment))]
        is_first = bool(contact_pairs)
        if is_first:
            first.append(fragment_index)
        by_fragment.append({"solvent_index": fragment_index,
                            "atom_indices": list(fragment),
                            "shell": "first" if is_first else "pending",
                            "solute_contacts": contact_pairs})
    for fragment_index, fragment in enumerate(solvent_fragments):
        if by_fragment[fragment_index]["shell"] == "first":
            continue
        network_contacts = []
        for first_index in first:
            contact_list = pairs(list(solvent_fragments[first_index]), list(fragment))
            for contact in contact_list:
                network_contacts.append({"first_shell_solvent": first_index,
                                         "first_shell_atom": contact["atom_a"],
                                         "solvent_atom": contact["atom_b"],
                                         "distance_angstrom": contact["distance_angstrom"]})
        by_fragment[fragment_index]["shell"] = "second" if network_contacts else "outer"
        by_fragment[fragment_index]["first_shell_contacts"] = sorted(
            network_contacts, key=lambda item: (item["distance_angstrom"], item["first_shell_solvent"],
                                                 item["first_shell_atom"], item["solvent_atom"]))

    sites = []
    for name, indices in request.site_groups:
        occupied = []
        per_solvent = []
        for record, fragment in zip(by_fragment, solvent_fragments):
            contact_pairs = [{"site_atom": pair["atom_a"], "solvent_atom": pair["atom_b"],
                              "distance_angstrom": pair["distance_angstrom"]}
                             for pair in pairs(list(indices), list(fragment))]
            if contact_pairs:
                occupied.append(record["solvent_index"])
            per_solvent.append({"solvent_index": record["solvent_index"],
                                "contacting": bool(contact_pairs), "contacts": contact_pairs})
        sites.append({"name": name, "solute_atom_indices": list(indices),
                      "occupancy_count": len(occupied), "occupied_by_solvents": occupied,
                      "solvent_contacts": per_solvent})
    return {
        "status": "measured",
        "model": "explicit-atom-distance-v1",
        "distance_cutoff_angstrom": cutoff,
        "interpretation": "geometric contacts; not evidence of bond formation or thermodynamic shell membership",
        "solute_connectivity": {
            "status": "changed" if original_pairs != observed_pairs else "unchanged",
            "model": "covalent-radii-distance-factor-1.3-v1",
            "interpretation": "connectivity warning heuristic only; not bond-order or reaction evidence",
            "original_pairs": [list(pair) for pair in sorted(original_pairs)],
            "observed_pairs": [list(pair) for pair in sorted(observed_pairs)],
            "added_pairs": [list(pair) for pair in sorted(observed_pairs - original_pairs)],
            "removed_pairs": [list(pair) for pair in sorted(original_pairs - observed_pairs)],
        },
        "whole_solute_first_shell_count": len(first),
        "whole_solute_second_shell_count": sum(item["shell"] == "second" for item in by_fragment),
        "whole_solute_outer_count": sum(item["shell"] == "outer" for item in by_fragment),
        "solvent_molecules": by_fragment,
        "sites": sites,
    }
