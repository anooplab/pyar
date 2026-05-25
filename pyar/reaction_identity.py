"""Reaction identity helpers for product detection and restart bookkeeping."""

from __future__ import annotations

from pyar.interface import babel


def write_disconnected_reference(molecule, path, separation=100.0):
    """Write an unreacted identity reference with fragments far apart."""
    reference = molecule.copy()
    if len(reference.fragments) == 2:
        reference.coordinates[reference.fragments[1], 0] += separation
    reference.mol_to_xyz(path)


def molecule_identity_from_xyz(xyzfile):
    """Return the OpenBabel identity for an XYZ file."""
    identity = {
        "inchi": babel.make_inchi_string_from_xyz(xyzfile),
        "smiles": babel.make_smile_string_from_xyz(xyzfile),
    }
    if not identity["inchi"] or not identity["smiles"]:
        raise ValueError(f"Could not determine complete product identity from {xyzfile}")
    return identity


def same_molecular_identity(first, second):
    """Return whether two identities represent the same molecular product.

    InChI is used as the stable identity key. SMILES is retained for readable
    reports, but equivalent products need not serialize to the same SMILES.
    """
    return first["inchi"] == second["inchi"]
