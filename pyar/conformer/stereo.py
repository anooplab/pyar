"""Stereochemistry bookkeeping for conformer search.

Conformers of one molecule must share its configuration. ETKDG enforces the
stereo elements a SMILES specifies but would leave unspecified ones random, so
the workflow refuses such input; a backend optimisation can in principle still
invert a centre. These helpers perceive stereo from 3D coordinates and express
it as a canonical isomeric SMILES so records can be compared, filtered and
reported.
"""

import logging

import numpy as np

stereo_logger = logging.getLogger('pyar.conformer.stereo')


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except ImportError:
        return None, None
    return Chem, Point3D


def describe_input_stereo(rdkit_molecule):
    """Return the input's stereo elements and, if fully specified, its reference signature.

    The result is a JSON-safe mapping with ``available`` (False when RDKit or a
    real molecule is unavailable), ``unspecified`` (list of ``{type, atom, label}``)
    and ``reference`` (canonical isomeric SMILES, or None when any stereo
    element is unspecified).
    """
    Chem, _ = _rdkit()
    if Chem is None:
        return {'available': False, 'unspecified': [], 'reference': None}
    try:
        elements = list(Chem.FindPotentialStereo(Chem.Mol(rdkit_molecule)))
        unspecified = [
            {'type': str(element.type), 'atom': int(element.centeredOn),
             'label': _element_label(rdkit_molecule, element, Chem)}
            for element in elements
            if element.specified != Chem.StereoSpecified.Specified
        ]
        reference = None if unspecified else Chem.MolToSmiles(Chem.RemoveHs(rdkit_molecule))
    except Exception as exc:  # fakes in tests, or molecules RDKit cannot canonicalize
        stereo_logger.debug('Stereo perception unavailable: %s', exc)
        return {'available': False, 'unspecified': [], 'reference': None}
    return {'available': True, 'unspecified': unspecified, 'reference': reference}


def _element_label(molecule, element, Chem):
    """Name a stereo element by 1-based atom numbers, e.g. ``C2 centre`` or ``C2=C3 bond``."""
    def atom(index):
        return f"{molecule.GetAtomWithIdx(int(index)).GetSymbol()}{int(index) + 1}"
    if element.type == Chem.StereoType.Bond_Double:
        bond = molecule.GetBondWithIdx(int(element.centeredOn))
        return f"{atom(bond.GetBeginAtomIdx())}={atom(bond.GetEndAtomIdx())} double bond"
    return f"{atom(element.centeredOn)} centre"


def stereo_signature(template, coordinates):
    """Return the canonical isomeric SMILES perceived from ``coordinates``, or None.

    ``template`` supplies the bond graph; ``coordinates`` must follow its atom
    order, as PyAR conformer molecules do.
    """
    Chem, Point3D = _rdkit()
    if Chem is None or template is None or coordinates is None:
        return None
    try:
        positions = np.asarray(coordinates, dtype=float)
        work = Chem.Mol(template)
        if positions.shape != (work.GetNumAtoms(), 3):
            return None
        work.RemoveAllConformers()
        conformer = Chem.Conformer(work.GetNumAtoms())
        for index, (x, y, z) in enumerate(positions):
            conformer.SetAtomPosition(index, Point3D(float(x), float(y), float(z)))
        conformer.Set3D(True)
        work.AddConformer(conformer, assignId=True)
        Chem.AssignStereochemistryFrom3D(work)
        return Chem.MolToSmiles(Chem.RemoveHs(work))
    except Exception as exc:
        stereo_logger.debug('Stereo signature unavailable: %s', exc)
        return None


def conformer_coordinates(rdkit_molecule, conformer_id):
    """Return the Cartesian coordinates of one RDKit conformer, or None."""
    try:
        return np.asarray(rdkit_molecule.GetConformer(int(conformer_id)).GetPositions(), dtype=float)
    except Exception:
        return None
