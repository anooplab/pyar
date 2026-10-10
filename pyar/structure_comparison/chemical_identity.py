"""Optional, in-memory XYZ chemical perception using RDKit's xyz2mol approach.

Charge is supplied or explicitly assumed, never searched to obtain a result.
This chemical graph is independent of PyAR's coordinate-adjacency graph.
"""

from numbers import Integral


INSTALL_HINT = 'pip install "pyar-chem[identity]"'
PERCEPTION_METHOD = 'rdkit-determine-bonds'


def resolve_charge(charge=None, *, known_charge=None):
    """Resolve explicit, reliably known, or assumed charge in that order.

    Callers must only pass ``known_charge`` from authoritative input/state;
    the Molecule constructor's default neutral state is not such evidence.
    """
    value, source = ((charge, 'explicit') if charge is not None else
                     (known_charge, 'known') if known_charge is not None else (0, 'assumed'))
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError('Chemical-perception charge must be an integer')
    return int(value), source


def _rdkit_api():
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdDetermineBonds
    return Chem, rdDetermineBonds, rdBase


def perceive_chemical_identity(molecule, charge=None, *, known_charge=None):
    """Return canonical isomeric SMILES and diagnostics, without file output.

    Geometry must already be validated. Optional dependency/perception failures
    are represented in the result, not substituted with guessed chemistry.
    """
    charge_used, charge_source = resolve_charge(charge, known_charge=known_charge)
    result = {
        'canonical_smiles': None, 'charge_used': charge_used,
        'charge_source': charge_source, 'perception_success': False,
        'perception_method': PERCEPTION_METHOD, 'status': 'unavailable',
        'reason': None, 'installation_hint': None,
    }
    try:
        Chem, rdDetermineBonds, rdBase = _rdkit_api()
    except (ImportError, OSError) as exc:
        result.update(reason=f'RDKit bond perception is unavailable: {exc}',
                      installation_hint=INSTALL_HINT)
        return result
    # Fixed-point coordinates avoid RDKit XYZ readers rejecting exponent notation.
    xyz_block = f'{len(molecule.atoms_list)}\n\n' + ''.join(
        f'{symbol} {x:.17f} {y:.17f} {z:.17f}\n'
        for symbol, (x, y, z) in zip(molecule.atoms_list, molecule.coordinates))
    try:
        # Preserve failure diagnostics in the result, without RDKit log noise.
        with rdBase.BlockLogs():
            perceived = Chem.MolFromXYZBlock(xyz_block)
            if perceived is None:
                raise ValueError('RDKit could not parse the XYZ geometry')
            rdDetermineBonds.DetermineBonds(perceived, charge=charge_used)
            smiles = Chem.MolToSmiles(Chem.RemoveHs(perceived), canonical=True,
                                     isomericSmiles=True)
            if not smiles:
                raise ValueError('RDKit returned an empty canonical SMILES')
        result.update(canonical_smiles=smiles, perception_success=True, status='success')
    except Exception as exc:
        # This boundary deliberately catches RDKit's varied chemistry exceptions.
        # Input validation and the independent coordinate comparison remain outside.
        result.update(status='failed', reason=f'Bond-order perception failed: {exc}')
    return result
