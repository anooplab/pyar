"""One isolated attempt; parent owns publication, budgets, locks and retries."""
from pathlib import Path
import sys

from pyar.core.molecule import Molecule
from pyar.optimization_service import optimize_once
from pyar.persistence import read_json, write_json


def main():
    value = read_json('input.json')['request']
    settings = value['settings']
    molecule = Molecule(value['symbols'], value['coordinates_angstrom'], name=value['name'],
        charge=settings['charge'], multiplicity=settings['multiplicity'])
    molecule.scftype = settings['scftype']
    molecule.fragments = settings.get('fragments', [])
    result = optimize_once(molecule, settings)
    write_json('result.json', result.to_dict(), overwrite=False)


if __name__ == '__main__':
    main()
