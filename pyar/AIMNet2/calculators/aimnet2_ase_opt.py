import os
import re

import numpy as np
from pyar.optional_dependencies import optional_dependency_error

pybel = None


def optimize(atoms, prec=1e-3, steps=1000, traj=None):
    try:
        import torch
        from ase.optimize import LBFGS
    except ImportError as exc:
        module_name = getattr(exc, "name", None) or "torch"
        raise optional_dependency_error(module_name, feature="AIMNet2 optimizer", extra="aimnet2") from exc
    with torch.jit.optimized_execution(False):
        opt = LBFGS(atoms, trajectory=traj)
        opt.run(prec, steps)


def pybel2atoms(mol):
    try:
        import ase
    except ImportError as exc:
        raise optional_dependency_error("ase", feature="AIMNet2 optimizer", extra="aimnet2") from exc
    coord = np.array([a.coords for a in mol.atoms])
    numbers = np.array([a.atomicnum for a in mol.atoms])
    atoms = ase.Atoms(positions=coord, numbers=numbers)
    return atoms


def update_mol(mol, atoms, align=True):
    # make copy
    mol_old = pybel.Molecule(pybel.ob.OBMol(mol.OBMol))
    # update coord
    for i, c in enumerate(atoms.get_positions()):
        mol.OBMol.GetAtom(i+1).SetVector(*c.tolist())
    # align
    if align:
        aligner = pybel.ob.OBAlign(False, False)
        aligner.SetRefMol(mol_old.OBMol)
        aligner.SetTargetMol(mol.OBMol)
        aligner.Align()
        rmsd = aligner.GetRMSD()
        aligner.UpdateCoords(mol.OBMol)
        print(f'RMSD: {rmsd:.2f} Angs')


def guess_pybel_type(filename):
    assert '.' in filename
    return os.path.splitext(filename)[1][1:]


def guess_charge(mol):
    m = re.search(r'charge: (-?\d+)', mol.title)
    if m:
        charge = int(m.group(1))
    else:
        charge = mol.charge
    return charge


def main():
    import argparse
    global pybel
    parser = argparse.ArgumentParser()
    parser.add_argument('--charge', '-c', type=int, default=None,
                        help='Molecular charge (default: check molecule title for "charge: {int}" or get total charge from OpenBabel).')
    parser.add_argument('--traj', help='Trajectory file', type=str, default=None)
    parser.add_argument('--fmax', type=float, default=5e-3, help='Optimization threshold.')
    parser.add_argument('model')
    parser.add_argument('in_file')
    parser.add_argument('out_file')
    args = parser.parse_args()

    # PyAR writes XYZ and supplies an explicit charge. This route needs only
    # the pip-installable AIMNet2 stack; retain format/charge perception below
    # for direct callers requesting other formats or inferred charges.
    if guess_pybel_type(args.in_file) == 'xyz' and args.charge is not None:
        try:
            import torch
            from ase.io import iread, write
            from pyar.AIMNet2.calculators.aimnet2ase import AIMNet2Calculator
        except ImportError as exc:
            error = optional_dependency_error(
                getattr(exc, "name", None) or "torch",
                feature="AIMNet2 optimizer", extra="aimnet2",
            )
            raise SystemExit(str(error)) from None

        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model = torch.jit.load(args.model, map_location=device)
        with open(args.out_file, 'w') as output:
            for atoms in iread(args.in_file, index=':'):
                atoms.calc = AIMNet2Calculator(model, charge=args.charge)
                optimize(atoms, prec=args.fmax, steps=2000, traj=args.traj)
                write(output, atoms, format=guess_pybel_type(args.out_file))
        return

    try:
        import torch
        from openbabel import pybel as openbabel_pybel
        from pyar.AIMNet2.calculators.aimnet2ase import AIMNet2Calculator
    except ImportError as exc:
        module_name = getattr(exc, "name", None) or "openbabel"
        if module_name == "openbabel":
            raise SystemExit(
                'openbabel is required for pyar-aimnet2-ase-opt. '
                'Install with `mamba install -c conda-forge openbabel`.'
            ) from exc
        error = optional_dependency_error(
            module_name,
            feature="AIMNet2 optimizer",
            extra="aimnet2",
        )
        raise SystemExit(str(error)) from None
    pybel = openbabel_pybel

    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print('Loading AIMNet2 model from file', args.model)
    model = torch.jit.load(args.model, map_location=device)
    calc = AIMNet2Calculator(model)
    in_format = guess_pybel_type(args.in_file)
    out_format = guess_pybel_type(args.out_file)
    with open(args.out_file, 'w') as f:
        for mol in pybel.readfile(in_format, args.in_file):
            atoms = pybel2atoms(mol)
            charge = args.charge if args.charge is not None else guess_charge(mol)
            calc.do_reset()
            calc.set_charge(charge)
            atoms.set_calculator(calc)
            optimize(atoms, prec=args.fmax, steps=2000, traj=args.traj)
            update_mol(mol, atoms, align=False)
            f.write(mol.write(out_format))
            f.flush()


if __name__ == '__main__':
    main()
