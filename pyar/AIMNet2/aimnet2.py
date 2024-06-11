from rdkit import Chem
import torch
torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False



# device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
device = torch.device('cpu')
print(device)
in_path = '/scratch/20cy91r19/bitbucket/pyatomgen/pyar/AIMNet2/data/examples.sdf'

aimnet2 = torch.jit.load('/scratch/20cy91r19/bitbucket/pyatomgen/pyar/AIMNet2/models/aimnet2_wb97m-d3_0.jpt', map_location=device)
# aimnet2_0 = torch.jit.load('/scratch/20cy91r19/bitbucket/afirxn/afirxn/AIMNet2/models/aimnet2_wb97m-d3_0.jpt', map_location=device)

def sdf2aimnet_input(sdf: str, device=torch.device('cpu')) -> dict:
    """Converts sdf to aimnet input, assuming the sdf has only 1 conformer."""
    mol = next(Chem.SDMolSupplier(sdf, removeHs=False))
    conf = mol.GetConformer()
    coord = torch.tensor(conf.GetPositions(), device=device).unsqueeze(0)
    numbers = torch.tensor([atom.GetAtomicNum() for atom in mol.GetAtoms()], device=device).unsqueeze(0)
    charge = torch.tensor([Chem.GetFormalCharge(mol)], device=device, dtype=torch.float)
    return dict(coord=coord, numbers=numbers, charge=charge)


dct = sdf2aimnet_input(in_path, device=device)
dct['coord'].requires_grad_(True)
aimnet2_out = aimnet2(dct)
print('aimnet2 energy: ', aimnet2_out['energy'])  # there is no gradient for energy

for key, val in aimnet2_out.items():
    print(key, val.shape)

print('aimnet2 energy: ', aimnet2_out['energy'])

"""## aimnet2 single model"""

# aimnet2_0out = aimnet2_0(dct)
# print('aimnet2_0 energy: ', aimnet2_0out['energy'])  # there is gradient function for energy

# forces = -torch.autograd.grad(aimnet2_0out['energy'], dct['coord'], create_graph=True)[0]
# print(forces)
# print(forces.shape)

# using pytorch function for calculating the hessian
def func(coord, numbers=dct['numbers'], charge=dct['charge']):
    dct = dict(coord=coord, numbers=numbers, charge=charge)
    return aimnet2(dct)['energy']

hess2 = torch.autograd.functional.hessian(func, dct['coord'])
print(hess2)
print(hess2.shape)

hess3 = hess2.view(20, 3, 20, 3)
print(hess2)