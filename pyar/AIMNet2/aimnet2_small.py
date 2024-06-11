import torch


model = torch.jit.load('/scratch/20cy91r19/bitbucket/afirxn/afirxn/AIMNet2/models/aimnet2_wb97m-d3_ens.jpt')

# Water molecule: H-O-H
coord = torch.as_tensor([   (0.7494, 0.0000,  0.4424),
                            (0.0000, 0.0000, -0.1654), 
                            (-0.7494, 0.0000,  0.4424)]).to(torch.float).unsqueeze(0)
numbers = torch.as_tensor([1, 8, 1]).to(torch.long).unsqueeze(0)
charge = torch.as_tensor([0.0]).to(torch.float)

d = dict(coord=coord, numbers=numbers, charge=charge)

out = model(d)

ret = dict(energy=out['energy'].item())
print(ret)