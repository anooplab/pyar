#!/usr/bin/bash
python ./calculators/aimnet2_ase_opt.py  ./models/aimnet2_wb97m-d3_ens.jpt  --charge 0 --traj test-1.traj test-1.xyz test-1_min.xyz

# Convert test-1.traj to output.xyz with ASE
$ase convert test-1.traj output.xyz