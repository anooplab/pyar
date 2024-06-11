#!/bin/bash
#SBATCH -J gpu-aimnet2  # name of the job
#SBATCH -p gpu     # name of the partition
#SBATCH -n 8          # no of processes or tasks
#SBATCH -N 1

#SBATCH --gres=gpu:1               # request gpu card: it should be either 1 or 2

#SBATCH --cpus-per-task=4          # no of threads per process or task
#SBATCH -t 00:30:00   # walltime in HH:MM:SS, Max value 72:00:00

module try-add python/3.7.4 #list of modules you want to use, for example
module load pymlgen/2.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#pymlgen-cli -s hcn.xyz hcn.xyz --software aimnet_2 -ss 20  -N 16 -c 0 0 -m 1 1  
python /scratch/20cy91r19/bitbucket/pyatomgen/pymlgen/AIMNet2/calculators/aimnet2_ase_opt.py  /scratch/20cy91r19/bitbucket/pyatomgen/pymlgen/AIMNet2/models/aimnet2_wb97m-d3_ens.jpt  --charge 0 --traj test-1.traj /scratch/20cy91r19/bitbucket/pyatomgen/pymlgen/AIMNet2/test-1.xyz /scratch/20cy91r19/bitbucket/pyatomgen/pymlgen/AIMNet2/test-1_gpu_min.xyz
