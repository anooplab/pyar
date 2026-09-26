# ORCA bond-scan endpoint validation

Run on 2026-09-25 with ORCA 6.1.1, BP86/def2-SVP, four scan points, one
orientation, and one ORCA process per job. The five association cases use
CH3+CH3 (C-C), CH3+NH2 (C-N), CH3+OH (C-O), OH+H (O-H), and CH3+Cl (C-Cl).
Both fragments were neutral doublets; the merged calculations were singlets.
Each endpoint was the selected atoms' covalent-radius sum multiplied by the
factor in the table. All 15 jobs completed the scan and unconstrained
relaxation. Temporary outputs were kept under `/tmp/pyar-scan-benchmark-*`.

| Target | Radius sum (A) | Factor | Endpoint (A) | Relaxed distance (A) | New identity | Contact after relax | Max scan point | Relative max (kcal/mol) | Max internal |
|---|---:|---:|---:|---:|---|---|---:|---:|---|
| C-C | 1.500 | 1.0 | 1.500 | 1.531 | yes | yes | 1/4 | 0.000 | no |
| C-C | 1.500 | 0.9 | 1.350 | 1.531 | yes | yes | 1/4 | 0.000 | no |
| C-C | 1.500 | 0.8 | 1.200 | 1.531 | yes | yes | 4/4 | 10.715 | no |
| C-N | 1.460 | 1.0 | 1.460 | 1.462 | yes | yes | 1/4 | 0.000 | no |
| C-N | 1.460 | 0.9 | 1.314 | 1.461 | yes | yes | 1/4 | 0.000 | no |
| C-N | 1.460 | 0.8 | 1.168 | 1.462 | yes | yes | 4/4 | 0.107 | no |
| C-O | 1.380 | 1.0 | 1.380 | 1.414 | yes | yes | 1/4 | 0.000 | no |
| C-O | 1.380 | 0.9 | 1.242 | 1.414 | yes | yes | 1/4 | 0.000 | no |
| C-O | 1.380 | 0.8 | 1.104 | 1.414 | yes | yes | 4/4 | 18.155 | no |
| O-H | 0.950 | 1.0 | 0.950 | 0.975 | yes | yes | 1/4 | 0.000 | no |
| O-H | 0.950 | 0.9 | 0.855 | 0.975 | yes | yes | 1/4 | 0.000 | no |
| O-H | 0.950 | 0.8 | 0.760 | 0.975 | yes | yes | 1/4 | 0.000 | no |
| C-Cl | 1.740 | 1.0 | 1.740 | 1.797 | yes | yes | 1/4 | 0.000 | no |
| C-Cl | 1.740 | 0.9 | 1.566 | 1.797 | yes | yes | 1/4 | 0.000 | no |
| C-Cl | 1.740 | 0.8 | 1.392 | 1.798 | yes | yes | 1/4 | 0.000 | no |

All five pairs relaxed to the same product identity across the three factors;
the relaxed bond distances also agree within 0.001 A per pair. The scan maximum
was at an endpoint in every case, so none provides an internal transition-state
candidate. The 0.8 endpoint raises the sampled energy above the starting point
for C-C, C-N, and C-O, while O-H and C-Cl remain highest at the initial point.
These coarse four-point scans are endpoint sensitivity checks, not a calibrated
chemical benchmark and not evidence that 0.8 is universally optimal. Retain
the explicit `--scan-end` override for system-specific work.

A separate public-CLI HCN + HCN run (`tests/data/neb/hcn.xyz`, selected C-N
atoms, default 0.8 endpoint of 1.168 A, four points) also completed scan and
free relaxation. Its relaxed fragments separated (target distance 3.178 A,
contact absent) and canonical identity was unchanged. The scan maximum was the
compressed endpoint, not an internal maximum. This confirms the workflow
reports a failed-to-form product rather than treating the constrained endpoint
as product evidence.

The HCN acceptance command was:

```bash
python -m pyar.scripts.scan_bond tests/data/neb/hcn.xyz tests/data/neb/hcn.xyz \
  --atoms 1 2 -N 1 --software orca --method BP86 --basis def2-SVP \
  --nprocs 1 --scf-cycles 1000 --opt-cycles 100 --scan-points 4 \
  --output /tmp/pyar-scan-benchmark-hcn-default-np1
```

Two processes were not usable in this environment: ORCA MPI startup failed
because no network interface was available. The validated run used one process.
