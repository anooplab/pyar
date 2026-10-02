# MPCONF196GEN Benchmark
*Molecular Conformer Generation Benchmark*

<img src='img/macrocycle_fail.jpg' width=350>

A conformer generation benchmark created from the MPCONF196 benchmark energy dataset [1]. This evaluates a method's ability to generate conformers of macrocycles, a difficult task for molecular conformer generation.

We generated reference conformers with CREST iMTD-GC [2] using GFN2-xTB [3] and computed conformer energies using Egret-1 [4] neural network potential.

## Usage
- `SMILES.txt` contains the number of rotateable bonds, canonical SMILES strings and identifier for each molecule.
-  Each `molecules/{molecule_id}/` directory contains a `crest_conformers.xyz` file with conformers and energies for each molecule.

## Creating the Benchmark

All packages necessary for making this benchmark can be installed by creating a conda environment from `env.yml`.

We ran CREST on each `MPCONF196_data/{molecule_id}.xyz` file, outputs and conformers can be found in `molecules/{molecule_id}/crest.out` and `molecules/{molecule_id}/crest_conformers.xyz` respectively. The Egret-1 checkpoint and scripts we used to calculate conformer energies and create `SMILES.txt` can be found in `tools/`. We generated SMILES for each molecule from the `MPCONF196_data/{molecule_id}.xyz` file with [rowansci.com](https://rowansci.com) [5], except for `MPCONF196_data/YIVNOG_I.xyz` which we made by hand.

### Settings

- **Conformer search:** CREST 3.0.2, default iMTD-GC workflow at the GFN2-xTB level (`crest MPCONF196_data/{molecule_id}.xyz -T 16`). All other settings were left at their defaults (300 K Fermi temperature, 6 kcal/mol final energy window).
- **Conformer energies:** Single-point energies on the CREST ensemble computed with the Egret-1 neural network potential (`tools/get_egret_energies.py`, `mace_off` with `tools/EGRET_1.model`, `float64`), reported in kcal/mol.

## License

MPCONF196GEN Benchmark © 2026 Rowan Scientific Corporation. This benchmark (data and code) is released under the [Creative Commons Attribution 4.0 International License (CC-BY 4.0)](LICENSE).

## Citation

If you use the MPCONF196GEN benchmark in your research, please cite:

```bibtex
@software{mpconf196gen,
  title = {MPCONF196GEN: A Molecular Conformer Generation Benchmark},
  year = {2026},
  url = {https://github.com/rowansci/mpconf196-gen-benchmark}
}
```

## References

[1] Řezáč, J.; Bím, D.; Gutten, O.; Rulíšek, L. *Toward Accurate Conformational Energies of Smaller Peptides and Medium-Sized Macrocycles: MPCONF196 Benchmark Energy Data Set.* J. Chem. Theory Comput. 2018, 14 (3), 1254–1266. https://doi.org/10.1021/acs.jctc.7b01074

[2] Pracht, P.; Grimme, S.; Bannwarth, C.; Bohle, F.; Ehlert, S.; Feldmann, G.; Gorges, J.; Müller, M.; Neudecker, T.; Plett, C.; Spicher, S.; Steinbach, P.; Wesołowski, P. A.; Zeller, F. *CREST—A Program for the Exploration of Low-Energy Molecular Chemical Space.* J. Chem. Phys. 2024, 160 (11), 114110. https://doi.org/10.1063/5.0197592

[3] Bannwarth, C.; Ehlert, S.; Grimme, S. *GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions.* J. Chem. Theory Comput. 2019, 15 (3), 1652–1671. https://doi.org/10.1021/acs.jctc.8b01176

[4] Mann, E. L.; Wagen, C. C.; Vandezande, J. E.; Wagen, A. M.; Schneider, S. C. *Egret-1: Pretrained Neural Network Potentials for Efficient and Accurate Bioorganic Simulation.* arXiv:2504.20955, 2025. https://doi.org/10.48550/arXiv.2504.20955

[5] *Rowan Scientific.* https://www.rowansci.com

