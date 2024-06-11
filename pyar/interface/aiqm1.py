import random
import pyar.tabu
import ase
import torch
import numpy as np
import pandas as pd
from ase.optimize import BFGS
import os
os.system('pip install ase')
import matplotlib.pyplot as plt
from tqdm import tqdm
import sklearn
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import QED

# Tabu search parameters
max_iterations = 1000
tabu_list_size = 50
aspiration_criteria = 0.8 

# ANI-1x potential 

import torchani

model = torch.nn.Sequential(torchani.models.ANI())
# Reward functions
from rdkit.Chem import Descriptors, Lipinski, MolFromSmiles
from rdkit.Chem.Descriptors import MoleculeDescriptors

def druglikeness(mol):
    descriptors = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors.descList])
    descs = descriptors.CalcDescriptors(mol)
    logP = descs[6]
    SA = descs[11]
    weight = descs[12] 
    rings = descs[34]
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    
    score = 0
    if -0.4 <= logP <= 5.6:
        score += 1
    if 75 <= SA <= 130:
        score += 1 
    if 160 <= weight <= 480:
        score += 1
    if rings >= 1:
        score += 1
    if rot_bonds <= 10:
        score += 1
        
    return score / 5

# Example usage:
smiles = 'CCOc1ccc2nc(S(N)(=O)=O)sc2c1'
mol = MolFromSmiles(smiles)
print(druglikeness(mol))

def synthetic_availability(mol):
    # Define chembl_dataset
    chembl_dataset = pd.read_csv('chembl.csv')
    
    # Train synthetic accessibility classifier on ChEMBL dataset
    chembl_mols = [Chem.MolFromSmiles(s) for s in chembl_dataset['SMILES']]
    chembl_sas = chembl_dataset['SyntheticAccessibility']
    X = np.array([molecule_to_fingerprint(m) for m in chembl_mols])
    y = np.array(chembl_sas.replace('Isomer', False).astype(float))
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    
    model = MLPRegressor(hidden_layer_sizes=(100,100,100), max_iter=500)
    model.fit(X_train, y_train)

    # Predict for input molecule
    fp = molecule_to_fingerprint(mol)
    sa_score = model.predict([fp])[0]
    
    # Access X_test and y_test
    model.score(X_test, y_test)
    
    return sa_score

# Generate initial population
population = []
for i in range(100):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def generate_random_molecule():
        mol = Chem.MolFromSmiles('C')
        mol = AllChem.AddHs(mol)
        builder = AllChem.EditableMol(mol)
        for i in range(5):
            builder.AddAtom(Chem.Atom('C'))
            builder.AddBond(i, i+1, Chem.BondType.SINGLE)
        mol = builder.GetMol()
        AllChem.EmbedMolecule(mol)
        return mol
    population.append(mol)

for i in tqdm(range(max_iterations)):
    # Get current best molecule
    current_best = max(population, key=lambda x: reward(x))
    
    # Generate neighboring molecules
    neighbors = []
    for j in range(10):
        def mutate_molecule(mol):
            """Mutate a molecule by randomly modifying one of its atoms or bonds."""
            # Copy the molecule to avoid modifying the original
            mol = Chem.Mol(mol)
            
            # Randomly choose an atom or bond to modify
            atom_or_bond = random.choice(['atom', 'bond'])
            if atom_or_bond == 'atom':
                atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
                new_atom = Chem.Atom(random.choice(['C', 'N', 'O', 'S']))
                mol.ReplaceAtom(atom_idx, new_atom)
            else:
                bond_idx = random.randint(0, mol.GetNumBonds() - 1)
                new_bond_type = random.choice([Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE])
                mol.GetBondWithIdx(bond_idx).SetBondType(new_bond_type)
            
            return mol
        # Generate neighboring molecules
        neighbors = []
        for j in range(10):
            neighbor = mutate_molecule(current_best)
            neighbors.append(neighbor)
            
        # Get best neighbor
        best_neighbor = max(neighbors, key=lambda x: reward(x))

        # Add to population, respecting tabu list
    
    # Get best neighbor
    best_neighbor = max(neighbors, key=lambda x: reward(x))
    
    # Add to population, respecting tabu list
    from rdkit.Chem import QED
    import torch
    import random
    from tqdm import tqdm
    from rdkit.Chem import Descriptors, Lipinski, MolFromSmiles
    from rdkit.Chem.Descriptors import MoleculeDescriptors
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Tabu search parameters
    max_iterations = 1000
    tabu_list_size = 50
    aspiration_criteria = 0.8 

    # ANI-1x potential 
    model = torch.nn.Sequential()

    # Reward functions
    def druglikeness(mol):
        pass

    def synthetic_availability(mol):
        pass

    def reward(mol):
        qed = QED.qed(mol)
        logp = Descriptors.MolLogP(mol)
        sa = synthetic_availability(mol)
        dl = druglikeness(mol)
        return qed - logp - sa - dl

    def generate_random_molecule():
        # Generate random SMILES string
        chars = 'BCNOFPSIHKLMgtnospbclsi'
        smi = ''.join([random.choice(chars) for _ in range(100)])
        
        # Convert to RDKit molecule
        mol = Chem.MolFromSmiles(smi)
        
        # Perform sanity check
        try:
            Chem.SanitizeMol(mol)
            return mol
        except ValueError:
            return None

    def add_to_tabu_list(mol):
        """Add a molecule to the tabu list."""
        tabu_list[mol] = len(tabu_list)

    def find_most_tabu(population):
        """Find the molecule in the population that has been in the tabu list the longest."""
        oldest_mol = None
        oldest_age = -1
        for mol in population:
            age = tabu_list.get(mol, 0)
            if age > oldest_age:
                oldest_mol = mol
                oldest_age = age
        return oldest_mol

    # Generate initial population
    population = []
    for i in range(100):
        mol = generate_random_molecule()
        if mol is not None:
            population.append(mol)

    # Initialize tabu list
    tabu_list = {}

    for i in tqdm(range(max_iterations)):
        # Get current best molecule
        current_best = max(population, key=lambda x: reward(x))
        
        # Generate neighboring molecules
        neighbors = []
        for j in range(10):
            def mutate_molecule(mol):
                """Mutate a molecule by randomly modifying one of its atoms or bonds."""
                # Copy the molecule to avoid modifying the original
                mol = Chem.Mol(mol)
                
                # Randomly choose an atom or bond to modify
                atom_or_bond = random.choice(['atom', 'bond'])
                if atom_or_bond == 'atom':
                    atom_idx = random.randint(0, mol.GetNumAtoms() - 1)
                    new_atom = Chem.Atom(random.choice(['C', 'N', 'O', 'S']))
                    mol.ReplaceAtom(atom_idx, new_atom)
                else:
                    bond_idx = random.randint(0, mol.GetNumBonds() - 1)
                    new_bond_type = random.choice([Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE])
                    mol.GetBondWithIdx(bond_idx).SetBondType(new_bond_type)
                
                return mol
            neighbors.append(mutate_molecule(current_best))
        
        # Get best neighbor
        best_neighbor = max(neighbors, key=lambda x: reward(x))
        
        # Add to population, respecting tabu list
        if best_neighbor not in tabu_list or reward(best_neighbor) > reward(find_most_tabu(population)):
            population.append(best_neighbor)
        
        # Update tabu list
        add_to_tabu_list(current_best)
        if len(tabu_list) > tabu_list_size:
            tabu_list.pop(0)

    # Final results
    print('Best molecule found:')
    print(population[0])
    print('Reward:', reward(population[0]))
    population.append(best_neighbor)
    
    # Update tabu list
    add_to_tabu_list(current_best)
    if len(tabu_list) > tabu_list_size:
        tabu_list.pop(0)

# Final results
print('Best molecule found:')
print(population[0])
print('Reward:', reward(population[0]))

import random
from rdkit import Chem

def generate_random_molecule():
    # Define the characters to use for generating the SMILES string
    chars = 'BCNOFPSIHKLMgtnospbclfiBr'
    
    # Generate random SMILES string
    smi = ''.join([random.choice(chars) for _ in range(100)])
    
    # Convert to RDKit molecule
    mol = Chem.MolFromSmiles(smi)
    
    # Perform sanity check
    try:
        Chem.SanitizeMol(mol)
        return mol
    except ValueError:
        return generate_random_molecule() # Retry
            
def mutate_molecule(mol):
    # Randomly perform one of the following mutations:
    # - Remove/add atom
    # - Remove/add bond
    # - Alter functional group
    
    mutation = random.choice(['remove_atom', 'add_atom', 'remove_bond', 
                              'add_bond', 'alter_functional_group'])
    
    new_mol = Chem.RWMol(mol) # Mutable molecule
    
    if mutation == 'remove_atom':
        new_mol.RemoveAtom(random.randrange(0, new_mol.GetNumAtoms()))
    elif mutation == 'add_atom':
        import random
        from rdkit import Chem

        atom_types = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I']

        def generate_random_molecule():
            # Define the characters to use for generating the SMILES string
            chars = ''.join(atom_types)
            
            # Generate random SMILES string
            smi = ''.join([random.choice(chars) for _ in range(100)])
            
            # Convert to RDKit molecule
            mol = Chem.MolFromSmiles(smi)
            
            # Perform sanity check
            try:
                Chem.SanitizeMol(mol)
                return mol
            except ValueError:
                return generate_random_molecule() # Retry
    elif mutation == 'remove_bond':
        bond_idx = random.randrange(0, new_mol.GetNumBonds())
        new_mol.RemoveBond(bond_idx) 
    elif mutation == 'add_bond':
        atom1 = random.randrange(0, new_mol.GetNumAtoms())
        atom2 = random.randrange(0, new_mol.GetNumAtoms())
        if atom1 != atom2:
            new_mol.AddBond(atom1, atom2, Chem.BondType.SINGLE)
    elif mutation == 'alter_functional_group':
        # Replace functional group with another viable option
        RGROUPS = [
            Chem.MolFromSmarts('[#6X3:1](=[OX1])[OX2H0:2]'), 
            Chem.MolFromSmarts('[#6X3:1]-[#8X2H1:2]'),
        ]
        
        for i, rg in enumerate(RGROUPS):
            for match in new_mol.GetSubstructMatches(rg):
                new_rg = Chem.MolFromSmarts(RGROUPS[(i+1)%len(RGROUPS)].GetSmarts()) 
                new_mol.ReplaceSubstructs(match, new_rg, True)
                break
        
    new_mol = new_mol.GetMol()
    
    try:
        Chem.SanitizeMol(new_mol)
        return new_mol 
    except ValueError:
        return mutate_molecule(mol) # Retry with different mutation
            
def reward(mol):
    # Reward function balances drug-likeness and synthetic accessibility
    dl_score = druglikeness(mol)
    sa_score = synthetic_availability(mol)
    return 0.5*dl_score + 0.5*sa_score

def molecule_to_fingerprint(mol):
    # Convert molecule to fingerprint using RDKit
    fingerprinter = Chem.RDKFingerprint(maxPath=5, fpSize=2048)
    return fingerprinter.GetFingerprint(mol)
    
def find_most_tabu(population):
    # Return molecule from population that is most tabu
    tabu_values = {m: 0 for m in population}
    for m in population:
        for t in tabu_list:
            if molecule_similarity(m, t) > 0.4:
                tabu_values[m] += 1
    return max(tabu_values, key=tabu_values.get) 

def add_to_tabu_list(mol):
    # Add molecule to tabu list
    tabu_list.append(mol)
    
def molecule_similarity(a, b):
    # Calculate similarity between two molecules
    return Chem.rdMolAlign.AlignMol(a, b).AlignMol(a, b).GetFitScore()
    
if __name__ == '__main__':
    import random
    from rdkit import Chem
    from tqdm import tqdm

    def main():
        # Generate initial population
        population = []
        for i in range(100):
            mol = generate_random_molecule()
            if mol is not None:
                population.append(mol)

        # Initialize tabu list
        tabu_list = []

        max_iterations = 100
        tabu_list_size = 10

        for i in tqdm(range(max_iterations)):
            # Get current best molecule
            current_best = max(population, key=lambda x: reward(x))
            
            # Generate neighboring molecules
            neighbors = []
            for j in range(10):
                neighbors.append(mutate_molecule(current_best))
            
            # Get best neighbor
            best_neighbor = max(neighbors, key=lambda x: reward(x))
            
            # Add to population, respecting tabu list
            if best_neighbor not in tabu_list or reward(best_neighbor) > reward(find_most_tabu(population)):
                population.append(best_neighbor)
            
            # Update tabu list
            add_to_tabu_list(current_best)
            if len(tabu_list) > tabu_list_size:
                tabu_list.pop(0)

        # Final results
        print('Best molecule found:')
        print(population[0])
        print('Reward:', reward(population[0]))

    def generate_random_molecule():
        # Define the characters to use for generating the SMILES string
        elements = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I']
        functional_groups = ['C=O', 'C#N', 'C=C', 'C#C', 'N-H', 'N=C', 'N#C', 'N#N', 'O-H', 'O=C', 'O-C', 'S-H', 'S-C', 'F', 'Cl', 'Br', 'I']
        chars = elements + functional_groups + ['(', ')', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '-', '[', ']']

        # Generate random SMILES string
        smiles = ''.join(random.choices(chars, k=random.randint(5, 20)))

        # Convert SMILES string to molecule
        mol = Chem.MolFromSmiles(smiles)

        return mol

    def mutate_molecule(mol):
        # Randomly perform one of the following mutations:
        # - Remove/add atom
        # - Remove/add bond
        # - Alter functional group
        
        # Retry with different mutation
        return generate_random_molecule()

    def reward(mol):
        # Reward function balances drug-likeness and synthetic accessibility
        dl_score = druglikeness(mol)
        sa_score = synthetic_availability(mol)
        return 0.5*dl_score + 0.5*sa_score

    def molecule_to_fingerprint(mol):
        # Convert molecule to fingerprint using RDKit
        fingerprinter = Chem.RDKFingerprint(maxPath=5, fpSize=2048)
        return fingerprinter.GetFingerprint(mol)
        
    def find_most_tabu(population):
        # Return molecule from population that is most tabu
        tabu_values = {m: 0 for m in population}
        for m in population:
            for t in tabu_list:
                if molecule_similarity(m, t) > 0.4:
                    tabu_values[m] += 1
        return max(tabu_values, key=tabu_values.get) 

    def add_to_tabu_list(mol):
        # Add molecule to tabu list
        tabu_list.append(mol)
        
    def molecule_similarity(a, b):
        # Calculate similarity between two molecules
        return Chem.rdMolAlign.AlignMol(a, b).AlignMol(a, b).GetFitScore()

    if __name__ == '__main__':
        main()

