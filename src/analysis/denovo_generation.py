import os, random
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors

os.makedirs("data/ligands/gpr35/denovo", exist_ok=True)
seed_smiles = []
with open("data/ligands/gpr35/known_gpr35_ligands.smi") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2: seed_smiles.append(parts[0])

all_fragments = set()
for smiles in seed_smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol: all_fragments.update(BRICS.BRICSDecompose(mol))

novel_mols = []
for i in range(5000):
    chosen_frags = random.sample(list(all_fragments), random.choice([2, 3]))
    
    for mol in BRICS.BRICSBuild([Chem.MolFromSmiles(f) for f in chosen_frags]):
        if mol is None: continue
        
        # FIX: Sanitize the molecule to fix valences before calculating descriptors
        try:
            Chem.SanitizeMol(mol)
        except:
            continue # Skip if the molecule is chemically invalid
            
        if 150 < Descriptors.MolWt(mol) < 550 and Descriptors.MolLogP(mol) < 6:
            smi = Chem.MolToSmiles(mol)
            if smi not in [m[0] for m in novel_mols]:
                novel_mols.append((smi, f"DeNovo_{len(novel_mols)+1}"))
                
        if len(novel_mols) >= 20: break
    if len(novel_mols) >= 20: break

with open("data/ligands/gpr35/denovo/denovo_gpr35_ligands.smi", "w") as f:
    for smi, name in novel_mols: f.write(f"{smi}\t{name}\n")
print(f"Successfully generated {len(novel_mols)} novel molecules.")
