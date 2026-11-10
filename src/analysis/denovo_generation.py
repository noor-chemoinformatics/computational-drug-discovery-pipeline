import os, random
from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

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
unique_scaffolds = set()

print("Generating molecules with scaffold-level diversity...")

for i in range(5000):
    chosen_frags = random.sample(list(all_fragments), random.choice([2, 3]))
    
    for mol in BRICS.BRICSBuild([Chem.MolFromSmiles(f) for f in chosen_frags]):
        if mol is None: continue
        
        try:
            Chem.SanitizeMol(mol)
        except:
            continue 
            
        # 1. Calculate Lipinski properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        # Flag high-LogP and nitro-containing candidates as requested in your CV
        # (We will filter them out to keep only the best drug-like hits)
        has_nitro = "[N+](=O)[O-]" in Chem.MolToSmiles(mol)
        
        if 150 < mw < 550 and logp < 5 and not has_nitro:
            
            # 2. Extract the Murcko Scaffold for diversity
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            
            # 3. Only keep the molecule if its scaffold is NEW
            if scaffold not in unique_scaffolds:
                smi = Chem.MolToSmiles(mol)
                if smi not in [m[0] for m in novel_mols]:
                    unique_scaffolds.add(scaffold)
                    novel_mols.append((smi, f"DeNovo_{len(novel_mols)+1}"))
                
        if len(novel_mols) >= 20: break
    if len(novel_mols) >= 20: break

with open("data/ligands/gpr35/denovo/denovo_gpr35_ligands.smi", "w") as f:
    for smi, name in novel_mols: f.write(f"{smi}\t{name}\n")

print(f"Successfully generated {len(novel_mols)} molecules with {len(unique_scaffolds)} unique scaffolds.")
