import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# Read the generated De Novo SMILES
ligands = {}
with open("data/ligands/gpr35/denovo/denovo_gpr35_ligands.smi") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2: ligands[parts[1]] = parts[0]

# Our top scores from the docking run
docking_scores = {
    "DeNovo_20": -10.57, "DeNovo_19": -10.56, "DeNovo_18": -10.01,
    "DeNovo_7": -9.047, "DeNovo_12": -8.851, "DeNovo_15": -8.642,
    "DeNovo_3": -8.2, "DeNovo_16": -8.324
}

data = []
for name, score in docking_scores.items():
    mol = Chem.MolFromSmiles(ligands[name])
    mw, logp, tpsa = Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    
    penalty = 1.0
    if mw > 500 or logp > 5 or hbd > 5 or hba > 10: penalty = 0.5
    if tpsa > 140: penalty -= 0.2
    
    # Normalize docking score (-12 is best, -5 is worst)
    norm_dock = max(0.0, min(1.0, (score + 5) / -7))
    mpo_score = (0.7 * norm_dock) + (0.3 * penalty)
    
    data.append({
        "Ligand": name, "SMILES": ligands[name], "Docking (kcal/mol)": score, 
        "MW": round(mw, 1), "LogP": round(logp, 2), "MPO Score": round(mpo_score, 3)
    })

df = pd.DataFrame(data).sort_values(by="MPO Score", ascending=False)
df.to_csv("results/gpr35/denovo_mpo_ranked.csv", index=False)
print(df.to_string(index=False))
