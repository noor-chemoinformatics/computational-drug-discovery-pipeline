import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

ligands = {}
with open("data/ligands/gpr35/denovo/denovo_gpr35_ligands.smi") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2: ligands[parts[1]] = parts[0]

# Your new docking scores
docking_scores = {
    "DeNovo_1": -7.316, "DeNovo_2": -6.721, "DeNovo_3": -9.567,
    "DeNovo_4": -10.89, "DeNovo_5": -7.289, "DeNovo_6": -8.486,
    "DeNovo_7": -8.308, "DeNovo_8": -9.36, "DeNovo_9": -10.06,
    "DeNovo_10": -8.959, "DeNovo_11": -8.945, "DeNovo_12": -11.16,
    "DeNovo_13": -9.339, "DeNovo_14": -9.396
}

data = []
for name, score in docking_scores.items():
    mol = Chem.MolFromSmiles(ligands[name])
    mw, logp, tpsa = Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    
    # 1. Calculate Ligand Efficiency (LE)
    heavy_atoms = mol.GetNumHeavyAtoms()
    le = abs(score) / heavy_atoms
    
    # 2. Calculate MPO Score
    penalty = 1.0
    if tpsa > 140: penalty -= 0.2
    norm_dock = max(0.0, min(1.0, (score + 5) / -7))
    mpo_score = (0.7 * norm_dock) + (0.3 * penalty)
    
    data.append({
        "Ligand": name, "Docking (kcal/mol)": score, "LE": round(le, 2),
        "MW": round(mw, 1), "LogP": round(logp, 2), "MPO Score": round(mpo_score, 3)
    })

df = pd.DataFrame(data).sort_values(by="MPO Score", ascending=False)
df.to_csv("results/gpr35/denovo_mpo_ranked.csv", index=False)

print("\n--- Top De Novo Hits (Ranked by MPO) ---")
print(df.to_string(index=False))
