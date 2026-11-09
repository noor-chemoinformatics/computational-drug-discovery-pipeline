import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

docking_scores = {
    "Kynurenic_acid": -8.405, "Zaprinast": -6.91, "Lodoxamide": -8.773,
    "Salicylic_acid": -6.121, "Benzanilide": -7.769, "Phthalimide": -6.926,
    "Chlorobenzoic_acid": -6.176, "Fluorobenzoic_acid": -6.219,
    "Benzamide": -6.139, "Nitrobenzoic_acid": -6.497
}

ligands = {}
with open("data/ligands/gpr35/known_gpr35_ligands.smi") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2: ligands[parts[1]] = parts[0]

data = []
for name, score in docking_scores.items():
    mol = Chem.MolFromSmiles(ligands[name])
    mw, logp, tpsa = Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    
    penalty = 1.0
    if mw > 500 or logp > 5 or hbd > 5 or hba > 10: penalty = 0.5
    if tpsa > 140: penalty -= 0.2
    norm_dock = max(0.0, min(1.0, (score + 5) / -5))
    mpo_score = (0.7 * norm_dock) + (0.3 * penalty)
    
    data.append({
        "Ligand": name, "Docking (kcal/mol)": score, "MW": round(mw, 1),
        "LogP": round(logp, 2), "TPSA": round(tpsa, 1), "HBD": hbd, "HBA": hba, "MPO Score": round(mpo_score, 3)
    })

df = pd.DataFrame(data).sort_values(by="MPO Score", ascending=False)
df.to_csv("results/gpr35/mpo_ranked_hits.csv", index=False)
print(df.to_string(index=False))
