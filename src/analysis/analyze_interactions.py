import prolif as plf
from rdkit import Chem
import warnings, subprocess
warnings.filterwarnings("ignore")

print("Loading receptor...")
prot_mol = Chem.MolFromPDBFile("data/structures/gpr35/gpr35_receptor_clean.pdb", removeHs=False, sanitize=False)
prot_mol = plf.Molecule(prot_mol)

print("Loading docked ligand...")
subprocess.run(["obabel", "results/gpr35/lodoxamide_docked.pdbqt", "-O", "docked_clean.pdb"], check=True)
lig_mol = Chem.MolFromPDBFile("docked_clean.pdb", removeHs=False, sanitize=False)
lig_mol = plf.Molecule(lig_mol)

print("Calculating interactions...")
fp = plf.Fingerprint()
fp.run_from_iterable([lig_mol], prot_mol)
df = fp.to_dataframe()
interactions = df.iloc[0]
active_interactions = interactions[interactions > 0]

print("\n--- Key Interactions for GPR35-Lodoxamide ---")
for interaction, count in active_interactions.items():
    print(f"{interaction[2]} between {interaction[0]} and {interaction[1]}")
