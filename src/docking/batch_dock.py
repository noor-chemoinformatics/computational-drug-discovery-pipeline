import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

ligand_dir = "data/ligands/gpr35"
receptor = "data/structures/gpr35/gpr35_receptor.pdbqt"
results_dir = "results/gpr35"
os.makedirs(results_dir, exist_ok=True)

ligands = []
with open(f"{ligand_dir}/known_gpr35_ligands.smi") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            ligands.append((parts[0], parts[1]))

for smiles, name in ligands:
    print(f"--- Processing {name} ---")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: continue
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    pdb_file = f"{ligand_dir}/{name}.pdb"
    Chem.MolToPDBFile(mol, pdb_file)
    
    pdbqt_file = f"{ligand_dir}/{name}.pdbqt"
    subprocess.run(["obabel", pdb_file, "-O", pdbqt_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    out_file = f"{results_dir}/{name}_docked.pdbqt"
    vina_cmd = [
        "vina", "--receptor", receptor, "--ligand", pdbqt_file,
        "--center_x", "145.081", "--center_y", "139.934", "--center_z", "162.850",
        "--size_x", "18.263", "--size_y", "14.263", "--size_z", "17.820",
        "--out", out_file, "--exhaustiveness", "8"
    ]
    result = subprocess.run(vina_cmd, capture_output=True, text=True)
    best_score = None
    for line in result.stdout.split("\n"):
        if line.strip().startswith("1 "):
            best_score = float(line.split()[1])
            break
    if best_score:
        print(f">>> {name}: {best_score} kcal/mol\n")
