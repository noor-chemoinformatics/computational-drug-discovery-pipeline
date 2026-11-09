import os, subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

ligand_dir = "data/ligands/gpr35/denovo"
receptor = "data/structures/gpr35/gpr35_receptor.pdbqt"
results_dir = "results/gpr35/denovo"
os.makedirs(results_dir, exist_ok=True)

ligands = []
with open(f"{ligand_dir}/denovo_gpr35_ligands.smi") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2: ligands.append((parts[0], parts[1]))

for smiles, name in ligands:
    print(f"--- Processing {name} ---")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        print(f"  ERROR: RDKit could not parse SMILES: {smiles}")
        continue
        
    mol = Chem.AddHs(mol)
    # Try to generate 3D coordinates
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        print(f"  ERROR: Could not embed 3D coordinates for {name}")
        continue
        
    AllChem.MMFFOptimizeMolecule(mol)
    
    pdb_file = f"{ligand_dir}/{name}.pdb"
    Chem.MolToPDBFile(mol, pdb_file)
    
    pdbqt_file = f"{ligand_dir}/{name}.pdbqt"
    obabel_res = subprocess.run(["obabel", pdb_file, "-O", pdbqt_file], capture_output=True, text=True)
    
    out_file = f"{results_dir}/{name}_docked.pdbqt"
    vina_cmd = ["vina", "--receptor", receptor, "--ligand", pdbqt_file,
        "--center_x", "145.081", "--center_y", "139.934", "--center_z", "162.850",
        "--size_x", "18.263", "--size_y", "14.263", "--size_z", "17.820",
        "--out", out_file, "--exhaustiveness", "8"]
    
    result = subprocess.run(vina_cmd, capture_output=True, text=True)
    
    best_score = None
    for line in result.stdout.split("\n"):
        if line.strip().startswith("1 "):
            best_score = float(line.split()[1])
            break
            
    if best_score: 
        print(f"  >>> SUCCESS: {best_score} kcal/mol\n")
    else:
        print(f"  FAILED. Vina Error Output:\n{result.stderr}")
        print(f"  Vina Standard Output:\n{result.stdout[:200]}...")
        print()
