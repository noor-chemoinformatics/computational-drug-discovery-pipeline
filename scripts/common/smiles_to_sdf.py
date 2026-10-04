"""
Convert SMILES file to SDF with 3D coordinates.

Usage:
    python smiles_to_sdf.py <input.smi> <output.sdf>
"""
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def convert(input_smi, output_sdf):
    writer = Chem.SDWriter(output_sdf)
    count = 0
    with open(input_smi) as f:
        for line in f:
            if not line.strip():
                continue
            smi = line.strip().split()[0]
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                continue
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                pass
            writer.write(mol)
            count += 1
    writer.close()
    print(f"Converted {count} molecules to {output_sdf}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python smiles_to_sdf.py <input.smi> <output.sdf>")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
