"""
Batch prepare ligands from a SMILES or SDF file.

Usage:
    python batch_prepare_ligands.py <input.smi> <output_dir> [max_molecules]
"""
import sys
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_sdf(smiles_file, output_sdf, max_molecules=100):
    """Convert SMILES file to 3D SDF with hydrogens."""
    writer = Chem.SDWriter(output_sdf)
    count = 0
    skipped = 0

    with open(smiles_file, 'r') as f:
        for line in f:
            if count >= max_molecules:
                break

            parts = line.strip().split()
            if not parts:
                continue

            smi = parts[0]
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                skipped += 1
                continue

            # Add hydrogens and generate 3D coordinates
            mol = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                skipped += 1
                continue

            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                pass

            writer.write(mol)
            count += 1

    writer.close()
    print(f"Wrote {count} molecules to {output_sdf} (skipped {skipped})")
    return count


def sdf_to_pdbqt_batch(input_sdf, output_dir):
    """Convert SDF to individual PDBQT files using Meeko."""
    os.makedirs(output_dir, exist_ok=True)

    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    count = 0

    for i, mol in enumerate(supplier):
        if mol is None:
            continue

        # Write individual SDF
        sdf_file = os.path.join(output_dir, f"lig_{i:04d}.sdf")
        writer = Chem.SDWriter(sdf_file)
        writer.write(mol)
        writer.close()

        # Convert to PDBQT
        pdbqt_file = os.path.join(output_dir, f"lig_{i:04d}.pdbqt")
        cmd = ["mk_prepare_ligand.py", "-i", sdf_file, "-o", pdbqt_file]

        try:
            subprocess.run(cmd, check=True, capture_output=True)
            count += 1
        except subprocess.CalledProcessError as e:
            print(f"Failed to prepare ligand {i}: {e}")
        finally:
            # Clean up intermediate SDF
            if os.path.exists(sdf_file):
                os.remove(sdf_file)

    print(f"Prepared {count} PDBQT files in {output_dir}")
    return count


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python batch_prepare_ligands.py <input.smi> <output_dir> [max_molecules]")
        sys.exit(1)

    input_smi = sys.argv[1]
    output_dir = sys.argv[2]
    max_mols = int(sys.argv[3]) if len(sys.argv) > 3 else 100

    # Create a temporary SDF from SMILES
    tmp_sdf = "data/processed/batch_ligands.sdf"
    count = smiles_to_sdf(input_smi, tmp_sdf, max_mols)

    if count > 0:
        sdf_to_pdbqt_batch(tmp_sdf, output_dir)
