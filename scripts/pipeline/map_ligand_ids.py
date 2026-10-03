"""
Map ligand IDs (lig_XXXX.pdbqt) back to SMILES.

Usage:
    python map_ligand_ids.py <smiles_file> <output_csv> [max_molecules]
"""
import sys
import pandas as pd
from rdkit import Chem


def map_ids(smiles_file, output_csv, max_molecules=100):
    """Map sequential ligand IDs to their SMILES."""
    records = []
    count = 0

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
                continue

            records.append({
                "ligand_id": f"lig_{count:04d}",
                "smiles": smi,
            })
            count += 1

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"Mapped {len(df)} ligands to SMILES")
    print(f"Saved to {output_csv}")
    return df


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python map_ligand_ids.py <smiles_file> <output_csv> [max_molecules]")
        sys.exit(1)

    smiles_file = sys.argv[1]
    output_csv = sys.argv[2]
    max_mols = int(sys.argv[3]) if len(sys.argv) > 3 else 100

    map_ids(smiles_file, output_csv, max_mols)
